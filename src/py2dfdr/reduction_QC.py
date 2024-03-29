import os
import numpy as np
from astropy.io import fits
from astropy.visualization import simple_norm
from matplotlib import pyplot as plt
import logging

flux_cmap = plt.get_cmap('magma').copy()
flux_cmap.set_under('fuchsia')
flux_cmap.set_over('red')

def check_exists(path):
    """Check if a file exists.
    Params
    ------
    path: (str) Path to file or directory
    
    Returns
    -------
    exists: (bool)
    """
    return os.path.exists(path)

def check_exists_decorator(func):
    """Check if file exists before calling function."""
    def wrapper(*args, **kwargs):
        if check_exists(args[0]):
            return func(*args, **kwargs)
        else:
            logging.warning('[QC] WARNING: Input path {args[0]} does not exist')
            return
    return wrapper

@check_exists_decorator
def check_image(path, percentiles=None, save_dir=None, title=None):
    """blah."""
    logging.info('[QC] · QC plot for:\n   {}'.format(path))
    if percentiles is None:
        percentiles = [1, 5, 16, 50, 84, 95, 99]
    master = fits.getdata(path)
    percents = np.nanpercentile(master.flatten(), percentiles)
    if not np.isfinite(percents).any():
        logging.warning('[QC] WARNING: ALL PIXELS HAVE NAN COUNTS, EXITING PLOT')
        return

    fig = plt.figure(figsize=(10, 10))
    if title is not None:
        if title == 'auto':
            title = fits.getval(path, 'OBJECT')
        fig.suptitle(title)
    else:
        fig.suptitle(path)
    ####################################
    ax = fig.add_subplot(221)
    ax.set_title("Counts histogram")
    ax.hist(master.flatten(), range=[percents[0], percents[-1]],
            bins=master.flatten().size//1000, log=True, color='k')
    ax.set_xlabel('counts')
    ax.set_ylabel('# pixels')
    for pcnt, label in zip(percents, percentiles):
        ax.axvline(pcnt, ls='-', label=label)
    ax.legend(loc='center', bbox_to_anchor=(0.5, 1.2), ncol=3)
    ####################################
    ax = fig.add_subplot(222)
    ax.set_title("2D reduced file")
    mappable = ax.imshow(master, vmin=percents[0], vmax=percents[-1],
                            cmap=flux_cmap, aspect='auto',
                            origin='lower')
    plt.colorbar(mappable, label='counts', extend='both')
    ####################################
    central_fibre = master.shape[0] // 2
    ax = fig.add_subplot(223)
    ax.set_title("Fiber {} spectra".format(central_fibre))
    ax.plot(master[central_fibre], c='k')
    ax.set_ylim(percents[0], percents[-1])
    ####################################
    central_column = master.shape[1] // 2
    ax = fig.add_subplot(224)
    ax.set_title("Column {} spectra".format(central_column))
    ax.plot(master[:, central_column], c='k')
    ax.set_ylim(percents[0], percents[-1])

    if save_dir is not None:
        fig.savefig(save_dir, bbox_inches='tight')
        logging.info('[QC] Plot saved as:\n {}'.format(save_dir))
    plt.clf()
    plt.close()


@check_exists_decorator
def check_saturated(path, sat_level=65500, log=True, plot=True):
    """Compute the fraction of saturated (nan/inf) pixels."""

    data = fits.getdata(path)
    finite_values = np.isfinite(data) & (data < sat_level)
    frac_sat = finite_values[~finite_values].size / finite_values.size
    
    pixel_value = np.sort(data[np.isfinite(data)].flatten())
    cumulative_distrib = np.arange(1, pixel_value.size + 1, 1)
    cumulative_fraction = cumulative_distrib / cumulative_distrib[-1]

    if log:
        logging.info('[QC] · Checking saturation levels for\n   {}'.format(path))
        logging.info('[QC] · Fraction of saturated pixels={:.3f}'.format(frac_sat))
    if plot:
        title = fits.getval(path, 'OBJECT')
        fig = plt.figure(constrained_layout=True)
        plt.suptitle(title)
        gs = fig.add_gridspec(2, 5)
        ax = fig.add_subplot(gs[:, :4])
        norm = simple_norm(data, 'sqrt', max_cut=sat_level)
        mappable = ax.imshow(data, origin='lower', norm=norm, cmap=flux_cmap,
        interpolation='none')
        plt.colorbar(mappable, aspect=40)
        ax.set_xlabel('Column pixel')
        ax.set_ylabel('Row pixel')

        ax = fig.add_subplot(gs[0, -1])
        ax.plot(pixel_value, cumulative_fraction)
        ax.set_xlim(1, sat_level)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Pixel value')
        ax.set_ylabel('Cum. fraction')
        ax.grid(visible=True)

        ax = fig.add_subplot(gs[1, -1])
        ax.plot(pixel_value, cumulative_fraction)
        p10_idx = np.searchsorted(cumulative_fraction, .1)
        p90_idx = np.searchsorted(cumulative_fraction, .9)
        ax.set_xlim(pixel_value[p10_idx], pixel_value[p90_idx])
        ax.set_ylim(0.1, 0.9)
        ax.set_xlabel('Pixel value')
        ax.set_ylabel('Cum. fraction')
        ax.grid(visible=True)

        fig_path = path.replace(".fits", "_qc_sat.png")
        fig.savefig(fig_path, bbox_inches='tight', dpi=200)


    return frac_sat


@check_exists_decorator
def get_keyword(path, keyword, hdu_index=0):
    """Return the keyword value of a fits file."""
    with fits.open(path) as f:
        try:
            val = f[hdu_index].header[keyword]
        except Exception:
            val = None
    return val


@check_exists_decorator
def check_tramline(path, plot=True):
    """blah."""
    logging.info('[QC] · Test for Tramline:\n   {}'.format(path))
    with fits.open(path) as f:
        median_fwhm = f[0].header['MWIDTH']
        fibrepos = f[0].data
        mean_amplitude = f[1].data
    # Dispersion in pixels along the spectral axis
    fibrepos_disp = np.std(fibrepos, axis=1)
    intensity_disp = np.std(mean_amplitude, axis=0)
    # Quality assesment
    fibre_separation = (median_fwhm > 2) & (median_fwhm < 3)
    smoothness = ((fibrepos_disp.max() < 10) & (fibrepos_disp.min() > 0)
                  & (intensity_disp.max() < 10) & (intensity_disp.min() > 0))
    if fibre_separation & smoothness:
        bad_tramline = False
    else:
        bad_tramline = True
    # Quality control plots
    if plot:
        output = path.replace('.fits', '.png')
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(6, 6),
                                gridspec_kw=dict(wspace=0.5, hspace=0.5))
        fig.suptitle(path, fontsize=8)
        ax = axs[0, 0]
        ax.set_title(r'$\mu(\lambda)$')
        mappable = ax.imshow(fibrepos, aspect='auto', origin='lower')
        plt.colorbar(mappable, ax=ax)
        ax.set_ylabel('Fibre')
        ax.set_xlabel(r'$\lambda (pix)$')
        ax = axs[0, 1]
        ax.set_title('Median FWHM\nbetween fibres: {:.2f} pix'.format(
            median_fwhm))
        ax.plot(fibrepos_disp)
        ax.set_ylabel(r'$\sigma(\mu)$')
        ax.set_xlabel('Fibre')
        ax = axs[1, 0]
        ax.set_title(r'$I(\lambda)$')
        mappable = ax.imshow(mean_amplitude, aspect='auto', origin='lower')
        plt.colorbar(mappable, ax=ax)
        ax.set_ylabel('Fibre')
        ax.set_xlabel(r'$\lambda (pix)$')
        ax = axs[1, 1]
        ax.plot(intensity_disp)
        ax.set_ylabel(r'$\sigma(I)$')
        ax.set_xlabel(r'$\lambda (pix)$')
        fig.savefig(output, bbox_inches='tight')
        plt.clf()
        plt.close(fig)
        logging.info('[QC] ·  Plot saved as:\n {}'.format(output))
    logging.info('[QC] ·  Bad tramline: {}'.format(bad_tramline))
    return bad_tramline


@check_exists_decorator
def clean_nan(path):
    """blah..."""
    with fits.open(path, mode='update') as hdul:
        nans_mask = np.isfinite(hdul[0].data)
        fits_shape = nans_mask.shape
        nans_pos = np.where(nans_mask == False)
        if nans_pos[0].size > 0:
            for entry in zip(nans_pos[0], nans_pos[1]):
                hdul[0].data[entry[0], entry[1]] = np.nanmean([
                    hdul[0].data[np.clip(entry[0]-1, a_min=0, a_max=None),
                                 entry[1]],
                    hdul[0].data[entry[0],
                                 np.clip(entry[1]-1, a_min=0, a_max=None)],
                    hdul[0].data[np.clip(entry[0]+1, a_min=None,
                                         a_max=fits_shape[0]-1),
                                 entry[1]],
                    hdul[0].data[entry[0],
                                 np.clip(entry[1]+1, a_min=None,
                                         a_max=fits_shape[1]-1)]])
            hdul.flush()
        hdul.close()

# Mr Krtxo \(ﾟ▽ﾟ)/

