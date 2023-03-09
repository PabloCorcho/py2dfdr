"""
2dfdr wrapper module
Created by Pablo Corcho-Caballero
Last modification: 09/03/2023
"""
import os
import subprocess
import logging


def aaorun_cleanup(log=False):
    """Clean temporary files that may impact on 2dfdr performance
    
    Description
    -----------
    This method removes scratch files located at
     - $home/imp_scratch
     - /tmp
    """
    process = subprocess.run('cleanup', shell=True, timeout=60, stdout=subprocess.PIPE, text=True)
    if log:
        logging.info('[aaorun] Cleaning')


def aaorun_command(command, file, options=None, output=None,
                   idx_file='koala.idx',
                   aaorun='aaorun', wdir=None, timeout=900, log=False):
    """Run an aaorun command.

    Description
    -----------
    This method executes any aaorun command.

    Link to 2dfdr repo:
    https://dev.aao.org.au/rds/2dfdr

    To see the commands available for aaorun, type
    "aaorun help" on the terminal.
    For examples, type 
    "aaorun examples"

    Params
    ------
    - command: (str) 2dfdr command to be executed.
    - file: (str) abs_path to target file.
    - options: (str, optional, default=None) string containing the additional
        parameters to be passed.
    - output: (str, optional, default=None) file name where the output will
        be stored.
    - idx_file: (str) path to the idx file used for reducing data.
        Default is koala.idx
    - aaorun: (str, optional, default="aaorun") path to the binary executable file.
    - wdir: (str, optional, default=None) working directory. If None,
        the parent directory of "file" will be used as working dir.
    - timeout: (int, optional, default=900) Timeout in seconds.
    - log: (bool, optional, default=False) if True, a log file will be provided
        containing the 2dfdr output.

    Example
    -------
    # TODO
    """
    # Initialise
    aaorun_print('[aaorun] · Initialising AAORUN process: ' + command, log, level='INFO')
    # Set working directory
    if wdir is None:
        wdir = os.path.dirname(file)
    aaorun_print('[aaorun] · Working directory: ' + wdir, log, level='INFO')
    # Set additional options
    if options is None:
        cmd_options = ' '.join(['-idxfile %s' % idx_file, '-wdir %s' % wdir])
    else:
        extra_options = [opt for opt in options]
        extra_options.append('-idxfile %s' % idx_file)
        extra_options.append('-wdir %s' % wdir)
        cmd_options = ' '.join(extra_options)
    # Combine all command arguments
    aaorun_cmd = ' '.join([aaorun, command, file, cmd_options])
    # Run the command
    aaorun_print('[aaorun] · Running command: ' + aaorun_cmd, log, level='INFO')
    if output is None:
        try:
            process = subprocess.run(aaorun_cmd, shell=True,
                                     timeout=timeout, stdout=subprocess.PIPE,
                                     text=True)
            aaorun_cleanup(log=log)
            if process.returncode != 0:
                aaorun_print('[aaorun] · WARNING: Command \n {} \n FAILED! \n {}'.
                                    format(aaorun_cmd, process.stderr), log, level='warning')
                return 1
            else:
                return 0
        except subprocess.TimeoutExpired:
            aaorun_print('[aaorun] · WARNING: Command \n {} \n  Ran too long (>{:.1f})'
                         .format(aaorun_cmd, timeout),
                         log, level='warning')
            aaorun_cleanup(log=log)
            return 1
    else:
        with open(output, 'w') as outfile:
            try:
                process = subprocess.run(aaorun_cmd, shell=True,
                                         timeout=timeout, stdout=outfile,
                                         text=True)
                aaorun_cleanup(log=log)
                if process.returncode != 0:
                    aaorun_print('[aaorun] · WARNING: Command \n {} \n FAILED! \n {}'.
                                    format(aaorun_cmd, process.stderr), log, level='warning')
                    return 1
                else:
                    return 0
            except subprocess.TimeoutExpired:
                aaorun_print('[aaorun] · WARNING: Command \n {} \n  Ran too long (>{:.1f})'
                         .format(aaorun_cmd, timeout),
                         log, level='warning')
                aaorun_cleanup(log=log)
                return 1


def aaorun_print(msg, log, level='INFO'):
    """Print a message on the logger."""
    if log:
        if level == 'INFO':
            logging.info(msg)
        elif level == 'warning':
            logging.warning(msg)
        elif level == 'error':
            logging.error(msg)
    else:
        print("[aaorun] " + msg)

# Mr Krtxo
