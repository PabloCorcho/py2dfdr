import os
import logging
import yaml


def missing_idx(file_name):
    logging.error(file_name + ' not provided!')
    raise NameError(file_name + ' not provided!')


def missing_yml(file_name):
    logging.error('· [ERROR] Unable to load yml file')
    logging.error(yaml.YAMLError)
    raise NameError(yaml.YAMLError)


def missing_master(file_name):
    logging.error(
        '--> [ERROR] MASTER FILE *NOT* found at:\n  {}'.format(
            file_name))
    raise NoMasterFileError(file_name)


# -------------------------------------------------------
# Logging outputs
# -------------------------------------------------------
def log_header(file_name):
    """Make a pretty header line for a process."""
    logging.info('-'*50 + '\n    {}    \n'.format(file_name)
                 + '-'*50)


# -------------------------------------------------------
# Exceptions
# -------------------------------------------------------
class NoMasterDict(Exception):
    """Exception raised when a master file is not in memory"""
    def __init__(self, master):
        self.message = 'master-{} not in memory'.format(master)
        super().__init__(self.message)


class NoFileError(Exception):
    """Exception raised when a file is not found."""
    def __init__(self, file):
        self.message = '{} file not found'.format(file)
        super().__init__(self.message)


class NoIDXfileError(NoFileError):
    """Exception raised when a IDX file is not found."""
    def __init__(self, file):
        self.message = '{} IDX '.format(file)
        super().__init__(self.message)


class NoMasterFileError(NoFileError):
    """Exception raised when a Master file is not found."""
    def __init__(self, file):
        self.message = '{} master'.format(file)
        super().__init__(self.message)

# Mr Krtxo \(ﾟ▽ﾟ)/
