import os
import logging
import yaml

def missing_idx(file_name):
    logging.error(file_name + ' not provided!')
    raise NameError(file_name + ' not provided!')

def missing_yml(file_name):
    logging.error('· [ERROR] Unable to load yml file\n')
                  logging.error(exc)
    raise NameError(yaml.YAMLError)
    
# -------------------------------------------------------
# Loggin outputs
# -------------------------------------------------------

def reducing_file(file_name):
    logging.info('-'*50 + '\n    {}    \n'.format(file_name)
                 + '-'*50 + '\n')
