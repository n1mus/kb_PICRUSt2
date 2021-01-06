import time
import logging
import subprocess
import sys
import gzip
import shutil


class NonZeroReturnException(Exception): pass



####################################################################################################
####################################################################################################
def run_check(cmd: str, shell=True):
    logging.info('Running cmd `%s`' % cmd)
    t0 = time.time() 
    
    completed_proc = subprocess.run(cmd, shell=shell, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stderr)

    logging.info('Cmd took %.2fmin' % ((time.time() - t0)/60))

    if completed_proc.returncode != 0:
        raise NonZeroReturnException(
            "Command `%s` exited with non-zero return code `%d`. "
            "Check logs for more details" %
            (cmd, completed_proc.returncode)
        )



####################################################################################################
####################################################################################################
def gunzip(read, write):
    '''
    Gunzip from `read` to `write`
    '''
    with gzip.open(read, 'rb') as fh_read:
        with open(write, 'wb') as fh_write:
            shutil.copyfileobj(fh_read, fh_write)

