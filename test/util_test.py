import unittest
from unittest.mock import patch
from configparser import ConfigParser
import uuid
import pandas as pd
import numpy as np
import itertools
import shutil

from kb_PICRUSt2.kb_PICRUSt2Server import MethodContext
from kb_PICRUSt2.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_PICRUSt2.util.debug import dprint
from kb_PICRUSt2.util.cli import run_check, NonZeroReturnException
from kb_PICRUSt2.util.file import get_numbered_duplicate
from mock import *
import config

class kb_PICRUSt2Test(config.BaseTest):

####################################################################################################
####################################################################################################
    def test_run_check(self):
        with self.assertRaisesRegex(NonZeroReturnException, '`127`'):
            run_check('set -o pipefail && tmp |& tee tmp')

        with self.assertRaisesRegex(NonZeroReturnException, '`127`'):
            run_check('set -o pipefail && echo hi |& tmp')

        run_check('set -o pipefail && echo hi |& tee tmp') # run correctly


####################################################################################################
####################################################################################################
    def test_get_numbered_duplicate(self):
        # test numbering system
        q = 'the_attr'
        names = ['null']
        assert get_numbered_duplicate(names, q) == 'the_attr'
        names = ['the_attr']
        assert get_numbered_duplicate(names, q) == 'the_attr (1)'
        names = ['the_attr', 'the_attr (1)', 'the_attr (2)']
        assert get_numbered_duplicate(names, q) == 'the_attr (3)', get_numbered_duplicate(names, q)
        names = ['the_attr (1)', 'the_attr (2)']
        assert get_numbered_duplicate(names, q) == 'the_attr'
        names = ['the_attr', 'the_attr (1)', 'the_attr (5)']
        assert get_numbered_duplicate(names, q) == 'the_attr (2)'
        names = ['the_attr (0)', 'the_attr (1)', 'the_attr (2)']
        assert get_numbered_duplicate(names, q) == 'the_attr'
        names = ['the_attr (-1)', 'the_attr (1)', 'the_attr (2)']
        assert get_numbered_duplicate(names, q) == 'the_attr'

        # test internal regexing
        q = 'the[0-9]_attr'
        names = [q]
        assert get_numbered_duplicate(names, q) == q + ' (1)'


