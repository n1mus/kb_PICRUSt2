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

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2 
from kb_PICRUSt2.impl.kbase_obj import AmpliconMatrix, AttributeMapping
from kb_PICRUSt2.impl import appfile
from kb_PICRUSt2.impl.config import Var
from kb_PICRUSt2.impl.error import * # Exceptions
from kb_PICRUSt2.util.debug import dprint
from mock import *
import config


run_dir = os.path.join(config.shared_folder, 'test_appfile_' + str(uuid.uuid4()))
os.mkdir(run_dir)


####################################################################################################
####################################################################################################
#@patch.dict('kb_PICRUSt2.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('dummy_10by8')})
def test_parse_picrust2_traits():
    '''
    '''
    logging.info('Testing with test_appfile')


    ## Test `appfile.parse_picrust2_traits` ##
    
    flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/PICRUSt2_output/pathways_out/path_abun_predictions.tsv.gz'

    id2traits_d = appfile.parse_picrust2_traits(flpth)

    ans = {
       "amplicon_id_0": "N10-formyl-tetrahydrofolate biosynthesis,4-hydroxyphenylacetate degradation,aerobactin biosynthesis,superpathway of chorismate metabolism,homolactic fermentation",
       "amplicon_id_1":                                          "4-hydroxyphenylacetate degradation,aerobactin biosynthesis,superpathway of chorismate metabolism,homolactic fermentation",
       "amplicon_id_2":                                                                             "aerobactin biosynthesis,superpathway of chorismate metabolism,homolactic fermentation",
       "amplicon_id_3":                                                                                                     "superpathway of chorismate metabolism,homolactic fermentation",
       "amplicon_id_4":                                                                                                                                           "homolactic fermentation",
       "amplicon_id_5":                                                                                                                                                                  "",
       "amplicon_id_6": "N10-formyl-tetrahydrofolate biosynthesis,4-hydroxyphenylacetate degradation,aerobactin biosynthesis,superpathway of chorismate metabolism,homolactic fermentation",
       "amplicon_id_7": "N10-formyl-tetrahydrofolate biosynthesis,4-hydroxyphenylacetate degradation,aerobactin biosynthesis,superpathway of chorismate metabolism",
       "amplicon_id_8": "N10-formyl-tetrahydrofolate biosynthesis,4-hydroxyphenylacetate degradation,aerobactin biosynthesis",
       "amplicon_id_9": "N10-formyl-tetrahydrofolate biosynthesis,4-hydroxyphenylacetate degradation"
    }

    assert id2traits_d == ans


