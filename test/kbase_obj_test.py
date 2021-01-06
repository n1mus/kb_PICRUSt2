import unittest
from unittest.mock import patch
import uuid
import pandas as pd
import numpy as np
import itertools
import shutil
from pytest import raises

from kb_PICRUSt2.kb_PICRUSt2Server import MethodContext
from kb_PICRUSt2.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2 
from kb_PICRUSt2.impl.kbase_obj import AmpliconMatrix, AttributeMapping
from kb_PICRUSt2.impl import appfile
from kb_PICRUSt2.impl.report import do_heatmap, HTMLReportWriter
from kb_PICRUSt2.impl.config import Var
from kb_PICRUSt2.impl.error import * # Exceptions
from kb_PICRUSt2.util.debug import dprint
from kb_PICRUSt2.util.validate import ValidationException
from mock import *





####################################################################################################
####################################################################################################
@patch.dict('kb_PICRUSt2.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('dummy_10by8')})
def test_AmpliconMatrix_validation():
    '''
    Test validation of amplicon table in AmpliconMatrix
    Should be (1) count data,  missing (None) allowed
    Can assume that obj holds data in list of lists of numeric/None
    '''
    
    logging.info('Testing with test_AmpliconMatrix_validation')

    amp_mat = AmpliconMatrix(dummy_10by8_AmpMat) # these values have been truncated to ints
    amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [0.0, 0.0, 1319.0, 1.0] # float
    amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [0, 0, 1319, 1] # int
    amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, 0., 0., 1319., 1.] # float, with missing
    amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, 0, 0, 1319, 1] # int, with missing
    amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, 0, -0., 1319.0, 1] # int/float, with missing
    amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, None, 0, -0, 0.0, 0, 0.] # 0s, with missing
    amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, 0.999999999] # close enough
    amp_mat.validate_amplicon_abundance_data() 

    amp_mat.obj['data']['values'] = [None, -0.0000000001] # close enough
    amp_mat.validate_amplicon_abundance_data() 

    amp_mat.obj['data']['values'] = [0.9]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data() 

    amp_mat.obj['data']['values'] = [-1]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data() 

    amp_mat.obj['data']['values'] = [None, None, None]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data() 

    amp_mat.obj['data']['values'] = [None]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, -1]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, None, 1.00001]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [-1.0, 0, 1319]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, 0, 1, 2, 3, 4.5]
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data()

    amp_mat.obj['data']['values'] = [None, 0.0, 1.0, 2.0, 3.0, 4.00001] # 4.00001 would pass with np.allclose default rtol
    with raises(ValidationException): amp_mat.validate_amplicon_abundance_data()




####################################################################################################
####################################################################################################
@patch.dict('kb_PICRUSt2.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('dummy_10by8'), })
def test_AttributeMapping():
    '''
    Mostly writing attributes
    '''
    # TODO test with non-identical row_mapping

    logging.info('Testing with test_AttributeMapping')

    amp_mat = AmpliconMatrix(dummy_10by8_AmpMat)
    attr_map = AttributeMapping(dummy_10by8_AttrMap, amp_mat)

    ##
    ## write new attribute/source
    ind_0, name_0 = attr_map.add_attribute_slot('cloud type', 'testing')
    assert ind_0 == 2 
    assert name_0 == 'cloud type'
    assert len(attr_map.obj['attributes']) == len(list(attr_map.obj['instances'].values())[0]) 

    attr_map.map_update_attribute(ind_0, {
        "amplicon_id_0": "dummy0",
        "amplicon_id_1": "dummy0",
        "amplicon_id_2": "dummy0",
        "amplicon_id_3": "dummy0",
        "amplicon_id_4": "dummy0",
        "amplicon_id_5": "dummy0",
        "amplicon_id_6": "dummy0",
        "amplicon_id_7": "dummy0",
        "amplicon_id_8": "dummy0",
        "amplicon_id_9": "dummy0"
    })

    assert(attr_map.obj['instances']['amplicon_id_4'][ind_0] == 'dummy0')

    ##
    ## overwrite attribute/source
    ind_1, name_1 = attr_map.add_attribute_slot('celestial body', 'testing')
    assert ind_1 == 3, ind_1
    assert name_1 == 'celestial body (1)', name_1
    assert len(attr_map.obj['attributes']) == len(list(attr_map.obj['instances'].values())[0]) 

    attr_map.map_update_attribute(ind_1, {
        "amplicon_id_0": "dummy1",
        "amplicon_id_1": "dummy1",
        "amplicon_id_2": "dummy1",
        "amplicon_id_3": "dummy1",
        "amplicon_id_4": "dummy1",
        "amplicon_id_5": "dummy1",
        "amplicon_id_6": "dummy1",
        "amplicon_id_7": "dummy1",
        "amplicon_id_8": "dummy1",
        "amplicon_id_9": "dummy1"
    })

    ##
    ## all same length
    num_attr = len(attr_map.obj['attributes'])
    for attr_l in attr_map.obj['instances'].values():
        assert(len(attr_l) == num_attr)

    ## 
    ## check did not add dummy attribute to wrong slot
    ind_lftvr = list(set(range(num_attr)) - {ind_0, ind_1})

    for attr_l in attr_map.obj['instances']:
        for ind in ind_lftvr:
            assert('dummy' not in attr_l[ind])


