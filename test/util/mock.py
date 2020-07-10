from unittest.mock import create_autospec
import os
import sys
from shutil import rmtree, copytree
import logging
import json

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace

from kb_PICRUSt2.kb_PICRUSt2Impl import run_check
from kb_PICRUSt2.util.dprint import dprint
from kb_PICRUSt2.util.config import var
from .upa import *


##################################
##################################
testData_dir = '/kb/module/test/data'
##################################
##################################





def get_mock_dfu(dataset):
    '''
    Avoid lengthy `get_objects` and `save_objects``:w
    '''
    # validate
    if dataset not in ['17770', 'secret', 'dummy_10by8']:
        raise NotImplementedError('Input dataset not recognized')

    mock_dfu = create_autospec(DataFileUtil, instance=True)

    ##
    ## mock `save_objects`
    def mock_dfu_save_objects(params):
        params_str = str(params)
        if len(params_str) > 100: params_str = params_str[:100] + ' ...'
        logging.info('Mocking `dfu.save_objects` with `params=%s`' % params_str)

        return [['-1111', 1, 2, 3, '-1111', 5, '-1111']] # UPA made from pos 6/0/4
    
    mock_dfu.save_objects.side_effect = mock_dfu_save_objects

    ##
    ## mock `get_objects`
    def mock_dfu_get_objects(params):
        logging.info('Mocking `dfu.get_objects` with `params=%s`' % str(params))

        upa = params['object_refs'][0]
        flnm = {
            _17770: 'get_objects_AmpliconSet.json',
            _17770_AmpMat: 'get_objects_AmpliconMatrix.json',
            _17770_AttrMap: 'get_objects_AttributeMapping.json',
            secret: 'get_objects_AmpliconSet.json',
            secret_AmpMat: 'get_objects_AmpliconMatrix.json',
            dummy_10by8: 'get_objects_AmpliconSet.json',
            dummy_10by8_AmpMat: 'get_objects_AmpliconMatrix.json',
            dummy_10by8_AttrMap: 'get_objects_AttributeMapping.json',
            }[upa]
        flpth = os.path.join(testData_dir, 'by_dataset_input', dataset, 'get_objects', flnm)

        with open(flpth) as f:
            obj = json.load(f)

        return obj

    mock_dfu.get_objects.side_effect = mock_dfu_get_objects

    return mock_dfu


def get_mock_run_check(dataset):
    '''
    Avoid expensive runs of tool
    Copy over `var.out_dir`
    '''
    if dataset not in ['17770', 'secret', 'dummy_10by8']:
        raise NotImplementedError()

    mock_run_check = create_autospec(run_check)

    # side effect
    def mock_run_check_(cmd):
        logging.info('Mocking running cmd `%s`' % cmd)

        # test data
        src_flpth = os.path.join(testData_dir, 'by_dataset_input', dataset, 'return/PICRUSt2_output')

        # check if it already exists
        # since `check_run` is called twice in this app
        # and calling `copytree` twice would except
        if not os.path.exists(var.out_dir):
            copytree(src_flpth, var.out_dir)


    mock_run_check.side_effect = mock_run_check_

    return mock_run_check


def get_mock_kbr(dataset=None): 
    '''
    Avoid lengthy `create_extended_report`

    Does not use input currently
    '''

    mock_kbr = create_autospec(KBaseReport, instance=True) 

    # mock `create_extended_report`
    def mock_create_extended_report(params):
        logging.info('Mocking `kbr.create_extended_report`')

        return {
            'name': 'kbr_mock_name',
            'ref': 'kbr/mock/ref',
        }

    mock_kbr.create_extended_report.side_effect = mock_create_extended_report
    
    return mock_kbr



