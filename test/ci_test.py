# -*- coding: utf-8 -*-
import os
import time
import logging
import unittest
from unittest.mock import patch
from configparser import ConfigParser
import uuid
import pandas as pd
import numpy as np
import itertools
import shutil
import tracemalloc

from kb_PICRUSt2.kb_PICRUSt2Server import MethodContext
from kb_PICRUSt2.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2 
from kb_PICRUSt2.impl.kbase_obj import AmpliconMatrix, AttributeMapping
from kb_PICRUSt2.impl import appfile
from kb_PICRUSt2.impl.report import do_heatmap, HTMLReportWriter
from kb_PICRUSt2.impl.config import Var
from kb_PICRUSt2.util.debug import dprint
from mock import *
import config


######################################
######################################
######### TOGGLE PATCH ###############
######################################
###################################### 
do_patch = True

if do_patch:
    patch_ = patch
    patch_dict_ = patch.dict

else:
    patch_ = lambda *args, **kwargs: lambda f: f
    patch_dict_ = lambda *args, **kwargs: lambda f: f
######################################
######################################
######################################
######################################


class kb_PICRUSt2Test(config.BaseTest):
 
####################################################################################################
####################################################################################################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30'))  # ?
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma50by30'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma50by30'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_has_row_AttributeMapping_create_all(self):
        ret = config.get_serviceImpl().run_picrust2_pipeline(
            config.ctx, {
                **config.get_ws(),
                'amplicon_matrix_upa': enigma50by30,
                'output_name': 'an_output_name',
            }
        )

        self.assertTrue(len(Var.objects_created) == 8, Var.objects_created) 


####################################################################################################
####################################################################################################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *args: get_mock_dfu('enigma50by30_noAttrMaps_noSampleSet'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma50by30_noAttrMaps_noSampleSet'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma50by30_noAttrMaps_noSampleSet'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_has_no_row_AttributeMapping(self):
        '''
        Will not update/create row AttributeMapping if already none
        Default is to create all the FPs
        '''
        ret = config.get_serviceImpl().run_picrust2_pipeline(
            config.ctx, {
                **config.get_ws(),
                'amplicon_matrix_upa': enigma50by30_noAttrMaps_noSampleSet,
                'output_name': 'an_output_name',
            }
        )
        
        self.assertTrue(len(Var.objects_created) == 6, Var.objects_created)

####################################################################################################
####################################################################################################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30')) # ?
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma50by30'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma50by30'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_FP_options(self):
        '''
        Creating amplicon FPs logic is independent of creating sample FPs logic
        '''

        with self.subTest():
            ret = config.get_serviceImpl().run_picrust2_pipeline(
                config.ctx, {
                    **config.get_ws(),
                    'amplicon_matrix_upa': enigma50by30,
                    'fp_options': {
                        'create_amplicon_fps': False,
                        'create_sample_fps': False,
                    },
                    'output_name': 'an_output_name',
                }
            )

            self.assertTrue(len(Var.objects_created) == 2, Var.objects_created) 


        with self.subTest():
            ret = config.get_serviceImpl().run_picrust2_pipeline(
                config.ctx, {
                    **config.get_ws(),
                    'amplicon_matrix_upa': enigma50by30,
                    'fp_options': {
                        'create_amplicon_fps': True,
                        'create_sample_fps': True,
                    },
                    'output_name': 'an_output_name',
                }
            )

            self.assertTrue(len(Var.objects_created) == 8, Var.objects_created) 


####################################################################################################
####################################################################################################
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma17770by511'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma17770by511'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma17770by511'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_large_dataset(self):
        '''
        '''

        ret = config.get_serviceImpl().run_picrust2_pipeline(
            config.ctx, {
                **config.get_ws(),
                'amplicon_matrix_upa': enigma17770by511,
                'fp_options': {
                    'create_amplicon_fps': True,
                    'create_sample_fps': True,
                },
                'output_name': 'an_output_name',
            }
        )

        self.assertTrue(len(Var.objects_created) == 8, Var.objects_created) 
        

####################################################################################################
####################################################################################################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda u: get_mock_dfu('userTest'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('userTest'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('userTest'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda u: get_mock_kbr())
    def test_userTest_data(self):
        ret = config.get_serviceImpl().run_picrust2_pipeline(
            config.ctx, {
                **config.get_ws(),
                'amplicon_matrix_upa': userTest,
                'fp_options': {
                    'create_amplicon_fps': True,
                    'create_sample_fps': True,
                },
                'output_name': 'an_output_name',
            })
        self.assertTrue(len(Var.objects_created) == 6, Var.objects_created) 

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    @classmethod
    def list_tests(cls):
        return [key for key, value in cls.__dict__.items() if type(key) == str and key.startswith('test') and callable(value)]

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        dec = '!!!' * 300
        print(dec, "DON'T FORGET TO SEE HTML(S)", dec)
        skipped_tests = list(set(all_tests) - set(cls.list_tests()))
        print('* All tests (%d): %s' % (len(all_tests), all_tests))
        print('* Tests skipped (%d): %s' % (len(skipped_tests), skipped_tests))
        print('* Tests run (%d): %s' % (len(cls.list_tests()), cls.list_tests()))


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!! select what to run !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
all_tests = []
for key, value in kb_PICRUSt2Test.__dict__.copy().items():
    if key.startswith('test') and callable(value):
        all_tests.append(key)

integration_tests = [
    'test_has_row_AttributeMapping_create_all', 'test_has_no_row_AttributeMapping',
    'test_FP_options'
    'test_large_dataset',
]
large_tests = [
    'test_large_dataset', 'test_large_heatmap',
]
run_tests = [
    'test_userTest_data',
]

for test in all_tests:
        if test not in run_tests:
            delattr(kb_PICRUSt2Test, test)
            pass





