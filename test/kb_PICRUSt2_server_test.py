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

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2, run_check, OutfileWrangler
from kb_PICRUSt2.util.config import var
from kb_PICRUSt2.util.dprint import dprint
from kb_PICRUSt2.util.kbase_obj import AmpliconMatrix, AttributeMapping, Report
from kb_PICRUSt2.util.report import do_heatmap, HTMLReportWriter
from kb_PICRUSt2.util.error import *
from util.mock import *
from util.upa import *

#tracemalloc.start(25)

######################################
######################################
######### TOGGLE PATCH ###############
######################################
###################################### 
do_patch = False 

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



class kb_PICRUSt2Test(unittest.TestCase):

####################################################################################################
####################################################################################################
############################# UNIT TESTS ###########################################################
####################################################################################################
####################################################################################################

    ####################
    ####################
    def test_run_check(self):
        '''
        Test `run_check` which runs PICRUSt2 executable
        ''' # TODO test return value of commands chained with &&
        logging.info('Testing with test_run_check')

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && ;s |& tee tmp')
            self.assertTrue('`2`') in str(cm.exception) # return code 2

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && tmp |& tee tmp')
            self.assertTrue('`127`') in str(cm.exception) # return code 127

        with self.assertRaises(NonZeroReturnException) as cm:
            run_check('set -o pipefail && echo hi |& tmp')
            self.assertTrue('`127`') in str(cm.exception) # return code 127

        run_check('set -o pipefail && echo hi |& tee tmp') # run correctly


    ####################
    ####################
    @patch.dict('kb_PICRUSt2.util.kbase_obj.var', values={'dfu': get_mock_dfu('dummy_10by8')})
    def test_OutfileWrangler(self):
        '''
        '''
        logging.info('Testing with test_OutfileWrangler')

        run_dir = os.path.join('/kb/module/work/tmp', 'test_OutfileWranger_' + str(uuid.uuid4()))
        os.mkdir(run_dir)

        ## Test `OutfileWrangler.parse_picrust2_traits` ##
        with self.subTest('Test OutfileWrangler.parse_picrust2_traits'):
        
            flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/PICRUSt2_output/pathways_out/path_abun_predictions.tsv.gz'

            id2traits_d = OutfileWrangler.parse_picrust2_traits(flpth)

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

            self.assertTrue(id2traits_d == ans)

        ## Test `OutfileWrangler.pad_0_vecs` ##
        with self.subTest('Test OutfileWrangler.pad_0_vecs'):

            amp_mat = AmpliconMatrix(dummy_10by8_AmpMat)

            # amplicon x func
            flpth0 = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/PICRUSt2_output/pathways_out/path_abun_predictions.tsv'
            flpth1 = os.path.join(run_dir, os.path.basename(flpth0))

            df0 = pd.read_csv(flpth0, sep='\t', index_col=0, header=0)
            df_drop = df0.drop(index=df0.index[np.where(df0.sum(axis=1) == 0)][0])
            df_drop.to_csv(flpth1, sep='\t')

            OutfileWrangler.pad_0_vecs(flpth1, amp_mat)
            df1 = pd.read_csv(flpth1, sep='\t', index_col=0, header=0)

            assert np.allclose(
                df0.values,
                df1.values
            )
            
            # func x sample    
            flpth1 = os.path.join(run_dir, os.path.basename(flpth0))

            df0 = pd.read_csv(flpth0, sep='\t', index_col=0, header=0)
            df_drop = df0.drop(index=df0.index[np.where(df0.sum(axis=1) == 0)][0])
            df_drop.to_csv(flpth1, sep='\t')

            OutfileWrangler.pad_0_vecs(flpth1, amp_mat)
            df1 = pd.read_csv(flpth1, sep='\t', index_col=0, header=0)

            assert np.allclose(
                df0.values,
                df1.values
            )


    ####################
    ####################
    @patch.dict('kb_PICRUSt2.util.kbase_obj.var', values={'dfu': get_mock_dfu('dummy_10by8')})
    def test_AmpliconSet_and_AmpliconMatrix(self):
        '''
        Combine AmpliconSet and AmpliconMatrix since they mostly write input files
        '''
        logging.info('Testing with test_AmpliconSet_and_AmpliconMatrix')

        # set up `run_dir`
        var.run_dir = os.path.join(self.shared_folder, 'test_AmpliconSet_and_AmpliconMatrix_' + str(uuid.uuid4()))
        os.mkdir(var.run_dir)

        #amp_set = AmpliconSet(dummy_10by8)
        amp_mat = AmpliconMatrix(dummy_10by8_AmpMat)

        # write

        seq_flpth = os.path.join(var.run_dir, 'study_seqs.fna')
        seq_abundance_table_flpth = os.path.join(var.run_dir, 'study_seqs.tsv')

        #amp_set.to_fasta(seq_flpth)
        amp_mat.to_seq_abundance_table(seq_abundance_table_flpth)

        # compare
        
        seq_ref_flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.fna'
        seq_abundance_table_ref_flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.tsv'
        
        #with open(seq_flpth) as f1:
        #    with open(seq_ref_flpth) as f2:
        #        self.assertTrue(f1.read() == f2.read())

        with open(seq_abundance_table_flpth) as f1:
            with open(seq_abundance_table_ref_flpth) as f2:
                self.assertTrue(f1.read() == f2.read())



    ####################
    ####################
    @patch.dict('kb_PICRUSt2.util.kbase_obj.var', values={'dfu': get_mock_dfu('dummy_10by8'), 'warnings': []})
    def test_AttributeMapping(self):
        '''
        Mostly writing attributes
        '''
        logging.info('Testing with test_AttributeMapping')

        amp_mat = AmpliconMatrix(dummy_10by8_AmpMat)
        attr_map = AttributeMapping(dummy_10by8_AttrMap, amp_mat)

        ##
        ## write new attribute/source
        ind_0 = attr_map.get_attribute_slot_warn('biome', 'testing')
        self.assertTrue(ind_0 == 2)
        self.assertTrue(len(var.warnings) == 0)

        attr_map.update_attribute(ind_0, {
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

        self.assertTrue(attr_map.obj['instances']['amplicon_id_4'][ind_0] == 'dummy0')

        ##
        ## overwrite attribute/source
        ind_1 = attr_map.get_attribute_slot_warn('celestial body', 'upload')
        self.assertTrue(ind_1 == 0)
        self.assertTrue(len(var.warnings) == 1)

        attr_map.update_attribute(ind_1, {
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
            self.assertTrue(len(attr_l) == num_attr)

        ## 
        ## check did not add dummy attribute to wrong slot
        ind_lftvr = list(set(range(num_attr)) - {ind_0, ind_1})

        for attr_l in attr_map.obj['instances']:
            for ind in ind_lftvr:
                self.assertTrue('dummy' not in attr_l[ind])


    ####################
    ####################
    @patch.dict('kb_PICRUSt2.util.report.var', values={'warnings': []})
    def test_large_heatmap(self):
        '''
        Test largest possible heatmap
        '''
        logging.info('Testing with test_large_heatmap')

        ##
        def write_random_tsv(flpth, dim, max):
            values = (np.random.random((dim, dim)) * max).round(decimals=2)
            df = pd.DataFrame(
                    values, 
                    index=['dummy_ind_%d' % i for i in range(dim)], 
                    columns=['dummy_col_%d' % i for i in range(dim)]
            )
            df.to_csv(flpth, sep='\t', compression='gzip')

        ##
        def has_n_htmls(dir, n):
            return len([flnm for flnm in os.listdir(report_dir) if flnm.endswith('.html')]) == n
    
        ## make `run_dir`
        run_dir = os.path.join(self.shared_folder, 'test_report_' + str(uuid.uuid4()))
        os.mkdir(run_dir)

        ###
        ### heatmap random large
        
        report_dir = os.path.join(run_dir, 'report_randomLargeHeatmap')
        fig_dir = os.path.join(report_dir, 'fig')
        os.makedirs(fig_dir)

        tsvgz_flpth = os.path.join(report_dir, 'random.tsv.gz')
        png_flpth = os.path.join(fig_dir, 'heatmap_random.png')
        html_flpth = os.path.join(report_dir, 'heatmap_random.html')

        write_random_tsv(tsvgz_flpth, dim=7000, max=1500)

        do_heatmap(tsvgz_flpth, png_flpth, html_flpth)


    ####################
    ####################
    @patch.dict('kb_PICRUSt2.util.report.var', values={'warnings': []})
    def test_small_heatmap(self):
        '''
        '''
        logging.info('Testing with test_small_heatmap')

        ##
        def write_random_tsv(flpth, dim, max):
            values = (np.random.random((dim, dim)) * max).round(decimals=2)
            df = pd.DataFrame(
                    values, 
                    index=['dummy_ind_%d' % i for i in range(dim)], 
                    columns=['dummy_col_%d' % i for i in range(dim)]
            )
            df.to_csv(flpth, sep='\t', compression='gzip')

        ##
        def has_n_htmls(dir, n):
            return len([flnm for flnm in os.listdir(report_dir) if flnm.endswith('.html')]) == n
    
        ## make `run_dir`
        run_dir = os.path.join(self.shared_folder, 'test_report_' + str(uuid.uuid4()))
        os.mkdir(run_dir)

        ###
        ### heatmap random small
        
        report_dir = os.path.join(run_dir, 'report_randomSmallHeatmap')
        fig_dir = os.path.join(report_dir, 'fig')
        os.makedirs(fig_dir)

        tsvgz_flpth = os.path.join(report_dir, 'random.tsv.gz')
        png_flpth = os.path.join(fig_dir, 'heatmap_random.png')
        html_flpth = os.path.join(report_dir, 'heatmap_random.html')

        write_random_tsv(tsvgz_flpth, dim=500, max=1500)

        do_heatmap(tsvgz_flpth, png_flpth, html_flpth)




    ####################
    ####################
    @patch.dict('kb_PICRUSt2.util.report.var', values={'warnings': []})
    def test_report(self):
        '''
        Should make 6 `report_dir_*` subdirectories, 
        Check them (`cd test_local/workdir/tmp && firefox test_report_*/report_dir_*/*.html &)
        '''

        logging.info('Testing with test_report')


        ##
        def write_random_tsv(flpth, dim, max):
            values = (np.random.random((dim, dim)) * max).round(decimals=2)
            df = pd.DataFrame(
                    values, 
                    index=['dummy_ind_%d' % i for i in range(dim)], 
                    columns=['dummy_col_%d' % i for i in range(dim)]
            )
            df.to_csv(flpth, sep='\t', compression='gzip')

        ##
        def has_n_htmls(dir, n):
            return len([flnm for flnm in os.listdir(report_dir) if flnm.endswith('.html')]) == n
    
        ## make `run_dir`
        run_dir = os.path.join(self.shared_folder, 'test_report_' + str(uuid.uuid4()))
        os.mkdir(run_dir)



        ###
        ### heatmap enigma17770by511
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_enigma17770by511')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_enigma17770by511 = '/kb/module/test/data/by_dataset_input/enigma17770by511/return/PICRUSt2_output/pathways_out/path_abun_unstrat.tsv.gz'
            png_flpth = os.path.join(fig_dir, 'heatmap_enigma17770by511.png')
            html_flpth = os.path.join(report_dir, 'heatmap_enigma17770by511.html')

            do_heatmap(tsvgz_flpth_enigma17770by511, png_flpth, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### heatmap enigma50by30
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_enigma50by30')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_enigma50by30 = '/kb/module/test/data/by_dataset_input/enigma50by30/return/PICRUSt2_output/pathways_out/path_abun_unstrat.tsv.gz'
            png_flpth = os.path.join(fig_dir, 'heatmap_enigma50by30.png')
            html_flpth = os.path.join(report_dir, 'heatmap_enigma50by30.html')

            do_heatmap(tsvgz_flpth_enigma50by30, png_flpth, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### heatmap random
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_random_2k')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_random = os.path.join(report_dir, 'random.tsv.gz')
            png_flpth = os.path.join(fig_dir, 'heatmap_random.png')
            html_flpth = os.path.join(report_dir, 'heatmap_random.html')

            write_random_tsv(tsvgz_flpth_random, dim=2000, max=100)

            do_heatmap(tsvgz_flpth_random, png_flpth, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### HTMLReportWriter with 1 tsvgz
        with self.subTest():

            cmd_l = ['wingardium leviosa']

            tsvgz_flpth_l = [
                    tsvgz_flpth_enigma17770by511,

            ]

            report_dir = os.path.join(run_dir, 'report_HTMLReportWriter_1heatmap')

            HTMLReportWriter(cmd_l, tsvgz_flpth_l, report_dir).write()

            self.assertTrue(has_n_htmls(report_dir, 2))


        ###
        ### HTMLReportWriter with 3 tsvgz
        with self.subTest():

            cmd_l = ['expecto patronum', 'accio heatmap', 'luminos']

            report_dir = os.path.join(run_dir, 'report_HTMLReportWriter_3heatmaps')
            os.mkdir(report_dir)

            # copy tsvgzs into `report_dir` with unique names
            # since some have same filenames

            tsvgz_flpth_src_l = [
                    tsvgz_flpth_enigma17770by511,
                    tsvgz_flpth_enigma50by30,
                    tsvgz_flpth_random,
            ]

            tsvgz_flpth_dst_l = [os.path.join(report_dir, flnm) for flnm in [
                    'enigma17770by511.tsv.gz',
                    'enigma50by30.tsv.gz',
                    'random.tsv.gz',
            ]]

            for src, dst in zip(tsvgz_flpth_src_l, tsvgz_flpth_dst_l):
                shutil.copyfile(src, dst)

            HTMLReportWriter(cmd_l, tsvgz_flpth_dst_l, report_dir).write()

            self.assertTrue(has_n_htmls(report_dir, 4))


        ###
        ### HTMLReportWriter, heatmapping throws error
        with self.subTest():
            
            self.assertTrue(len(var.warnings) == 0)

            cmd_l = ['avada kedavra'] # not PG but signifies that this should fail

            tsvgz_flpth_l = [
                '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.fna', # not a tsvgz
            ]

            report_dir = os.path.join(run_dir, 'report_HTMLReportWriter_error')

            HTMLReportWriter(cmd_l, tsvgz_flpth_l, report_dir).write()

            self.assertTrue(len(var.warnings) == 1)
            self.assertTrue(has_n_htmls(report_dir, 1))






####################################################################################################
####################################################################################################
############################# INTEGRATION TESTS ####################################################
####################################################################################################
####################################################################################################

 
    ####################
    ####################
    #@patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30')) # ?
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma50by30'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma50by30'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_has_row_AttributeMapping(self):
        ret = self.serviceImpl.run_picrust2_pipeline(
            self.ctx, {
                **self.params_ws,
                'amplicon_matrix_upa': enigma50by30,
                'output_name': 'an_output_name',
            }
        )

        self.assertTrue(len(var.objects_created) == 8, var.objects_created) 


    ####################
    ####################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *args: get_mock_dfu('enigma50by30_noAttrMaps_noSampleSet'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma50by30_noAttrMaps_noSampleSet'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma50by30_noAttrMaps_noSampleSet'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_has_no_row_AttributeMapping(self):
        '''
        No taxonomy doesn't actually matter here
        The naming is just being consistent with kb_faprotax
        '''
        ret = self.serviceImpl.run_picrust2_pipeline(
            self.ctx, {
                **self.params_ws,
                'amplicon_matrix_upa': enigma50by30_noAttrMaps_noSampleSet,
                'output_name': 'an_output_name',
            }
        )

        self.assertTrue(len(var.objects_created) == 3, var.objects_created)

    ####################
    ####################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma50by30')) # ?
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma50by30'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma50by30'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_FP_options_wSampleSet(self):
        '''
        Creating amplicon FPs logic is independent of creating sample FPs logic
        '''

        with self.subTest():
            ret = self.serviceImpl.run_picrust2_pipeline(
                self.ctx, {
                    **self.params_ws,
                    'amplicon_matrix_upa': enigma50by30,
                    'fp_options': {
                        'create_amplicon_fps': False,
                        'create_sample_fps': False,
                    }
                }
            )

            self.assertTrue(len(var.objects_created) == 2, var.objects_created) 


        with self.subTest():
            ret = self.serviceImpl.run_picrust2_pipeline(
                self.ctx, {
                    **self.params_ws,
                    'amplicon_matrix_upa': enigma50by30,
                    'fp_options': {
                        'create_amplicon_fps': True,
                        'create_sample_fps': True,
                    }
                }
            )

            self.assertTrue(len(var.objects_created) == 8, var.objects_created) 


    ####################
    ####################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *args: get_mock_dfu('enigma50by30_noAttrMaps_noSampleSet'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma50by30_noAttrMaps_noSampleSet'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma50by30_noAttrMaps_noSampleSet'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_FP_options_noSampleSet(self):

        with self.subTest():
            ret = self.serviceImpl.run_picrust2_pipeline(
                self.ctx, {
                    **self.params_ws,
                    'amplicon_matrix_upa': enigma50by30_noAttrMaps_noSampleSet,
                    'fp_options': {
                        'create_amplicon_fps': False,
                        'create_sample_fps': False,
                    }
                }
            )

            self.assertTrue(len(var.objects_created) == 0, var.objects_created) 

        with self.subTest():
            ret = self.serviceImpl.run_picrust2_pipeline(
                self.ctx, {
                    **self.params_ws,
                    'amplicon_matrix_upa': enigma50by30_noAttrMaps_noSampleSet,
                    'fp_options': {
                        'create_amplicon_fps': False,
                        'create_sample_fps': True,
                    }
                }
            )

            self.assertTrue(len(var.objects_created) == 0, var.objects_created) 



    ####################
    ####################
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *a: get_mock_dfu('enigma17770by511'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.GenericsAPI', new=lambda *a, **k: get_mock_gapi('enigma17770by511'))
    @patch('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('enigma17770by511'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.FunctionalProfileUtil', new=lambda *a, **k: get_mock_fpu(''))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *a, **k: get_mock_kbr())
    def test_large_dataset(self):

        ret = self.serviceImpl.run_picrust2_pipeline(
            self.ctx, {
                **self.params_ws,
                'amplicon_matrix_upa': enigma17770by511,
                    'fp_options': {
                        'create_amplicon_fps': False,
                        'create_sample_fps': False,
                    },
            }
        )

        self.assertTrue(len(var.objects_created) == 2, var.objects_created) 
        

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    """
    #####
    #####
    def setUp(self):
        self.snapshot1 = tracemalloc.take_snapshot()
        top_stats = self.snapshot1.statistics('lineno')

        print("[ Top 10 ]")
        for stat in top_stats[:10]:
            print(stat)

    #####
    #####
    def tearDown(self):
        snapshot2 = tracemalloc.take_snapshot()

        ##
        top_stats = snapshot2.compare_to(self.snapshot1, 'lineno')

        print("[ Top 10 differences ]")
        for stat in top_stats[:10]:
            print(stat)

        ##
        top_stats = snapshot2.statistics('traceback')

        # pick the biggest memory block
        stat = top_stats[0]
        print("%s memory blocks: %.1f KiB" % (stat.count, stat.size / 1024))
        for line in stat.traceback.format():
            print(line)
    """

    #####
    #####
    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_PICRUSt2'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_PICRUSt2',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.wsName = 'kb_PICRUSt2_' + str(uuid.uuid4())                                                 
        cls.wsId = cls.wsClient.create_workspace({'workspace': cls.wsName})[0]                      
        cls.params_ws = {                                                                           
            'workspace_id': cls.wsId,                                                               
            'workspace_name': cls.wsName,                                                           
            }                                                                                       
        cls.serviceImpl = kb_PICRUSt2(cls.cfg)
        cls.shared_folder = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def list_tests(cls):
        return [key for key, value in cls.__dict__.items() if type(key) == str and key.startswith('test') and callable(value)]

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        print('Tests run:', cls.list_tests())
        dec = '!!!' * 300
        print(dec, "DON'T FORGET TO SEE HTML(S)", dec)

    def shortDescription(self):
        '''Override unittest using test*() docstrings in lieu of test*() method name in output summary'''
        return None

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!! select what to run !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
When you just want to run certain tests,
e.g., filter to tests in `run_tests`

Comment out parts like `delattr` to deactivate
'''
integration_tests = [
    'test_has_row_AttributeMapping', 'test_has_no_row_AttributeMapping',
    'test_FP_options_wSampleSet', 'test_FP_options_noSampleSet',
    'test_large_dataset',
]
unit_tests = [
    'test_run_check', 'test_OutfileWrangler', 
    'test_AmpliconSet_and_AmpliconMatrix', 'test_AttributeMapping',
    'test_report', 'test_large_heatmap', 'test_small_heatmap',
]
large_tests = [
    'test_large_dataset', 'test_large_heatmap',
]
run_tests = [
    'test_FP_options_noSampleSet',
]
norun_tests = [
    'test_large_dataset',
] 

for key, value in kb_PICRUSt2Test.__dict__.copy().items():
    if key.startswith('test') and callable(value):
        if key not in run_tests:
            delattr(kb_PICRUSt2Test, key)
            pass





