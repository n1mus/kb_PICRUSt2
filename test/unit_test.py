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
from kb_PICRUSt2.impl.error import * # Exceptions
from kb_PICRUSt2.util.debug import dprint
from kb_PICRUSt2.util.validate import ValidationException
from mock import *


class kb_PICRUSt2Test(unittest.TestCase):

####################################################################################################
####################################################################################################
    @patch.dict('kb_PICRUSt2.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('dummy_10by8')})
    def test_appfile(self):
        '''
        '''
        logging.info('Testing with test_appfile')

        run_dir = os.path.join('/kb/module/work/tmp', 'test_OutfileWranger_' + str(uuid.uuid4()))
        os.mkdir(run_dir)

        ## Test `appfile.parse_picrust2_traits` ##
        with self.subTest('Test appfile.parse_picrust2_traits'):
        
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

            self.assertTrue(id2traits_d == ans)



####################################################################################################
####################################################################################################
    @patch.dict('kb_PICRUSt2.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('dummy_10by8')})
    def test_AmpliconMatrix(self):
        '''
        Combine AmpliconSet and AmpliconMatrix since they mostly write input files
        '''
        logging.info('Testing with test_AmpliconMatrix')

        # set up `run_dir`
        Var.run_dir = os.path.join(self.shared_folder, 'test_AmpliconMatrix_' + str(uuid.uuid4()))
        os.mkdir(Var.run_dir)

        amp_mat = AmpliconMatrix(dummy_10by8_AmpMat)

        # write

        seq_flpth = os.path.join(Var.run_dir, 'study_seqs.fna')
        seq_abundance_table_flpth = os.path.join(Var.run_dir, 'study_seqs.tsv')

        amp_mat.to_seq_abundance_table(seq_abundance_table_flpth)

        # compare
        
        seq_ref_flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.fna'
        seq_abundance_table_ref_flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.tsv'

        with open(seq_abundance_table_flpth) as f1:
            with open(seq_abundance_table_ref_flpth) as f2:
                self.assertTrue(f1.read() == f2.read())

    
####################################################################################################
####################################################################################################
    @patch.dict('kb_PICRUSt2.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('dummy_10by8')})
    def test_AmpliconMatrix_validation(self):
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
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data() 

        amp_mat.obj['data']['values'] = [-1]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data() 

        amp_mat.obj['data']['values'] = [None, None, None]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data() 

        amp_mat.obj['data']['values'] = [None]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, -1]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, None, 1.00001]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [-1.0, 0, 1319]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0, 1, 2, 3, 4.5]
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data()

        amp_mat.obj['data']['values'] = [None, 0.0, 1.0, 2.0, 3.0, 4.00001] # 4.00001 would pass with np.allclose default rtol
        with self.assertRaises(ValidationException): amp_mat.validate_amplicon_abundance_data()




####################################################################################################
####################################################################################################
    @patch.dict('kb_PICRUSt2.impl.kbase_obj.Var', values={'dfu': get_mock_dfu('dummy_10by8'), 'warnings': []})
    def test_AttributeMapping(self):
        '''
        Mostly writing attributes
        '''
        # TODO test with non-identical row_mapping

        logging.info('Testing with test_AttributeMapping')

        amp_mat = AmpliconMatrix(dummy_10by8_AmpMat)
        attr_map = AttributeMapping(dummy_10by8_AttrMap, amp_mat)

        ##
        ## write new attribute/source
        ind_0, overwrite = attr_map.get_add_attribute_slot('biome', 'testing')
        self.assertTrue(ind_0 == 2)
        self.assertTrue(overwrite is False)

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

        self.assertTrue(attr_map.obj['instances']['amplicon_id_4'][ind_0] == 'dummy0')

        ##
        ## overwrite attribute/source
        ind_1, overwrite = attr_map.get_add_attribute_slot('celestial body', 'upload')
        self.assertTrue(ind_1 == 0)
        self.assertTrue(overwrite is True)

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
            self.assertTrue(len(attr_l) == num_attr)

        ## 
        ## check did not add dummy attribute to wrong slot
        ind_lftvr = list(set(range(num_attr)) - {ind_0, ind_1})

        for attr_l in attr_map.obj['instances']:
            for ind in ind_lftvr:
                self.assertTrue('dummy' not in attr_l[ind])


####################################################################################################
####################################################################################################
    @patch.dict('kb_PICRUSt2.impl.report.Var', values={'warnings': []})
    def test_large_heatmap(self):
        '''
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
        html_flpth = os.path.join(report_dir, 'heatmap_random.html')

        write_random_tsv(tsvgz_flpth, dim=7000, max=1500)

        do_heatmap(tsvgz_flpth, html_flpth)



####################################################################################################
####################################################################################################
    @patch.dict('kb_PICRUSt2.impl.report.Var', values={'warnings': []})
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
            html_flpth = os.path.join(report_dir, 'heatmap_enigma17770by511.html')

            do_heatmap(tsvgz_flpth_enigma17770by511, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### heatmap enigma50by30
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_enigma50by30')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_enigma50by30 = '/kb/module/test/data/by_dataset_input/enigma50by30/return/PICRUSt2_output/pathways_out/path_abun_unstrat.tsv.gz'
            html_flpth = os.path.join(report_dir, 'heatmap_enigma50by30.html')

            do_heatmap(tsvgz_flpth_enigma50by30, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### heatmap random
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_random_2k')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_random = os.path.join(report_dir, 'random.tsv.gz')
            html_flpth = os.path.join(report_dir, 'heatmap_random.html')

            write_random_tsv(tsvgz_flpth_random, dim=2000, max=100)

            do_heatmap(tsvgz_flpth_random, html_flpth)

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

            tsvgz_flpth_l = [
                    tsvgz_flpth_enigma17770by511,
                    tsvgz_flpth_enigma50by30,
                    tsvgz_flpth_random,
            ]

            HTMLReportWriter(cmd_l, tsvgz_flpth_l, report_dir).write()

            self.assertTrue(has_n_htmls(report_dir, 4))

        '''
        ###
        ### HTMLReportWriter, heatmapping throws error
        with self.subTest():
            
            self.assertTrue(len(Var.warnings) == 0)

            cmd_l = ['avada kedavra'] # not PG but signifies that this should fail

            tsvgz_flpth_l = [
                '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.fna', # not a tsvgz
            ]

            report_dir = os.path.join(run_dir, 'report_HTMLReportWriter_error')

            HTMLReportWriter(cmd_l, tsvgz_flpth_l, report_dir).write()

            self.assertTrue(len(Var.warnings) == 1)
            self.assertTrue(has_n_htmls(report_dir, 1))
        '''



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    @classmethod
    def setUpClass(cls):
        cls.shared_folder = '/kb/module/work/tmp/'

    @classmethod
    def list_tests(cls):
        return [key for key, value in cls.__dict__.items() if type(key) == str and key.startswith('test') and callable(value)]

    @classmethod
    def tearDownClass(cls):
        dec = '!!!' * 300
        print(dec, "DON'T FORGET TO SEE HTML(S)", dec)
        skipped_tests = list(set(all_tests) - set(cls.list_tests()))
        print('* All tests (%d): %s' % (len(all_tests), all_tests))
        print('* Tests skipped (%d): %s' % (len(skipped_tests), skipped_tests))
        print('* Tests run (%d): %s' % (len(cls.list_tests()), cls.list_tests()))

    def shortDescription(self):
        '''Override unittest using test*() docstrings in lieu of test*() method name in output summary'''
        return None

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!! select what to run !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
all_tests = []
for key, value in kb_PICRUSt2Test.__dict__.copy().items():
    if key.startswith('test') and callable(value):
        all_tests.append(key)

unit_tests = [
    'test_AmpliconMatrix_validation',
    'test_run_check', 'test_appfile', 
    'test_AmpliconMatrix', 'test_AttributeMapping',
    'test_report', 'test_large_heatmap', 'test_small_heatmap',
]
large_tests = [
    'test_large_dataset', 'test_large_heatmap',
]
run_tests = [
    'test_AmpliconMatrix_validation',
]

for test in all_tests:
        if test not in run_tests:
            #delattr(kb_PICRUSt2Test, test)
            pass





