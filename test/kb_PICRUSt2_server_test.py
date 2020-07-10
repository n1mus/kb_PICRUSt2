# -*- coding: utf-8 -*-
import os
import time
import unittest
from unittest.mock import patch
from configparser import ConfigParser
import uuid
import pandas as pd

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2, run_check, OutfileWrangler
from kb_PICRUSt2.kb_PICRUSt2Server import MethodContext
from kb_PICRUSt2.authclient import KBaseAuth as _KBaseAuth

from kb_PICRUSt2.util.config import var
from kb_PICRUSt2.util.dprint import dprint
from kb_PICRUSt2.util.kbase_obj import AmpliconSet, AmpliconMatrix, AttributeMapping, Report
from kb_PICRUSt2.util.error import *
from util.mock import *
from util.upa import *

from installed_clients.WorkspaceClient import Workspace


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



class kb_PICRUSt2Test(unittest.TestCase):

####################################################################################################
####################################################################################################
############################# UNIT TESTS ###########################################################
####################################################################################################
####################################################################################################

    def test_run_check(self):
        '''
        Test `run_check` which runs PICRUSt2 executable
        '''

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


    def test_OutfileWrangler(self):
        '''
        Test `OutfileWrangler.parse_picrust2_traits` mostly
        '''

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




    @patch.dict('kb_PICRUSt2.util.kbase_obj.var', values={'dfu': get_mock_dfu('dummy_10by8')})
    def test_AmpliconSet_and_AmpliconMatrix(self):
        '''
        Combine AmpliconSet and AmpliconMatrix since they mostly write input files
        '''
        # set up `run_dir`
        run_dir = os.path.join(self.shared_folder, 'test_AmpliconSet_and_AmpliconMatrix_' + str(uuid.uuid4()))
        os.mkdir(run_dir)

        amp_set = AmpliconSet(dummy_10by8)
        amp_mat = AmpliconMatrix(dummy_10by8_AmpMat)

        # write

        seq_flpth = os.path.join(run_dir, 'study_seqs.fna')
        seq_abundance_table_flpth = os.path.join(run_dir, 'study_seqs.tsv')

        amp_set.to_fasta(seq_flpth)
        amp_mat.to_seq_abundance_table(seq_abundance_table_flpth)

        # compare
        
        seq_ref_flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.fna'
        seq_abundance_table_ref_flpth = '/kb/module/test/data/by_dataset_input/dummy_10by8/return/study_seqs.tsv'
        
        with open(seq_flpth) as f1:
            with open(seq_ref_flpth) as f2:
                self.assertTrue(f1.read() == f2.read())

        with open(seq_abundance_table_flpth) as f1:
            with open(seq_abundance_table_ref_flpth) as f2:
                self.assertTrue(f1.read() == f2.read())



    @patch.dict('kb_PICRUSt2.util.kbase_obj.var', values={'dfu': get_mock_dfu('dummy_10by8'), 'warnings': []})
    def test_AttributeMapping(self):
        '''
        Mostly writing attributes
        '''
        attr_map = AttributeMapping(dummy_10by8_AttrMap)

        ##
        ## write new attribute/source
        ind_0 = attr_map.get_attribute_slot('biome', 'testing')
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
        ind_1 = attr_map.get_attribute_slot('celestial body', 'upload')
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

    def test_report(self):
        '''
        Stress test
        '''
        # make `run_dir`

        var.run_dir = os.path.join(self.shared_folder, 'test_report_' + str(uuid.uuid4()))
        os.mkdir(var.run_dir)



####################################################################################################
####################################################################################################
############################# INTEGRATION TESTS ####################################################
####################################################################################################
####################################################################################################

 
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *args: get_mock_dfu('17770'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('17770'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_has_row_AttributeMapping(self):
        ret = self.serviceImpl.run_picrust2_pipeline(
            self.ctx, {
                **self.params_ws,
                'amplicon_set_upa': _17770,
                'output_name': 'an_output_name',
            }
        )

        if do_patch:
            return

        rprt = Report(ret[0]['report_ref'])

        self.assertTrue(len(rprt.obj.objects_created) == 3) # AttributeMapping, AmpliconMatrix, and AmpliconSet
        self.assertTrue(len(rprt.obj.file_links) == 1)


    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.DataFileUtil', new=lambda *args: get_mock_dfu('secret'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.run_check', new=get_mock_run_check('secret'))
    @patch_('kb_PICRUSt2.kb_PICRUSt2Impl.KBaseReport', new=lambda *args, **kwargs: get_mock_kbr())
    def test_has_no_row_AttributeMapping(self):
        '''
        No taxonomy doesn't actually matter here
        The naming is just being consistent with kb_faprotax
        '''
        ret = self.serviceImpl.run_picrust2_pipeline(
            self.ctx, {
                **self.params_ws,
                'amplicon_set_upa': secret,
                'output_name': 'an_output_name',
            }
        )

        if do_patch:
            return

        rprt = Report(ret[0]['report_ref'])

        self.assertTrue(len(rprt.obj.objects_created) == 2) # AttributeMapping and AmpliconSet
        self.assertTrue(len(rprt.obj.file_links) == 1)



        

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
integration_tests = ['test_has_row_AttributeMapping', 'test_has_no_row_AttributeMapping']
unit_tests = ['test_run_check', 'test_OutfileWrangler', 'test_AmpliconSet_and_AmpliconMatrix', 'test_AttributeMapping']
run_tests = ['test_has_no_row_AttributeMapping']

for key, value in kb_PICRUSt2Test.__dict__.copy().items():
    if key.startswith('test') and callable(value):
        if key not in integration_tests:
            delattr(kb_PICRUSt2Test, key)
            pass





