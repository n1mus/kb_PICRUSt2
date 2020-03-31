# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser
import uuid
import pandas as pd

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2
from kb_PICRUSt2.kb_PICRUSt2Server import MethodContext
from kb_PICRUSt2.authclient import KBaseAuth as _KBaseAuth

from kb_PICRUSt2.util.varstash import Var
from kb_PICRUSt2.util.dprint import dprint
from kb_PICRUSt2.util.kbase_obj import AttributeMapping

from installed_clients.WorkspaceClient import Workspace


params_debug = {
    'skip_run': True,
    'skip_retFiles': True,
    }

TRAIT = 'PiCrust2 Traits'


enigma_amp_set_upa = "48255/26/3"
enigmaFirst50_amp_set_upa = '48402/6/2'

class kb_PICRUSt2Test(unittest.TestCase):

    def test(self):
        ret = self.serviceImpl.run_kb_PICRUSt2(
            self.ctx, {
                'amplicon_set_upa': enigma_amp_set_upa,
                **self.params_ws,
                #**params_debug,
                }
            )
        return

        row_attrmap = AttributeMapping(Var.objects_created[0])
        instances_d = row_attrmap.obj['instances']
        attribute_d_l = row_attrmap.obj['attributes']

        # find index in attribute list
        for i, attribute_d in enumerate(attribute_d_l):
            if attribute_d['attribute'] == 'PiCrust2 Traits':
                ind = i

        # id to attribute
        results_d = {id: attr_l[ind] for id, attr_l in instances_d.items()}

        # id to traits
        answers_d = self.parse_answers_file()

        html_l = []

        for id in results_d:
            assert id in answers_d

            res = results_d[id]
            ans = answers_d[id]

            if res != ans:
                res_l = res.split(':')
                ans_l = ans.split(':')

                subset = set(ans_l).issubset(res_l)

                all_l = list(set(res_l + ans_l))

                html = []

                for func in all_l:
                    if func not in ans_l:
                        func = '<b>' + func + '</b>'
                    elif func not in res_l:
                        func = '<i>' + func + '</i>'

                html = '<p>' + ':'.join(html) + '</p>'
                if subset:
                    html = '__ans<=res__' + html
                html_l.append(html)

        len_original = len(html_l)
        html_l = list(set(html_l))
        len_dedup = len(html_l)

        html_l.append(f'original num mismatches: {len_original}, dedup num mismatches: {len_dedup}')

        with open(f'/kb/module/work/tmp/{uuid.uuid4()}.html', 'w') as fp:
            fp.write('\n'.join(html_l))


    @staticmethod
    def parse_answers_file():
        answers_flpth = '/kb/module/test/data/OTUMetaData_reduced.tsv'
        answers_df = pd.read_csv(
            answers_flpth, sep='\t', header=0, index_col='#OTU ID', usecols=['#OTU ID', 'PiCrust2 Traits']).fillna('')
        answers_d = answers_df.to_dict(orient='index')
        answers_d = {key: value['PiCrust2 Traits'] for key, value in answers_d.items()}
        return answers_d      



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
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

