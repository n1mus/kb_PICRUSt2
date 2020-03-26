# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser
import uuid

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2
from kb_PICRUSt2.kb_PICRUSt2Server import MethodContext
from kb_PICRUSt2.authclient import KBaseAuth as _KBaseAuth

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
                'amplicon_set_upa': enigmaFirst50_amp_set_upa,
                **self.params_ws,
                **params_debug,
                }
            )

        row_attrmap = AttributeMapping(Var.objects_created[0])
        instances = row_attrmap.obj['instances']
        attributes = row_attrmap.obj['attributes']





        



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
        cls.wsName = 'kb_faprotax_' + str(uuid.uuid4())                                                 
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
