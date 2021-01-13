import json

from .config import Var
from ..util.debug import dprint



def flatten(d):
    '''
    At most 1 level nesting
    Removes any nested dict's key
    '''
    d1 = d.copy()
    for k, v in d.items():
        if isinstance(v, dict):
            for k1, v1 in d1.pop(k).items():
                d1[k1] = v1
    return d1


class Params:


    DEFAULTS = {
        'cog': 0,
        'ec': 1,
        'ko': 1,
        'pfam': 0,
        'tigrfam': 0,
        'pheno': 0,
        'metacyc': 1,
        'create_amplicon_fps': True,
        'create_sample_fps': True,
    }


    ALL = [
        'amplicon_matrix_upa',
        'output_name',
        #---
        'functions',
        'fp_options',
        #---
        'cog',
        'ec',
        'ko',
        'pfam',
        'tigrfam',
        'pheno',
        'metacyc',
        'create_amplicon_fps', 
        'create_sample_fps',
        #---
        'workspace_id',
        'workspace_name',
    ]

    def __init__(self, params):

        self._validate(params)
        params = flatten(params)
        
        ## Custom transformations to internal state ##
        for param in [
            'create_amplicon_fps',
            'create_sample_fps',
            'cog',
            'pfam',
            'tigrfam',
            'pheno',
        ]:
            self._rep_as_bool(params, param)
       
        ##
        self.params = params


    def _validate(self, params):
        '''
        None of params go directly to shell
        '''
        # TODO

        for p in params:
            if p not in self.ALL:
                raise Exception(p)

        for p in flatten(params):
            if p not in self.ALL:
                raise Exception(p)




    def __getitem__(self, key):
        '''
        For required params (e.g., input UPAs, workspace stuff)
        Should not use this for default-backed params
        as those can be left off params
        so use `getd` for those
        '''
        return self.params[key]

    def getd(self, key):
        '''
        For default-backed params (e.g., tunable numbers)
        Like `get`
        Return the user-supplied value, or the default value if none was supplied
        '''
        if key not in self.DEFAULTS:
            raise Exception('`params.getd(x)` only applicable to params with defaults')
        if key not in self.ALL:
            raise Exception('Not a valid param')

        return self.params.get(key, self.DEFAULTS[key])


    def __repr__(self) -> str:
        return 'Wrapper for params:\n%s' % (json.dumps(self.params, indent=4))


    @staticmethod
    def _rep_as_bool(params, param) -> None:
        if type(params.get(param)) is int:
            params[param] = True if params[param] == 1 else False
