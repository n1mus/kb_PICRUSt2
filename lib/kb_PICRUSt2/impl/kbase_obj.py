import shutil
import logging
import pandas as pd
import numpy as np
import os
import sys
import gzip
import json
import functools

from .config import Var
from .error import * # custom Exceptions
from ..util import validate as vd
from ..util.debug import dprint


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
class AmpliconMatrix:

####################################################################################################
####################################################################################################
    def __init__(self, upa):
        self.upa = upa
        self._get_obj()


####################################################################################################
####################################################################################################
    def _get_obj(self):
        logging.info('Loading AmpliconMatrix object')

        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.row_attrmap_upa = obj['data'][0]['data'].get('row_attributemapping_ref')
        self.obj = obj['data'][0]['data']
        
        if 'run_dir' in Var: # comment directory with AmpMat name. optional since unit tests may not have run_dir
            dprint('touch %s' % os.path.join(Var.run_dir, '#' + self.name), run='cli')


####################################################################################################
####################################################################################################
    def to_seq_abundance_table(self, flpth):
        '''
        Prerequisite: validate first with `validate_amplicon_abundance_data`
        '''
        logging.info(f"Writing sequence abundance table to %s" % flpth)

        def copy_replace(nested: list, rep=None):
            nested_ = []
            for l in nested:
                l_ = l.copy()
                nested_.append(l_)
                for i, e in enumerate(l_):
                    if e == '':
                        l_[i] = None
            return nested_

        data = np.array(copy_replace(self.obj['data']['values']), dtype=float) 
        row_ids = self.obj['data']['row_ids']
        col_ids = self.obj['data']['col_ids']

        data = pd.DataFrame(
            data, 
            index=row_ids, # Amplicon Ids 
            columns=col_ids # sample names
            )
        data.index.name = Var.amplicon_header_name
        data.to_csv(flpth, sep='\t', float_format='%g')

####################################################################################################
####################################################################################################
    def to_fasta(self, flpth):
        fetched_flpth = Var.gapi.fetch_sequence(self.upa)
        shutil.copyfile(fetched_flpth, flpth)

####################################################################################################
####################################################################################################
    def validate_amplicon_abundance_data(self):
        '''
        Can't be all missing
        Should be count data, missing allowed

        Data already restricted to int, float, None
        Because of KBase float types, which this is composed of,
        don't have to worry about complex, inf, etc.
        '''
        a = np.array(self.obj['data']['values'])

        # Can't be all missing
        if vd.get_num_missing(a) == a.size:
            raise vd.ValidationException(
                'Input AmpliconMatrix cannot have all missing matrix values'
            )

        a = vd.as_numeric(a, dtype=float) # casting as int will truncate floats into ints

        base_msg = (
            'Input AmpliconMatrix must have count data (missing values allowed) in matrix. '
        )

        # Integer
        if not vd.is_int_like(a):
            raise vd.ValidationException(
                base_msg + 'Non-integer detected'
            )

        # Gte 0
        if np.any(np.round(a) < 0): # allow for small negative deltas
            raise vd.ValidationException(
                base_msg + 'Negative value detected'
            )


####################################################################################################
####################################################################################################
    def _swap_ids(self, id2attr: dict, axis='row') -> dict:
        '''
        `id2attr` will be AmpliconMatrix ids to attribute
        Swap those ids out for the AttributeMapping ids
        '''

        if f'{axis}_mapping' not in self.obj:
            msg = (
                'Dude this object has a %s_attributemapping_ref '
                'and needs a %s_mapping. Letting it slide for now.'
                % (axis, axis)
            )
            logging.warning(msg)
            Var.warnings.append(msg)
            return id2attr

        id2attr = {
            self.obj[f'{axis}_mapping'][id]: attr
            for id, attr in id2attr.items()
        }

        return id2attr


####################################################################################################
####################################################################################################
    def save(self, name=None):
        logging.info('Saving AmpliconMatrix')

        upa_new = Var.gapi.save_object({
            'obj_type': 'KBaseMatrices.AmpliconMatrix', # TODO version
            'obj_name': name if name is not None else self.name,
            'data': self.obj,
            'workspace_id': Var.params['workspace_id'],
        })['obj_ref']

        return upa_new



####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
class AttributeMapping:

####################################################################################################
####################################################################################################
    def __init__(self, upa, amp_mat):
        '''
        Needs amp_mat for root ref and id mapping
        '''
        self.upa = upa
        self.amp_mat = amp_mat
        self._get_obj()

####################################################################################################
####################################################################################################
    def _get_obj(self):
        logging.info('Loading AttributeMapping object')

        obj = Var.dfu.get_objects({
            'object_refs': ['%s;%s' %(self.amp_mat.upa, self.upa)]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']


####################################################################################################
####################################################################################################
    def map_update_attribute(self, ind: int, id2attr: dict):
        '''
        Update attribute at index `ind` using mapping `id2attr`
        '''
        id2attr = self.amp_mat._swap_ids(id2attr)

        for id, attr in id2attr.items():
            self.obj['instances'][id][ind] = attr


####################################################################################################
####################################################################################################
    def get_add_attribute_slot(self, attribute, source) -> int:
        '''
        Get attribute slot matching both `attribute` and `source`
        
        Return the index of that slot, 
        which corresponds to both obj['attributes'] and each list in obj['instances']

        If slot matching both `attribute` and `source` does not exist
        * if `create=True` add slot for it
        * if `create=False` return -1
        '''
        
        # check if already exists
        for ind, attr_d in enumerate(self.obj['attributes']):
            if attr_d['attribute'] == attribute and attr_d['source'] == source:
                return ind, True

        # append slot to `attributes`
        self.obj['attributes'].append({
            'attribute': attribute,
            'source': source,
        })

        # append slots to `instances` 
        for attr_l in self.obj['instances'].values():
            attr_l.append(None)

        return len(attr_l) - 1, False


####################################################################################################
####################################################################################################
    def save(self):
        logging.info('Saving AttributeMapping')
        
        info = Var.dfu.save_objects(
            {'id': Var.params['workspace_id'],
             "objects": [{
                 "type": "KBaseExperiments.AttributeMapping", # TODO version
                 "data": self.obj,
                 "name": self.name,
                 "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new








####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
class Report:
    '''For facilitating testing'''

####################################################################################################
####################################################################################################
    def __init__(self, upa):
        self.upa = upa
        self._get_obj()

####################################################################################################
####################################################################################################
    def _get_obj(self):

        logging.info('Loading AttributeMapping object')

        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.obj = obj['data'][0]['data']


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
class KBaseObject:

####################################################################################################
####################################################################################################
    def __init__(self, upa):
        self.upa = upa
        self._get_obj()

####################################################################################################
####################################################################################################
    def _get_obj(self):
        logging.info('Loading KBase object data')
        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']


