import shutil
import logging
import pandas as pd
import numpy as np
import os
import sys
import gzip
import json

from .dprint import dprint
from .config import var


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)




####################################################################################################
####################################################################################################
####################################################################################################

class AmpliconMatrix:

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()


    def _get_obj(self):
        logging.info('Loading AmpliconMatrix object')

        obj = var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.row_attrmap_upa = obj['data'][0]['data'].get('row_attributemapping_ref')
        self.obj = obj['data'][0]['data']
        
        if 'run_dir' in var: # comment directory with AmpMat name. optional since unit tests may not have run_dir
            dprint('touch %s' % os.path.join(var.run_dir, '#' + self.name), run='cli')


    def to_seq_abundance_table(self, flpth):
        logging.info(f"Writing sequence abundance table to %s" % flpth)

        data = np.array(self.obj['data']['values'], dtype=float)
        row_ids = self.obj['data']['row_ids']
        col_ids = self.obj['data']['col_ids']

        data = pd.DataFrame(
            data, 
            index=row_ids, # ASV Ids 
            columns=col_ids # sample names
            )
        data.index.name = "ASV_Id"
        data.to_csv(flpth, sep='\t')

    def to_fasta(self, flpth):
        fetched_flpth = var.gapi.fetch_sequence(self.upa)
        shutil.copyfile(fetched_flpth, flpth)

    def save(self, name=None):
        logging.info('Saving AmpliconMatrix')

        info = var.dfu.save_objects(
            {'id': var.params['workspace_id'],
             "objects": [{
                 "type": "KBaseMatrices.AmpliconMatrix", # TODO version
                 "data": self.obj,
                 "name": name if name is not None else self.name,
                 "extra_provenance_input_refs": [self.upa]
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new



####################################################################################################
####################################################################################################
####################################################################################################

class AttributeMapping:

    def __init__(self, upa, amp_mat):
        self.upa = upa
        self.amp_mat = amp_mat
        self._get_obj()

    def _get_obj(self):
        logging.info('Loading AttributeMapping object')

        obj = var.dfu.get_objects({
            'object_refs': ['%s;%s' %(self.amp_mat.upa, self.upa)]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']


    def update_attribute(self, ind: int, id2attr_d: dict):
        '''
        Update attribute at index `ind` using mapping `id2attr_d`
        '''
        for id, attr in id2attr_d.items():
            self.obj['instances'][id][ind] = attr


    def get_attribute_slot_warn(self, attribute, source) -> int:
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
                msg = (
                    'Overwriting attribute `%s` with source `%s` '
                    'in row AttributeMapping with name `%s`'
                    % (attribute, source, self.name)
                )
                logging.warning(msg)
                var.warnings.append(msg)
                return ind

        # append slot to `attributes`
        self.obj['attributes'].append({
            'attribute': attribute,
            'source': source,
            })

        # append slots to `instances` 
        for attr_l in self.obj['instances'].values():
            attr_l.append('')

        return len(attr_l) - 1


    def save(self):
        logging.info('Saving AttributeMapping')
        
        info = var.dfu.save_objects(
            {'id': var.params['workspace_id'],
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
class Report:
    '''For facilitating testing'''

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()

    def _get_obj(self):

        logging.info('Loading AttributeMapping object')

        obj = var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.obj = obj['data'][0]['data']




