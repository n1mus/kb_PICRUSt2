import logging
import pandas as pd
import numpy as np
import os
import sys
import gzip

from .dprint import dprint
from .config import _globals




pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)


####################################################################################################
####################################################################################################
####################################################################################################

class AttributeMapping:

    def __init__(self, upa, mini_test=False):
        self.upa = upa
        self.mini_test = mini_test

        self._get_obj()

    def _get_obj(self):
        obj = _globals.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']


    def parse_picrust2_traits(self, id_x_code_tsv_gz_flpth) -> dict:
        '''
        id_x_code_tsv_flpth - created by picrust2_pipeline.py
                              rows are amplicon ids, columns are metacyc codes
                              (leaves out any ids or codes with no hits)
        '''

        ## parse and prep ds

        MAP_FLPTH = os.path.join(_globals.picrust2_pckg_dir, # metacyc code to description
            'default_files/description_mapfiles/metacyc_pathways_info.txt.gz') 

        map_df = pd.read_csv(MAP_FLPTH, sep='\t', header=None, index_col=0, compression='gzip')
        map_d = map_df[1].to_dict()

        id_x_code_df = pd.read_csv(id_x_code_tsv_gz_flpth, sep='\t', header=0, index_col='sequence', compression='gzip')
        id_x_desc_df = id_x_code_df.rename(columns=map_d)

        dprint('map_d', 'id_x_desc_df', max_lines=5, run=locals())

        # aggregate metacyc descriptions for each id

        desc_npArr = np.array(id_x_desc_df.columns.tolist())
        traits_l = []

        for sequence, row in id_x_desc_df.iterrows():
            abun_l = list(row) 
            nonzero_ind_npArr = np.nonzero(abun_l)[0]
            traits_l.append(':'.join(list(desc_npArr[nonzero_ind_npArr])))

        id_x_desc_df['traits'] = traits_l

        id2traits_d = id_x_desc_df['traits'].to_dict()

        return id2traits_d
        

        
    def update_attribute(self, id2attr_d, attribute, source):
        # find index of attribute
        for ind, attr_d in enumerate(self.obj['attributes']):
            if attr_d['attribute'] == attribute:
                attr_ind = ind
                break

        for id, attr_l in self.obj['instances'].items():
            attr_l[attr_ind] = id2attr_d.get(id, '')

        self.obj['attributes'][attr_ind]['source'] = source


    def add_attribute_slot(self, attribute):
        
        # check if already exists
        for attr_d in self.obj['attributes']:
            if attr_d['attribute'] == attribute:
                msg = 'Adding attribute slot %s to AttributeMapping with name %s, ' % (attribute, self.name) + \
                      'but that attribute already exists in object'
                logging.warning(msg)
                _globals.warnings.append(msg)
                return

        # append slot to `attributes`
        self.obj['attributes'].append({
            'attribute': attribute,
            })

        # append slots to `instances` 
        for _, attr_l in self.obj['instances'].items():
            attr_l.append('')



    def save(self):
        
        info = _globals.dfu.save_objects(
            {'id': _globals.params['workspace_id'],
             "objects": [{
                 "type": "KBaseExperiments.AttributeMapping",
                 "data": self.obj,
                 "name": self.name,
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new



####################################################################################################
####################################################################################################
####################################################################################################

class AmpliconSet:

    def __init__(self, upa, mini_test=False):
        self.upa = upa
        self.mini_test = mini_test

        self._get_obj()
        self._to_fasta()



    def _get_obj(self):
        obj = _globals.dfu.get_objects({
            'object_refs': [self.upa]
            })
        
        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']
        self.amp_mat_upa = self.obj['amplicon_matrix_ref']


    # TODO move to AttributeMapping
    def _to_fasta(self):
        seq_flpth = os.path.join(_globals.run_dir, 'study_seqs.fna')
        
        logging.info(f'Writing fasta to {seq_flpth}')

        amplicon_d = self.obj['amplicons']
        
        with open(seq_flpth, 'w') as fp:
            for i, (ASV_id, d) in enumerate(amplicon_d.items()):
                fp.write('>' + ASV_id + '\n')
                fp.write(d['consensus_sequence'] + '\n')

                if _globals.debug and self.mini_test and i > 20:
                    break
              
        self.seq_flpth = seq_flpth



    def update_amplicon_matrix_ref(self, amp_mat_upa_new):
        self.obj['amplicon_matrix_ref'] = amp_mat_upa_new


    def save(self, name=None):
        dprint('self.obj', run=locals())

        info = _globals.dfu.save_objects(
            {'id': _globals.params['workspace_id'],
             "objects": [{
                 "type": "KBaseExperiments.AmpliconSet",
                 "data": self.obj,
                 "name": name if name else self.name,
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new




####################################################################################################
####################################################################################################
####################################################################################################

class AmpliconMatrix:

    def __init__(self, upa):
        self.upa = upa

        self._get_obj()
        self._to_seq_abundance_table()


    def _get_obj(self):
        obj = _globals.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.row_attrmap_upa = obj['data'][0]['data']['row_attributemapping_ref']
        self.obj = obj['data'][0]['data']


    def _to_seq_abundance_table(self):

        logging.info(f"Parsing AmpliconMatrix and AmpliconSet data from object")

        data = np.array(self.obj['data']['values'], dtype=float)
        row_ids = self.obj['data']['row_ids']
        col_ids = self.obj['data']['col_ids']

        data = pd.DataFrame(
            data, 
            index=row_ids, # ASV Ids 
            columns=col_ids # sample names
            )
        data.index.name = "ASV_Id"

        self.seq_abundance_table_flpth = os.path.join(_globals.run_dir, 'study_seqs.tsv')

        data.to_csv(self.seq_abundance_table_flpth, sep='\t')


    def update_row_attributemapping_ref(self, row_attrmap_upa_new):
        self.obj['row_attributemapping_ref'] = row_attrmap_upa_new


    def save(self, name=None):
        dprint('self.obj', run=locals())

        info = _globals.dfu.save_objects(
            {'id': _globals.params['workspace_id'],
             "objects": [{
                 "type": "KBaseMatrices.AmpliconMatrix",
                 "data": self.obj,
                 "name": name if name else self.name,
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new



       


