import logging
import pandas as pd
import numpy as np
import os
import sys
import gzip

from .dprint import dprint
from .varstash import Var




pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)



class AttributeMapping:

    def __init__(self, upa):
        self.upa = upa
        self._get_obj()

    def _get_obj(self):
        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })

        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']


    def parse_picrust2_traits(self, id_x_func_tsv_gz_flpth, func_to_desc_tsv_filepath) -> dict:

        map_flpth = os.path.join(Var.picrust2_pckg_dir, 
            'default_files/description_mapfiles/metacyc_pathways_info.txt.gz')

        map_df = pd.read_csv(map_flpth, sep='\t', header=None, index_col=0, compression='gzip')
        map_d = map_df.to_dict(orient='index')
        map_d = {key: value[1] for key, value in map_d.items()}

        id_x_func_df = pd.read_csv(id_x_func_tsv_gz_flpth, sep='\t', header=0, index_col='sequence', compression='gzip')
        id_x_desc_df = id_x_func_df.rename(columns=map_d)

        dprint('map_d', 'id_x_desc_df', run=locals())

        desc_npArr = np.array(id_x_desc_df.columns.tolist())
        traits_l = []

        for sequence, row in id_x_desc_df.iterrows():
            abun_l = list(row)
            nonzero_ind_npArr = np.nonzero(abun_l)[0]
            traits_l.append(':'.join(list(desc_npArr[nonzero_ind_npArr])))

        id_x_desc_df['traits'] = traits_l

        id2traits_df = id_x_desc_df[['traits']] # TODO
        id2traits_d = id2traits_df.to_dict(orient='index')
        id2traits_d = {key: value['traits'] for key, value in id2traits_d.items()}

        return id2traits_d
        

        
    def add_attribute(self, id_to_attr_d, attribute='PiCrust2 Traits'):
        for ind, attr_d in enumerate(self.obj['attributes']):
            if attr_d['attribute'] == attribute:
                attr_ind = ind
                break

        for id, attr_l in self.obj['instances'].items():
            attr_l[attr_ind] = id_to_attr_d.get(id, '')

        dprint('self.obj["instances"]', run=locals())



    def save(self):
        
        info = Var.dfu.save_objects(
            {'id': Var.params['workspace_id'],
             "objects": [{
                 "type": "KBaseExperiments.AttributeMapping",
                 "data": self.obj,
                 "name": self.name + '.PICRUSt2',
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new




class AmpliconSet:

    def __init__(self, upa, test=False):
        self.upa = upa
        self.test = test

        self._get_obj()
        self._to_fasta()



    def _get_obj(self):
        obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })
        
        self.name = obj['data'][0]['info'][1]
        self.obj = obj['data'][0]['data']
        self.amp_mat_upa = self.obj['amplicon_matrix_ref']


    def _to_fasta(self):
        seq_flpth = os.path.join(Var.sub_dir, 'study_seqs.fna')
        
        logging.info(f'Writing fasta to {seq_flpth}')

        amplicon_d = self.obj['amplicons']
        
        with open(seq_flpth, 'w') as fp:
            for i, (ASV_id, d) in enumerate(amplicon_d.items()):
                fp.write('>' + ASV_id + '\n')
                fp.write(d['consensus_sequence'] + '\n')

                if self.test and i > 20:
                    break
              
        self.seq_flpth = seq_flpth



    def update_amplicon_matrix_ref(self, amp_mat_upa_new):
        self.obj['amplicon_matrix_ref'] = amp_mat_upa_new


    def save(self):
        dprint('self.obj', run=locals())

        info = Var.dfu.save_objects(
            {'id': Var.params['workspace_id'],
             "objects": [{
                 "type": "KBaseExperiments.AmpliconSet",
                 "data": self.obj,
                 "name": self.name + '.FAPROTAX',
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new





class AmpliconMatrix:

    def __init__(self, upa, amp_set: AmpliconSet):
        self.upa = upa
        self.amp_set = amp_set

        self._get_obj()
        self._to_seq_abundance_table()


    def _get_obj(self):
        obj = Var.dfu.get_objects({
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

        self.seq_abundance_table_flpth = os.path.join(Var.sub_dir, 'study_seqs.tsv')

        data.to_csv(self.seq_abundance_table_flpth, sep='\t')


    def update_row_attributemapping_ref(self, row_attrmap_upa_new):
        self.obj['row_attributemapping_ref'] = row_attrmap_upa_new


    def save(self):
        dprint('self.obj', run=locals())

        info = Var.dfu.save_objects(
            {'id': Var.params['workspace_id'],
             "objects": [{
                 "type": "KBaseMatrices.AmpliconMatrix",
                 "data": self.obj,
                 "name": self.name + '.FAPROTAX',
             }]})[0]

        upa_new = "%s/%s/%s" % (info[6], info[0], info[4])

        return upa_new



       


