import logging
import pandas as pd
import numpy as np

from .dprint import *
from .varstash import Var




pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)



class AmpliconSet:

    def __init__(self, upa, test=True):
        self.upa = upa
        self.test = test

        self._get_obj()
        self._to_fasta()



    def _get_obj(self):
        self.obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })['data'][0]['data']

        self.amp_mat_upa = self.obj['amplicon_matrix_ref']


    def _to_fasta(self):
        seq_flpth = os.path.join(Var.sub_dir, 'study_seqs.fna')
        
        logging.info(f'Writing fasta to {seq_flpth}')

        amplicon_d = self.obj['amplicons']
        
        with open(seq_flpth, 'w') as fp:
            for i, (ASV_id, d) in enumerate(amplicon_d.items()):
                fp.write('>' + ASV_id + '\n')
                fp.write(d['consensus_sequence'] + '\n')

                if self.test and i > 5:
                    break
              
        self.seq_flpth = seq_flpth


    def get_amplicon_matrix_upa(self):
        return self.amp_mat_upa




class AmpliconMatrix:

    def __init__(self, upa, amp_set: AmpliconSet, test=True):
        self.upa = upa
        self.amp_set = amp_set

        self._get_obj()
        self._to_seq_abundance_table()


    def _get_obj(self):
        self.obj = Var.dfu.get_objects({
            'object_refs': [self.upa]
            })['data'][0]['data']


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


       




class FakeData:

    tax_path_l = [
"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Thermomonas;s__fusca",
"D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli;D_3__Bacillales;D_4__Bacillaceae;D_5__Bacillus;D_6__Bacillus sp. YZ5",
"Bacteria; Chlorobi; Chlorobia; Chlorobiales; Chlorobiaceae; Chlorobaculum; Chlorobaculum thiosulfatiphilum DSM249T",   
"Bacteria; Chlorobi; Chlorobia; Chlorobiales; Chlorobiaceae; Chlorobaculum",
"Bacteria; Chlorobi; Chlorobia; Chlorobiales; Chlorobiaceae; Chlorobaculum; uncultured bacterium",
"Bacteria;Proteobacteria;Alphaproteobacteria;Rhodobacterales;Rhodobacteraceae;Paracoccus;Fervidicola"            
        ]



