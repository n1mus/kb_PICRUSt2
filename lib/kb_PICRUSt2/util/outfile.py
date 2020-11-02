import os
import itertools
import pandas as pd
import numpy as np
import logging
import time

from .dprint import dprint
from .config import var



class OutfileWrangler:
    '''
    Put these methods in a class since some call each other
    '''

    #####
    #####
    @staticmethod
    def check_dropped_sample_ids(tsv_flpth, amp_mat):

        index = var.tsvFlnm2Index[os.path.basename(tsv_flpth)]
        assert index == 'sample' 

        df_partial = pd.read_csv(tsv_flpth, sep='\t', index_col=0).T
        id_l_full = amp_mat.obj['data']['col_ids']

        assert sorted(df_partial.index) == sorted(id_l_full), '%d vs %d' % (len(df_partial.index), len(id_l_full))



    #####
    ##### TODO read/write GZ since IO takes like 2 min
    @staticmethod
    def check_pad_dropped_amplicon_ids(tsv_flpth, amp_mat):
        '''
        Pad dropped amplicons (due to not aligning enough or too large NSTI)
        with NaN


        Parameters
        ----------
        tsv_flpth - the thing to pad
        nsti_tsvgz_flpth - has info about dropped amplicons
        amp_mat - has info about all amplicons
        '''
 
        # not applicable to sample TSV
        index = var.tsvFlnm2Index[os.path.basename(tsv_flpth)]
        assert index == 'amplicon' 

        logging.info('Starting padding dropped amplicons for TSV %s' % tsv_flpth)
        t0 = time.time()

        #
        df_partial = pd.read_csv(tsv_flpth, sep='\t', index_col=0)
        id_l_full = amp_mat.obj['data']['row_ids']
 
        # no need to pad
        if df_partial.shape[0] == len(id_l_full):
            logging.info('No need for padding TSV %s' % tsv_flpth)
            return

        # check align/nsti padded
        dropped_align, dropped_nsti = OutfileWrangler._get_dropped_ids(amp_mat)
        dprint(
            'len(dropped_align)',
            'len(dropped_nsti)',
            'len(df_partial.index)',
            'len(id_l_full)',
            run=locals()
        )
        difference = sorted(set(id_l_full) - set(df_partial.index))
        assert difference == sorted(dropped_align) or difference == sorted(dropped_align + dropped_nsti)

        # full-sized DF of NaNs
        # in order of AmpliconMatrix
        df_full = pd.DataFrame(
            np.zeros((len(id_l_full), df_partial.shape[1])) * np.nan, 
            index=id_l_full,
            columns=df_partial.columns
        )

        # assign non-dropped values
        df_full.loc[df_partial.index, df_partial.columns] = df_partial.values

        #
        dprint(
            "tsv_flpth",
            "df_partial.shape",
            "len(set(id_l_full) - set(df_partial.index)) # num amp ids padded",
            "df_full.shape # padded TSV",
            "df_full.values.size # num el",
            "(df_full != 0).sum().sum() # num non-zero el",
            "(df_full != 0).sum().sum() / df_full.values.size # sparsity",
            run=locals()
        ) 

        # use %g to try to keep TSV/JSON small for integers
        df_full.to_csv(tsv_flpth, sep='\t', na_rep='', float_format='%g')
   
        logging.info('Finished padding TSV %s, took %.2fs' %(tsv_flpth, (time.time()-t0)))


    #####
    #####
    @staticmethod
    def _get_dropped_ids(amp_mat, nsti_max=2):
        '''
        Parse amplicon ids that were dropped due to not aligning to
        or being too distant from
        reference genomes
        '''
        # get this from globals because ...
        # passing files is messy and out_dir should be app global anyway
        nsti_flpth = os.path.join(var.out_dir, 'marker_predicted_and_nsti.tsv.gz')

        ids_all = amp_mat.obj['data']['row_ids']

        df = pd.read_csv(nsti_flpth, sep='\t', index_col=0, header=0)

        dropped_align = list(set(ids_all) - set(df.index))
        dropped_nsti = list(df.index[df['metadata_NSTI'] > nsti_max])

        # should be disjoint
        assert len(set(dropped_align) & set(dropped_nsti)) == 0
        assert len(set(dropped_nsti) & set(dropped_align)) == 0, len(set(dropped_nsti) - set(dropped_align))
        
        return dropped_align, dropped_nsti



    """
    #####
    #####
    @staticmethod
    def pad_0_vecs(tsv_flpth, amp_mat):
        '''
        PICRUSt2 drops ids/samples that are all 0s
        Restore them in original order here
        '''
        logging.info('Starting 0-padding TSV %s' % tsv_flpth)
        t0 = time.time()

        index = var.tsvFlnm2Index[os.path.basename(tsv_flpth)]
        df_partial = pd.read_csv(tsv_flpth, sep='\t', index_col=0)

        # orient PICRUSt2 output matrix as something vs. func
        if index == 'sample':
            df_partial = df_partial.T
            id_l_full = amp_mat.obj['data']['col_ids']
        elif index == 'amplicon':
            id_l_full = amp_mat.obj['data']['row_ids']
        else:
            raise Exception(index)

        # no need to pad
        if df_partial.shape[0] == len(id_l_full):
            logging.info('No need for 0-padding TSV %s' % tsv_flpth)
            return

        df_full = pd.DataFrame(
            np.zeros((len(id_l_full), df_partial.shape[1])), 
            index=id_l_full,
            columns=df_partial.columns
        )

        df_full.loc[df_partial.index, df_partial.columns] = df_partial.values

        dprint(
            "tsv_flpth",
            "df_partial.shape",
            "len(set(id_l_full) - set(df_partial.index)) # num %ss padded" % index,
            "df_full.shape # padded TSV",
            "df_full.values.size # num el",
            "(df_full != 0).sum().sum() # num non-zero el",
            run=locals()
        ) 

        # undo orient
        if index == 'sample':
            df_full = df_full.T

        #df_full.to_csv(tsv_flpth, sep='\t')
   
        logging.info('Finished 0-padding TSV %s, took %.2fs' %(tsv_flpth, (time.time()-t0)))

    """


    #####
    #####
    @staticmethod
    def do_code2desc(df: pd.DataFrame, code2desc_tsv_gz_flpth: str, code_in='col') -> pd.DataFrame:
        '''
        df - code is in index or column names
        code2desc_tsv_gz_flpth - two columns, code and description
        '''

        if code_in not in ['row', 'col']:
            raise Exception()

        if code_in == 'row':
            df = df.T

        # parse code2desc tsv
        code2desc_df = pd.read_csv(
                code2desc_tsv_gz_flpth,
                sep='\t', 
                names=['code', 'description'], 
                index_col='code', 
                compression='gzip'
        )

        # convert to dict
        code2desc_d = code2desc_df['description'].to_dict() # cast to Series first to avoid annoying nested dict

        df = df.rename(columns=code2desc_d)

        if code_in == 'row':
            df = df.T

        return df

    #####
    #####
    @staticmethod
    def parse_picrust2_traits(id_x_code_tsv_gz_flpth, dlm=',') -> dict:
        '''
        id_x_code_tsv_flpth - created by picrust2_pipeline.py
                              rows are amplicon ids, columns are MetaCyc codes
                              (leaves out any ids or codes with no hits)

        Also translated MetaCyc codes into descriptoins
        '''
                
        ##
        ## translate MetaCyc pathway codes to descriptions

        id_x_code_df = pd.read_csv(
                id_x_code_tsv_gz_flpth, 
                sep='\t', 
                header=0, 
                index_col='sequence', 
                compression='gzip'
        )

        id_x_desc_df = OutfileWrangler.do_code2desc(id_x_code_df, var.metacyc_pathway_code2desc_tsvgz)

        ##
        ## aggregate MetaCyc descriptions for each id

        desc_npArr = np.array(id_x_desc_df.columns.tolist())
        traitsStr_l = []

        for sequence, row in id_x_desc_df.iterrows():
            abun_l = list(row) 
            ind_nonzero_npArr = np.nonzero(abun_l)[0]
            traitsStr_l.append(dlm.join(list(desc_npArr[ind_nonzero_npArr])))

        id_x_desc_df['traits'] = traitsStr_l

        id2traits_d = id_x_desc_df['traits'].to_dict() # cast to Series first to avoid annoying nested dict

        return id2traits_d


