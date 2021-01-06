import os
import itertools
import pandas as pd
import numpy as np
import logging
import time

from .config import Var
from ..util.debug import dprint



####################################################################################################
####################################################################################################
def check_dropped_sample_ids(tsv_flpth, amp_mat):

    df_partial = pd.read_csv(tsv_flpth, sep='\t', index_col=0).T
    id_l_full = amp_mat.obj['data']['col_ids']

    assert sorted(df_partial.index) == sorted(id_l_full), '%d vs %d' % (len(df_partial.index), len(id_l_full))




####################################################################################################
####################################################################################################
def check_dropped_amplicon_ids(tsv_flpth, amp_mat):
    '''
    Check dropped amplicons (due to not aligning enough or too large NSTI)

    Parameters
    ----------
    tsv_flpth - the thing to check
    nsti_tsvgz_flpth - has info about dropped amplicons
    amp_mat - has info about all amplicons
    '''

    logging.info('Starting checking dropped amplicons for TSV %s' % tsv_flpth)

    #
    id_l_partial = pd.read_csv(tsv_flpth, sep='\t', index_col=0).index
    id_l_full = amp_mat.obj['data']['row_ids']

    # check align/nsti
    dropped_align, dropped_nsti = _get_dropped_ids(amp_mat)

    difference0 = sorted(set(id_l_partial) - set(id_l_full))
    difference1 = sorted(set(id_l_full) - set(id_l_partial))
    assert difference0 == []
    assert difference1 == sorted(dropped_align) or difference1 == sorted(dropped_align + dropped_nsti)



####################################################################################################
####################################################################################################
def _get_dropped_ids(amp_mat, nsti_max=2):
    '''
    Parse amplicon ids that were dropped due to not aligning to
    or being too distant from
    reference genomes
    '''
    # get this from globals because ...
    # passing files is messy and out_dir should be app global anyway
    nsti_flpth = os.path.join(Var.out_dir, 'marker_predicted_and_nsti.tsv.gz')

    df = pd.read_csv(nsti_flpth, sep='\t', index_col=0, header=0)

    ids_all = amp_mat.obj['data']['row_ids']
    dropped_align = list(set(ids_all) - set(df.index))
    dropped_nsti = list(df.index[df['metadata_NSTI'] > nsti_max])

    # should be disjoint
    assert len(set(dropped_align) & set(dropped_nsti)) == 0
    assert len(set(dropped_nsti) & set(dropped_align)) == 0, len(set(dropped_nsti) - set(dropped_align))
    
    return dropped_align, dropped_nsti





####################################################################################################
####################################################################################################
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


####################################################################################################
####################################################################################################
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

    id_x_desc_df = do_code2desc(id_x_code_df, Var.metacyc_pathway_code2desc_tsvgz)

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


