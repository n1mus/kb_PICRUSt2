#!/home/sumin/anaconda3/bin/python

import pandas as pd
import numpy as np
import os
import sys
import logging

RESULTS_FLPTH = '/home/sumin/kbsdk-workspace/kb_PICRUSt2/test/data/PICRUSt2_output2/pathways_out/path_abun_predictions.tsv.gz'
ANSWERS_FLPTH = '/home/sumin/kbsdk-workspace/kb_PICRUSt2/test/data/OTUMetaData_reduced.tsv'
MAP_FLPTH = '/home/sumin/anaconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz'

pwd = os.path.dirname(os.path.realpath(__file__)) 
sys.path.append(pwd)
from dprint import dprint


def parse_answers():
    answers_df = pd.read_csv(
        ANSWERS_FLPTH, sep='\t', header=0, index_col='#OTU ID', usecols=['#OTU ID', 'PiCrust2 Traits']).fillna('')
    answers_d = answers_df['PiCrust2 Traits'].to_dict()
    
    return answers_d


def parse_map():
    '''
    Metacyc codes to Metacyc descriptions
    '''
    map_df = pd.read_csv(MAP_FLPTH, sep='\t', header=None, index_col=0, compression='gzip')
    map_d = map_df.to_dict(orient='index')
    map_d = {key: value[1] for key, value in map_d.items()}
    return map_d


def parse_results(all_keys: list):
    '''
    Result id to traits
    '''
    map_d = parse_map()

    id_x_func_df = pd.read_csv(RESULTS_FLPTH, sep='\t', header=0, index_col='sequence', compression='gzip')
    id_x_desc_df = id_x_func_df.rename(columns=map_d)

    desc_npArr = np.array(id_x_desc_df.columns.tolist())
    traits_l = []

    for sequence, row in id_x_desc_df.iterrows():
        abun_l = list(row)
        nonzero_ind_npArr = np.nonzero(abun_l)[0]
        traits_l.append(':'.join(list(desc_npArr[nonzero_ind_npArr])))

    id_x_desc_df['traits'] = traits_l
    id2traits_d = id_x_desc_df['traits'].to_dict()
    
    for key in all_keys:
        if key not in id2traits_d:
            id2traits_d[key] = ''

    return id2traits_d






if __name__ == '__main__':

    out_flpth = sys.argv[1] # diff.html

    answers_d = parse_answers()
    results_d = parse_results(answers_d.keys())

    keys_ans = sorted(list(answers_d.keys()))
    keys_res = sorted(list(results_d.keys()))

    keys_all = sorted(list(set(keys_ans + keys_res)))

    assert set(keys_res).issubset(keys_all)
    assert keys_ans == keys_all


    html_l = []

    exact_match_ct = 0
    outOfOrder_ct = 0
    subset_res_of_ans_ct = 0
    subset_ans_of_res_ct = 0
    venn_diagram_ct = 0
    no_overlap_ct = 0

    traits_ans = [answers_d[key] for key in keys_ans]



    for id in keys_all:

        res = results_d.get(id, '')
        ans = answers_d[id]

        # Case: Exact Match
        if res == ans:
            exact_match_ct += 1
        
        else:
            res_l = res.split(':')
            ans_l = ans.split(':')
            
            # Case: Out of Order (match)
            if sorted(res_l) == sorted(ans_l):
                outOfOrder_ct += 1
                continue
            
            # Case: Res <= Ans
            if set(res_l).issubset(ans_l):
                subset_res_of_ans_ct += 1
            # Case: Ans <= Res
            elif set(ans_l).issubset(res_l):
                subset_ans_of_res_ct += 1
            elif set(ans_l).intersection(set(res_l)):
                venn_diagram_ct += 1
            # Case: No overlap
            else:
                no_overlap_ct += 1

            all_l = sorted(list(set(res_l + ans_l)))

            line = []

            for func in all_l:
                if func in res_l and func not in ans_l:
                    func = '<font color="blue"><b>' + func + '</b></font>'
                elif func in ans_l and func not in res_l:
                    func = '<font color="red"><b>' + func + '</b></font>'
                line.append(func)

            line = '<p>' + ':'.join(line) + '</p>'
            html_l.append(line)

    mismatch_len_original = len(html_l)
    html_l = list(set(html_l))
    mismatch_len_dedup = len(html_l)


    dprint('exact_match_ct', 'outOfOrder_ct', 'subset_res_of_ans_ct', 'subset_ans_of_res_ct', 
         'venn_diagram_ct', 'no_overlap_ct', run=locals())

    dprint('mismatch_len_original', 'mismatch_len_dedup', run=locals())

    logging.info(f"Writing to {out_flpth}")

    with open(out_flpth, 'w') as fp:
        fp.write('\n'.join(html_l))

