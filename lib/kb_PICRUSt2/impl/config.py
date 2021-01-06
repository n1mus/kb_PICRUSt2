from dotmap import DotMap
import pandas as pd


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)


_config = DotMap({
    'debug': True, # toggle for app-global debug behavior

#-------- my file stuff -------------------------------------------------------------------------------    

    'amplicon_header_name': 'Amplicon_Id', # amplicon index/header name for tables and TSV

#-------- PICRUSt2 files ---------------------------------------------------------------------------

    'metacyc_pathway_code2desc_tsvgz': 
        '/miniconda/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 
    'picrust2_pipeline_flpth': '/miniconda/envs/picrust2/bin/picrust2_pipeline.py',

#-------- these file names/paths should all be in corresponding order ------------------------------

    'tsvgz_relflpth_l': [
        'KO_predicted.tsv.gz',
        'EC_predicted.tsv.gz',
        'pathways_out/path_abun_predictions.tsv.gz',
        'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz',
        'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz',
        'pathways_out/path_abun_unstrat.tsv.gz',
    ],
    'axis_labels': [
        ('amplicon ID', 'KO'), # amplicon x func
        ('amplicon ID', 'EC'), # amplicon x func
        ('amplicon ID', 'MetaCyc pathway'), # amplicon x func
        ('KO', 'sample ID'), # func x sample
        ('EC', 'sample ID'), # func x sample
        ('MetaCyc pathway', 'sample ID'), # func x sample
    ],

    'id_l': [ 
        'amplicon_ko',
        'amplicon_ec',
        'amplicon_metacyc',
        'metagenome_ko',
        'metagenome_ec',
        'metagenome_metacyc',
    ]

})

Var = DotMap(_config) # app-wide globals container

def reset_Var():
    Var.clear()
    Var.update(_config)





'''
TSVs are:

'pathways_out/path_abun_unstrat.tsv', # most important to workflow
'pathways_out/path_abun_unstrat_per_seq.tsv',
'pathways_out/path_abun_predictions.tsv',
'EC_predicted.tsv', # 100M
'KO_predicted.tsv', # 358M (Ginormo)
'EC_metagenome/pred_metagenome_unstrat.tsv',
'KO_metagenome/pred_metagenome_unstrat.tsv',
'''



