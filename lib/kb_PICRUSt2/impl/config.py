from dotmap import DotMap # TODO make so fails when accessing something non-existent
import pandas as pd


pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 20)


_config = DotMap(
    debug=True, # toggle for app-global debug behavior

#-------- my file stuff -------------------------------------------------------------------------------    

    amplicon_header_name='Amplicon_Id', # amplicon index/header name for tables and TSV

#-------- PICRUSt2 files ---------------------------------------------------------------------------

    metacyc_pathway_code2desc_tsvgz=
        '/miniconda/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 
    picrust2_pipeline_flpth='/miniconda/envs/picrust2/bin/picrust2_pipeline.py',

#-------- these file names/paths should all be in corresponding order ------------------------------

    func_l=[ # controls order in FP creation and TSV viz
        'cog',
        'ec',
        'ko',
        'pfam',
        'tigrfam',
        'pheno',
        'metacyc',
    ],

    func_2_cfg=dict( 
        cog=dict(
            name='COG',
            title='COG',
            relfp=[
                'COG_predicted.tsv.gz',
                'COG_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            ],
        ),      

        ec=dict(
            name='EC',
            title='EC',
            relfp=[
                'EC_predicted.tsv.gz',
                'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            ],
        ),

        ko=dict(
            name='KO',
            title='KO',
            relfp=[
                'KO_predicted.tsv.gz',
                'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            ],
        ),

        pfam=dict(
            name='Pfam',
            title='Pfam',
            relfp=[
                'PFAM_predicted.tsv.gz',
                'PFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            ],
        ),

        tigrfam=dict(
            name='TIGRFAMs',
            title='TIGRFAMs',
            relfp=[
                'TIGRFAM_predicted.tsv.gz',
                'TIGRFAM_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            ],
        ),

        pheno=dict(
            name='IMG phenotype',
            title='IMG Phenotype',
            relfp=[
                'PHENO_predicted.tsv.gz',
                'PHENO_metagenome_out/pred_metagenome_unstrat.tsv.gz',
            ],
        ),

        metacyc=dict(
            name='MetaCyc',
            title='MetaCyc',
            relfp=[
                'pathways_out/path_abun_predictions.tsv.gz',
                'pathways_out/path_abun_unstrat.tsv.gz',
            ],
  
        ),



    ),

)

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
'EC_metagenome_out/pred_metagenome_unstrat.tsv',
'KO_metagenome_out/pred_metagenome_unstrat.tsv',






    'id_l': [ # this one holds ordering info
        'amplicon_ko',
        'amplicon_ec',
        'amplicon_metacyc',
        #'amplicon_cog',
        #'amplicon_pfam',
        #'amplicon_tigrfam',
        #'amplicon_pheno',
        'metagenome_ko',
        'metagenome_ec',
        'metagenome_metacyc',
        #'metagenome_cog',
        #'metagenome_pfam',
        #'metagenome_tigrfam',
        #'metagenome_pheno',
    ],

    'optional_l' = [
        'cog',
        'pfam',
        'tigrfam',
        'pheno',
    ],

    'tsvgz_relflpth_l': [
        'KO_predicted.tsv.gz',
        'EC_predicted.tsv.gz',
        'pathways_out/path_abun_predictions.tsv.gz',
        'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz',
        'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz',
        'pathways_out/path_abun_unstrat.tsv.gz',
    ],
    'axis_labels': [
        ('amplicon ID', 'KO'),
        ('amplicon ID', 'EC'),
        ('amplicon ID', 'MetaCyc pathway'),
        ('KO', 'sample ID'),
        ('EC', 'sample ID'),
        ('MetaCyc pathway', 'sample ID'),
    ],



    'id_2_axis_labels': {
        'amplicon_ko': ('amplicon ID', 'KO'),
        'amplicon_ec': ('amplicon ID', 'EC'),
        'amplicon_metacyc': ('amplicon ID', 'MetaCyc pathway'),
        'metagenome_ko': ('KO', 'sample ID'),
        'metagenome_ec': ('EC', 'sample ID'),
        'metagenome_metacyc': ('MetaCyc pathway', 'sample ID'),
    },











'''



