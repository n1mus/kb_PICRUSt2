from dotmap import DotMap

_config = DotMap({
    'debug': True, # toggle for app-global debug behavior
    'metacyc_pathway_code2desc_tsvgz': 
        '/miniconda/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 
    'picrust2_pipeline_flpth': '/miniconda/envs/picrust2/bin/picrust2_pipeline.py',

    'tsvgz_relflpth_l': [ # TSV GZs relative to var.out_dir, where PICRUSt2 outputs files
        'pathways_out/path_abun_unstrat.tsv.gz', # func x sample
        'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz', # func x sample
        'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', # func x sample
        'pathways_out/path_abun_predictions.tsv.gz', # id x func
        'EC_predicted.tsv.gz', # id x func
        'KO_predicted.tsv.gz', # id x func
    ],
    'tsvgzRelFlpth2axisLabels': { 
            'pathways_out/path_abun_unstrat.tsv.gz': ('MetaCyc pathway', 'sample'), # func x sample
            'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz': ('EC', 'sample'), # func x sample
            'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz': ('KO', 'sample'), # func x sample
            'pathways_out/path_abun_predictions.tsv.gz': ('amplicon ID', 'MetaCyc pathway'), # id x func
            'EC_predicted.tsv.gz': ('amplicon ID', 'EC'), # id x func
            'KO_predicted.tsv.gz': ('amplicon ID', 'KO'), # id x func
    },
    'tsvFlnm2index': {
        'path_abun_unstrat.tsv': 'sample', # func x sample
        'EC_pred_metagenome_unstrat.tsv': 'sample', # func x sample
        'KO_pred_metagenome_unstrat.tsv': 'sample', # func x sample
        'path_abun_predictions.tsv': 'amplicon', # amplicon x func
        'EC_predicted.tsv': 'amplicon', # amplicon x func
        'KO_predicted.tsv': 'amplicon', # amplicon x func
    },

})

var = DotMap(_config) # app-wide globals container

def reset_var():
    var.clear()
    var.update(_config)
