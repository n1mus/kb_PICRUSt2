from dotmap import DotMap

_config = DotMap({
    'debug': True, # toggle for global debug behavior
    'metacyc_pathway_code2desc_tsvgz': 
        '/miniconda/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz', 
    'picrust2_pipeline_flpth': '/miniconda/envs/picrust2/bin/picrust2_pipeline.py',
    #'picrust2_pckg_dir': '/miniconda/envs/picrust2/lib/python3.6/site-packages/picrust2',
    #'conda_483_flpth': '/miniconda/envs/picrust2/bin/conda',
})

var = DotMap(_config) # app-wide globals container

def reset_var():
    var.clear()
    var.update(_config)
