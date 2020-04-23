# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import subprocess
import uuid
import functools
from dotmap import DotMap

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil

from .util.kbase_obj import AmpliconSet, AmpliconMatrix, AttributeMapping
from .util.dprint import dprint
from .util.config import _globals, reset
from .util import report

#END_HEADER


class kb_PICRUSt2:
    '''
    Module Name:
    kb_PICRUSt2

    Module Description:
    A KBase module: kb_PICRUSt2
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        callback_url = os.environ['SDK_CALLBACK_URL']
        workspace_url = config['workspace-url']
        shared_folder = config['scratch']
        
        self._globals = { # shared by all API-method runs (?)
            'callback_url': callback_url,
            'shared_folder': config['scratch'], 
            'dfu': DataFileUtil(callback_url),
            'kbr': KBaseReport(callback_url),
            'picrust2_pckg_dir': '/miniconda/envs/picrust2/lib/python3.6/site-packages/picrust2',
            #'picrust2_pipeline_flpth': '/miniconda/envs/picrust2/bin/picrust2_pipeline.py',
            #'conda_483_flpth': '/miniconda/envs/picrust2/bin/conda',
            }
        
        #END_CONSTRUCTOR
        pass


    def run_kb_PICRUSt2(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_PICRUSt2

        dprint('params', run=locals())
        

        reset(_globals) # clear all fields but `debug`

        _globals.update({
            **self._globals,
            'params': params,
            'run_dir': os.path.join(self._globals['shared_folder'], str(uuid.uuid4())),
            'warnings': [],
            'objects_created': [],
            })

        os.mkdir(_globals.run_dir)
    


        #
        ##
        ###
        ####
        #####


        if params.get('skip_obj'):
            logging.info('Skipping obj')

            dummy = DotMap(
                seq_flpth = 'dummy',
                seq_abundance_table_flpth = 'dummy'
            )

            amp_set = dummy
            amp_mat = dummy

        else:
            
            logging.info('Loading AmpliconSet and AmpliconMatrix')

            amp_set = AmpliconSet(params['amplicon_set_upa'], mini_test=params.get('mini_test'))
            amp_mat = AmpliconMatrix(amp_set.amp_mat_upa) 
            row_attrmap = AttributeMapping(amp_mat.row_attrmap_upa)



        #
        ##
        ### args
        ####
        #####
        
        
        out_dir = os.path.join(_globals.run_dir, 'PICRUSt2_output')
        log_flpth = os.path.join(_globals.run_dir, 'log.txt')


        cmd_pipeline = ' '.join([
            'set -o pipefail &&',
            'source activate picrust2 &&',
            'picrust2_pipeline.py',
            '-s', amp_set.seq_flpth,
            '-i', amp_mat.seq_abundance_table_flpth,
            '-o', out_dir,
            '--per_sequence_contrib',
            '-p 6',
            '--verbose',
            '|& tee', log_flpth,
            ])
        
    

        cmd_description = ' \\\n'.join([
            'cd %s &&' % out_dir,
            'source activate picrust2 &&',
            'add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC',
            '                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz',
            '&&',
            'add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO',
            '                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz',
            '&&',
            'add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC',
            '                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz'
            ])



        with open(os.path.join(_globals.run_dir, 'cmd.txt'), 'w') as fp:
            fp.write(cmd_pipeline + '\n\n')
            fp.write(cmd_description)


        #
        ##
        ### run
        ####
        #####

        subproc_run = functools.partial(
            subprocess.run, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stderr)


        def run_check(cmd):
            logging.info('Running command `%s`' % cmd)
            completed_proc = subproc_run(cmd)
            if completed_proc.returncode != 0:
                raise Exception(
                    f"Command: `{cmd}` exited "
                    f"with non-zero return code: `{completed_proc.returncode}`. "
                    f"Check logs for details"
                    )


        if params.get('skip_run'):
            logging.info('Skip run')
            out_dir = '/kb/module/test/data/PICRUSt2_output'

        else:

            run_check(cmd_pipeline)
            run_check(cmd_description)


        #
        ##
        ### AttributeMapping
        ####
        #####


        if params.get('skip_obj'):
            logging.info('Skip obj')


        else:

            path_abun_predictions_tsv_gz_flpth = os.path.join(
                out_dir, 'pathways_out/path_abun_predictions.tsv.gz') 

            attribute = 'PICRUSt2 MetaCyc Predictions'
            source = 'kb_PICRUSt2'

            id2traits_d = row_attrmap.parse_picrust2_traits(path_abun_predictions_tsv_gz_flpth)
            row_attrmap.add_attribute_slot(attribute)
            row_attrmap.update_attribute(id2traits_d, attribute, source)
            row_attrmap_upa_new = row_attrmap.save()

            amp_mat.update_row_attributemapping_ref(row_attrmap_upa_new)
            amp_mat_upa_new = amp_mat.save()         

            amp_set.update_amplicon_matrix_ref(amp_mat_upa_new)
            amp_set_upa_new = amp_set.save(name=params.get('output_name'))
            

            _globals.objects_created = [
                {'ref': row_attrmap_upa_new, 'description': 'Added or updated attributes for `%s`' % attribute}, 
                {'ref': amp_mat_upa_new, 'description': 'Updated row AttributeMapping reference'},
                {'ref': amp_set_upa_new, 'description': 'Updated AmpliconMatrix reference'},
                ]        



        #
        ##
        ### html report
        ####
        #####

        tsv_flnm_l = [
            #'EC_predicted.tsv',
            #'KO_predicted.tsv',
            #'EC_metagenome/pred_metagenome_unstrat.tsv',
            #'KO_metagenome/pred_metagenome_unstrat.tsv',
            #'pathways_out/path_abun_predictions.tsv',
            #'pathways_out/path_abun_unstrat_per_seq.tsv',
            'pathways_out/path_abun_unstrat.tsv',
            ]

        tsv_dir = os.path.join(_globals.shared_folder, str(uuid.uuid4()))
        os.mkdir(tsv_dir)


        for tsv_flnm in tsv_flnm_l:
            tsv_flpth = os.path.join(out_dir, tsv_flnm)
            subproc_run(f'gunzip -k {tsv_flpth}.gz  && mv {tsv_flpth} {tsv_dir}')
            

        tsv_flpth_l = [os.path.join(tsv_dir, tsv_flnm) for tsv_flnm in os.listdir(tsv_dir)]

        dprint('tsv_flpth_l', run=locals())
        
        global report

        html_dir_l = report.get_html_dir_l(tsv_flpth_l)

        html_flnm_l = []

        for html_dir in html_dir_l:
            for flnm in os.listdir(html_dir):
                if flnm.endswith('.html'):
                    html_flnm_l.append(flnm)






        #
        ##
        ### return files
        ####
        #####

        if params.get('skip_retFiles'):
            return


        def dir_to_shock(dir_path, name, description):
            '''
            For regular directories or html directories
            
            name - for regular directories: the name of the flat (zip) file returned to ui
                   for html directories: the name of the html file
            '''
            dfu_fileToShock_ret = _globals.dfu.file_to_shock({
                'file_path': dir_path,
                'make_handle': 0,
                'pack': 'zip',
                })

            return {
                'shock_id': dfu_fileToShock_ret['shock_id'],
                'name': name,
                'description': description
                }



        shockInfo_retFiles = dir_to_shock(
            _globals.run_dir, 
            'picrust2_results.zip',
            'Input for and output generated by PICRUSt2'
            )

        shockInfo_htmlDir_l = []
        for html_flnm, html_dir in zip(html_flnm_l, html_dir_l):
            shockInfo_htmlDir_l.append(dir_to_shock(
                html_dir,
                html_flnm,
                'A heatmap'
                ))


        #
        ##
        ###
        ####
        #####
        
        
        params_report = {
            'warnings': _globals.warnings,
            'file_links': [shockInfo_retFiles],
            'html_links': shockInfo_htmlDir_l,
            'direct_html_link_index': 0,
            'report_object_name': 'kb_PICRUSt2_report',
            'workspace_name': params['workspace_name'],
            'objects_created': _globals.objects_created,
            }

        report = _globals.kbr.create_extended_report(params_report)

        output = {
            'report_name': report['name'],
            'report_ref': report['ref'],
        }

        #END run_kb_PICRUSt2

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_PICRUSt2 return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
