# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import subprocess
import uuid
import functools

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil

from .util.kbase_obj import AmpliconSet, AmpliconMatrix, AttributeMapping
from .util.dprint import dprint
from .util.varstash import Var

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
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)


        Var.update({
            'debug': True,
            'shared_folder': self.shared_folder,
            'callback_url': self.callback_url,
            'dfu': DataFileUtil(self.callback_url),
            'sub_dir': os.path.join(self.shared_folder, str(uuid.uuid4())),
            'suffix': '_' + str(uuid.uuid4()),
            'warnings': [],
            'objects_created': [],
            'picrust2_pckg_dir': '/miniconda/envs/picrust2/lib/python3.6/site-packages/picrust2',
            #'picrust2_pipeline_flpth': '/miniconda/envs/picrust2/bin/picrust2_pipeline.py',
            #'conda_483_flpth': '/miniconda/envs/picrust2/bin/conda',
            })

        os.mkdir(Var.sub_dir)
        
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
        
        Var.update({
            'ctx': ctx,
            'params': params
            }) 

        dprint('params', run=locals())
    


        #
        ##
        ###
        ####
        #####



        logging.info('Loading AmpliconSet and AmpliconMatrix')

        amp_set = AmpliconSet(params['amplicon_set_upa'])
        amp_mat = AmpliconMatrix(amp_set.amp_mat_upa, amp_set) 
        row_attrmap = AttributeMapping(amp_mat.row_attrmap_upa)



        #
        ##
        ### args
        ####
        #####
        
        
        out_dir = os.path.join(Var.sub_dir, 'PICRUSt2_output')
        log_flpth = os.path.join(Var.sub_dir, 'log.txt')



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
            '| tee', log_flpth,
            ])
        
    

        cmd_description = ' \\\n'.join([
            'source activate picrust2 &&',
            'add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC'
            '                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz',
            '&&',
            'add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO'
            '                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz',
            '&&',
            'add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC'
            '                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz'
            ])


        cmd_debug_stderr = ' '.join([
            'source activate picrust &&',
            'sdfasdfsfd'
            ])

        cmd_debug_stdout = ' '.join([
            'source activate picrust2 &&',
            'ping google.com',
            '|& tee x'
            ])


        #
        ##
        ### run
        ####
        #####

        subprocess.run = functools.partial(
            subprocess.run, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=subprocess.PIPE)


        def check(cmd, completed_proc):
            if completed_proc.returncode != 0:
                raise Exception(
                    f"Command: `{cmd}` exited "
                    f"with non-zero return code: {completed_proc.returncode} "
                    f"and error: {completed_proc.stderr.decode('utf-8')}"
                    )


        """
        logging.info(f"Running cmd {cmd_debug_stdout}")
        completed_proc = subprocess.run(cmd_debug_stdout, cwd='/kb/module/work/tmp')
        check(cmd_debug_stdout, completed_proc)


        logging.info(f"Running cmd {cmd_debug_stderr}")
        completed_proc = subprocess.run(cmd_debug_stderr)
        check(cmd_debug_stderr, completed_proc)
        """

        if not (Var.debug and params.get('skip_run')):
            logging.info(f'Running PICRUSt2 via command `{cmd_pipeline}`')
            
            completed_proc = subprocess.run(cmd_pipeline)
            check(cmd_pipeline, completed_proc)


            logging.info(f'Adding descriptions via command `{cmd_description}`')

            completed_proc = subprocess.run(cmd_description)
            check(cmd_description, completed_proc)



        #
        ##
        ### AttributeMapping
        ####
        #####




        if Var.debug and params.get('skip_run'):
            out_dir = '/kb/module/test/data/PICRUSt2_output'


        path_abun_predictions_tsv_gz_flpth = os.path.join(
            out_dir, 'pathways_out/path_abun_predictions.tsv.gz') 

        func2desc_flpth = os.path.join(
            Var.picrust2_pckg_dir, 'default_files/description_mapfiles/metacyc_pathways_info.txt.gz')
        
        id2traits_d = row_attrmap.parse_picrust2_traits(path_abun_predictions_tsv_gz_flpth, func2desc_flpth)
        row_attrmap.add_attribute(id2traits_d)
        row_attrmap_upa_new = row_attrmap.save()

        amp_mat.update_row_attributemapping_ref(row_attrmap_upa_new)
        amp_mat_upa_new = amp_mat.save()         

        amp_set.update_amplicon_matrix_ref(amp_mat_upa_new)
        amp_set_upa_new = amp_set.save()
        

        Var.objects_created = [row_attrmap_upa_new, amp_mat_upa_new, amp_set_upa_new]

        

        #
        ##
        ### return files
        ####
        #####

        if Var.debug and params.get('skip_retFiles'):
            return row_attrmap_upa_new


        def dir_to_shock(dir_path, name, description):
            '''
            For regular directories or html directories
            
            name - for regular directories: the name of the flat (zip) file returned to ui
                   for html directories: the name of the html file
            '''
            dfu_fileToShock_ret = Var.dfu.file_to_shock({
                'file_path': dir_path,
                'make_handle': 0,
                'pack': 'zip',
                })

            dir_shockInfo = {
                'shock_id': dfu_fileToShock_ret['shock_id'],
                'name': name,
                'description': description
                }

            return dir_shockInfo


        shockInfo_retFiles = dir_to_shock(
            Var.sub_dir, 
            'picrust2_results.zip',
            'Input for and output generated by PICRUSt2'
            )




        #
        ##
        ### report
        ####
        #####
        
        
        params_report = {
            'warnings': Var.warnings,
            'file_links': [shockInfo_retFiles],
            'report_object_name': 'kb_PICRUSt2_report',
            'workspace_name': params['workspace_name']
            }

        kbr = KBaseReport(self.callback_url)
        report = kbr.create_extended_report(params_report)

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