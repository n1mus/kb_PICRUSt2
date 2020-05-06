# -*- coding: utf-8 -*-
#BEGIN_HEADER
import time
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
from .util.error import *
from .util.dprint import dprint
from .util.config import _globals, reset
from .util.report import HTMLReportWriter

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


    def run_picrust2_pipeline(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_picrust2_pipeline

        dprint('params', run=locals())
        

        reset(_globals) # clear all fields but `debug`

        _globals.update({
            **self._globals,
            'params': params,
            'run_dir': os.path.join(self._globals['shared_folder'], str(uuid.uuid4())),
            'warnings': [],
            })

        os.mkdir(_globals.run_dir)
    


        #
        ##
        ###
        ####
        #####


        if params.get('skip_obj'):
            logging.info('Skip obj')

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
            if amp_mat.row_attrmap_upa:
                row_attrmap = AttributeMapping(amp_mat.row_attrmap_upa)
            else:
                msg = \
"Input AmpliconSet's associated AmpliconMatrix does not have associated AttributeMapping to assign traits to"
                _globals.warnings.append(msg)



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
                raise NonZeroReturnException(
                    f"Command: `{cmd}` exited "
                    f"with non-zero return code: `{completed_proc.returncode}`. "
                    f"Check logs for details"
                    )


        if params.get('skip_run'):
            logging.info('Skip run')
            out_dir = '/kb/module/test/data/PICRUSt2_output'
            dprint('touch %s' % os.path.join(_globals.run_dir, 'workaround'), run='cli')

        else:

            run_check(cmd_pipeline)
            run_check(cmd_description)


        #
        ##
        ### save objects 
        ####
        #####


        objects_created = []


        if params.get('skip_obj'):
            logging.info('Skip obj')


        else:

            path_abun_predictions_tsv_gz_flpth = os.path.join(
                out_dir, 'pathways_out/path_abun_predictions.tsv.gz') 

            attribute = 'PICRUSt2 MetaCyc Predictions'
            source = 'kb_PICRUSt2'

            if amp_mat.row_attrmap_upa: # if there's an AttributeMapping, update that and referencing chain
                id2traits_d = row_attrmap.parse_picrust2_traits(path_abun_predictions_tsv_gz_flpth)
                row_attrmap.add_attribute_slot(attribute)
                row_attrmap.update_attribute(id2traits_d, attribute, source)
                row_attrmap_upa_new = row_attrmap.save()

                amp_mat.update_row_attributemapping_ref(row_attrmap_upa_new)
                amp_mat_upa_new = amp_mat.save()         

                amp_set.update_amplicon_matrix_ref(amp_mat_upa_new)
                amp_set_upa_new = amp_set.save(name=params.get('output_name'))
                

                objects_created = [
                    {'ref': row_attrmap_upa_new, 'description': 'Added or updated attribute `%s`' % attribute}, 
                    {'ref': amp_mat_upa_new, 'description': 'Updated row AttributeMapping reference'},
                    {'ref': amp_set_upa_new, 'description': 'Updated AmpliconMatrix reference'},
                    ]        



        #
        ##
        ### heatmap html report
        ####
        #####

        html_links = []


        if params.get('skip_report'):
            logging.info('Skip report')

        else:

            tsv_flnm_l = [
                'pathways_out/path_abun_unstrat.tsv', #*
                #'pathways_out/path_abun_unstrat_per_seq.tsv',
                #'pathways_out/path_abun_predictions.tsv',
                'EC_predicted.tsv', #* bigger
                'KO_predicted.tsv', #* big
                #'EC_metagenome/pred_metagenome_unstrat.tsv',
                #'KO_metagenome/pred_metagenome_unstrat.tsv',
                ]

            tsvgz_flpth_l = [os.path.join(out_dir, tsv_flnm + '.gz') for tsv_flnm in tsv_flnm_l]
            

            t0 = time.time()
            report_dir, html_flpth = HTMLReportWriter([cmd_pipeline, cmd_description], tsvgz_flpth_l).write()
            t = time.time() - t0

            dprint('Done with heatmaps. Took %.2fmin' % (t/60))

            html_links.append({
                'path': report_dir,
                'name': os.path.basename(html_flpth),
                })



        #
        ##
        ### return files
        ####
        #####


        if params.get('skip_retFiles'):
            return

        file_links = [
            {
                'path': _globals.run_dir, 
                'name': 'PICRUSt2_results.zip', 
                'description': 'Input, output, intermediate files, logs'
            }
        ]

       
        
        params_report = {
            'warnings': _globals.warnings,
            'objects_created': objects_created,
            'file_links': file_links,
            'html_links': html_links,
            'direct_html_link_index': 0,
            'report_object_name': 'kb_PICRUSt2_report',
            'workspace_name': params['workspace_name'],
            }

        report = _globals.kbr.create_extended_report(params_report)

        output = {
            'report_name': report['name'],
            'report_ref': report['ref'],
        }

        #END run_picrust2_pipeline

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_picrust2_pipeline return value ' +
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
