# -*- coding: utf-8 -*-
#BEGIN_HEADER
import time
import logging
import os
import sys
import subprocess
import uuid
import functools
import pandas as pd
import numpy as np
import shutil
import json

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.FunctionalProfileUtilClient import FunctionalProfileUtil
from installed_clients.GenericsAPIClient import GenericsAPI

from .impl.kbase_obj import AmpliconMatrix, AttributeMapping
from .impl import appfile
from .impl.config import Var, reset_Var
from .impl import report
from .impl.params import Params
from .util.debug import dprint
from .util.cli import run_check, gunzip
from .util.file import gunzip_out


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

        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.workspace_url = config['workspace-url']
        self.shared_folder = config['scratch']
       
        
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
        # return Variables are: output
        #BEGIN run_picrust2_pipeline
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

        #
        ##
        ### params, app-globals, directories, etc
        ####
        #####

        logging.info('BEGINNING KB_PICRUST2. params: %s' % str(params))

        params = Params(params)

        dprint('params', run=locals())
        
        reset_Var() # clear all fields but `debug`

        Var.update(
            params=params,
            dfu=DataFileUtil(self.callback_url),
            kbr=KBaseReport(self.callback_url),
            fpu=FunctionalProfileUtil(self.callback_url, service_ver='dev'),
            gapi=GenericsAPI(self.callback_url, service_ver='dev'),
            shared_folder=self.shared_folder,
            run_dir=os.path.join(self.shared_folder, 'run_dir_picrust2_' + str(uuid.uuid4())),
            warnings=[],
            objects_created=[],
        )

        os.mkdir(Var.run_dir) # for this API-method run

        Var.update(
            return_dir=os.path.join(Var.run_dir, 'return'),
        )

        os.mkdir(Var.return_dir) # for return input/output/logs etc.

        if Var.debug:
            with open(os.path.join(Var.run_dir, '#params'), 'w') as fh:
                json.dump(params.params, fh)
    
        # TODO document `run_dir` structure

        #
        ##
        ### obj
        ####
        #####


        # instantiate

        amp_mat = AmpliconMatrix(params['amplicon_matrix_upa']) 
        if 'row_attributemapping_ref' in amp_mat.obj:
            row_attrmap = AttributeMapping(amp_mat.obj['row_attributemapping_ref'], amp_mat)
        else:
            msg = (
                "Input AmpliconMatrix "
                "does not have a row AttributeMapping to assign PICRUSt2 functions to."
            )
            logging.warning(msg)
            Var.warnings.append(msg)



        # validate input data

        amp_mat.validate_amplicon_abundance_data()


        # generate input files
        
        seq_flpth = os.path.join(Var.return_dir, 'study_seqs.fna')
        seq_abundance_table_flpth = os.path.join(Var.return_dir, 'study_seqs.tsv') 

        amp_mat.to_fasta(seq_flpth)
        amp_mat.to_seq_abundance_table(seq_abundance_table_flpth)


        # objs should be app globals
        Var.amp_mat = amp_mat


        #
        ##
        ### args
        ####
        #####
        
        
        Var.out_dir = os.path.join(Var.return_dir, 'PICRUSt2_output')
        log_flpth = os.path.join(Var.return_dir, 'cmd_log.txt')


        cmd_pipeline = ' '.join([
            'set -o pipefail &&',
            'source activate picrust2 &&',
            'picrust2_pipeline.py',
            '-s', seq_flpth,
            '-i', seq_abundance_table_flpth,
            '-o', Var.out_dir,
            '--per_sequence_contrib',
            '-p 4',
            '|& tee', log_flpth,
        ])
        
    

        cmd_description = ' \\\n'.join([
            'cd %s &&' % Var.out_dir,
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

        run_check(cmd_pipeline)
        run_check(cmd_description)



        #
        ##
        ### sanity checks
        ####
        #####

        if Var.debug:
            tsvgz_flpth_l = [
                os.path.join(Var.out_dir, tsvgz_relflpth) 
                for tsvgz_relflpth in list(Var.tsvgz_relflpth_l)
            ]
            for i, tsvgz_flpth in enumerate(tsvgz_flpth_l):
                if i in [3,4,5]:
                    # Check no samples dropped (debug)
                    appfile.check_dropped_sample_ids(tsvgz_flpth, amp_mat)
                elif i in [0,1,2]:
                    # Check dropped amplicons are the unaligned/distant ones (debug)
                    appfile.check_dropped_amplicon_ids(tsvgz_flpth, amp_mat)




        #
        ##
        ### update/save Amplicon workflow objects 
        ####
        #####


        path_abun_predictions_tsv_gz_flpth = os.path.join(
            Var.out_dir, 'pathways_out/path_abun_predictions.tsv.gz') 

        attribute = 'PICRUSt2 MetaCyc Pathway Predictions'
        source = 'kb_PICRUSt2/run_picrust2_pipeline'

        # if row AttributeMapping, 
        # update that and referencing objs
        if amp_mat.row_attrmap_upa is not None: 

            # update row AttributeMapping with traits
            id2attr = appfile.parse_picrust2_traits(path_abun_predictions_tsv_gz_flpth)
            ind, attribute = row_attrmap.add_attribute_slot(attribute, source)
            row_attrmap.map_update_attribute(ind, id2attr)
            row_attrmap_upa_new = row_attrmap.save()

            # update AmpliconMatrix which references row AttributeMapping
            amp_mat.obj['row_attributemapping_ref'] = row_attrmap_upa_new
            amp_mat_upa_new = amp_mat.save(name=params['output_name'])         

            Var.objects_created.extend([
                {
                    'ref': row_attrmap_upa_new, 
                    'description': 'Added attribute `%s`' % attribute,
                }, 
                {
                    'ref': amp_mat_upa_new, 
                    'description': 'Updated row AttributeMapping reference'
                },
            ])

    
       

        #
        ##
        ### prepare dir to give FPU decompressed TSVs
        ####
        #####

        # gunzip TSVs out to another directory
        tsv_dir = os.path.join(Var.run_dir, 'decompressed_tsv')
        tsvgz_flpth_l = [os.path.join(Var.out_dir, relflpth) for relflpth in Var.tsvgz_relflpth_l]
        tsv_flpth_l = gunzip_out(tsvgz_flpth_l, tsv_dir)

        # reorder TSVs
        p = [5,3,4,2,0,1]
        tsv_flpth_l = [tsv_flpth_l[i] for i in p]

        # look at TSVs 
        dprint(
            'ls -lh %s/*' % tsv_dir,
            #'file -i %s/*/*' % tsv_dir, 
            run='cli'
        )


        #
        ##
        ### save FunctionalProfile objects
        ####
        #####

        logging.info(
            'Starting saving FunctionalProfiles if any (%s and %s)'
            % (params.getd('create_amplicon_fps'), params.getd('create_sample_fps'))
        )

        if Var.debug:
            FP_amp_mat_ref = params['amplicon_matrix_upa']  # this makes mocking more flexible in case something makes a fake UPA
        else:
            FP_amp_mat_ref = amp_mat_upa_new # this AmpliconMatrix is new one with new AttributeMapping

        ## Metagenome FPs
        if params.getd('create_sample_fps') is True :

            Var.objects_created.append(dict(
                ref=Var.fpu.import_func_profile(dict(
                    workspace_id=Var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_path_abun_unstrat' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[0],
                    profile_type='mg',
                    profile_category='community',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='Metagenome MetaCyc abundance', 
                ))['func_profile_ref'],
                description='Metagenome MetaCyc abundance',
            ))

            Var.objects_created.append(dict(
                ref=Var.fpu.import_func_profile(dict(
                    workspace_id=Var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_EC_pred_metagenome_unstrat' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[1],
                    profile_type='mg',
                    profile_category='community',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='Metagenome EC abundance',
                ))['func_profile_ref'],
                description='Metagenome EC abundance',
            ))

            Var.objects_created.append(dict(
                ref=Var.fpu.import_func_profile(dict(
                    workspace_id=Var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_KO_pred_metagenome_unstrat' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[2],
                    profile_type='mg',
                    profile_category='community',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='Metagenome KO abundance',
                ))['func_profile_ref'],
                description='Metagenome KO abundance',
            ))

        ## Amplicon FPs ##
        if params.getd('create_amplicon_fps') is True:

            Var.objects_created.append(dict(
                ref=Var.fpu.import_func_profile(dict(
                    workspace_id=Var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_path_abun_predictions' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[3],
                    profile_type='amplicon',
                    profile_category='organism',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='Amplicon MetaCyc abundance',
                ))['func_profile_ref'],
                description='Amplicon MetaCyc abundance',
            ))

            Var.objects_created.append(dict(
                ref=Var.fpu.import_func_profile(dict(
                    workspace_id=Var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_EC_predicted' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[4],
                    profile_type='amplicon',
                    profile_category='organism',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='Amplicon EC abundance',
                ))['func_profile_ref'],
                description='Amplicon EC abundance',
            ))
     
            Var.objects_created.append(dict(
                ref=Var.fpu.import_func_profile(dict(
                    workspace_id=Var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_KO_predicted' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[5],
                    profile_type='amplicon',
                    profile_category='organism',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='Amplicon KO abundance',
                ))['func_profile_ref'],
                description='Amplicon KO abundance',
            ))




        #
        ##
        ### html report w/ heatmaps
        ####
        #####

        logging.info('Beginning report business')

        tsvgz_flpth_l = [
            os.path.join(Var.out_dir, tsvgz_relflpth) 
            for tsvgz_relflpth in list(Var.tsvgz_relflpth_l)
        ]

        ##
        ## report

        Var.report_dir = os.path.join(Var.run_dir, 'report')

        report_html_flpth = report.HTMLReportWriter(
            [cmd_pipeline, cmd_description], 
            tsvgz_flpth_l,
            Var.report_dir
        ).write()

        html_links = [{
            'path': Var.report_dir,
            'name': os.path.basename(report_html_flpth),
        }]



        #
        ##
        ### return files
        ####
        #####


        file_links = [{
                'path': Var.return_dir, 
                'name': 'PICRUSt2_results.zip', 
                'description': 'Input, output, cmd, intermediate files, log'
        }]
       
        params_report = {
            'warnings': Var.warnings,
            'objects_created': Var.objects_created,
            'file_links': file_links,
            'html_links': html_links,
            'direct_html_link_index': 0,
            'report_object_name': 'kb_PICRUSt2_report',
            'workspace_name': params['workspace_name'],
            'html_window_height': report.REPORT_HEIGHT,
        }

        Var.params_report = params_report

        obj = Var.kbr.create_extended_report(params_report)

        output = {
            'report_name': obj['name'],
            'report_ref': obj['ref'],
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
