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
import gzip
import shutil

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.FunctionalProfileUtilClient import FunctionalProfileUtil
from installed_clients.GenericsAPIClient import GenericsAPI

from .util.kbase_obj import AmpliconMatrix, AttributeMapping
from .util.error import *
from .util.dprint import dprint
from .util.config import var, reset_var
from .util.report import HTMLReportWriter


####################################################################################################
####################################################################################################
def run_check(cmd):
    '''Wrap tool-running method for patching'''

    logging.info(f'Running PICRUSt2 via command `{cmd}`')

    t0 = time.time()
    completed_proc = subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=sys.stdout, stderr=sys.stdout)

    logging.info('Completed in %.1f minutes' % ((time.time()-t0)/60))

    if completed_proc.returncode != 0:
        msg = (
            "PICRUSt2 command `%s` returned with non-zero return code `%d`. "
            "Please check logs for more details" % 
            (cmd, completed_proc.returncode)
        )
        raise NonZeroReturnException(msg)

####################################################################################################
####################################################################################################
class OutfileWrangler:

    #####
    #####
    @staticmethod
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

    #####
    #####
    @staticmethod
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

        id_x_desc_df = OutfileWrangler.do_code2desc(id_x_code_df, var.metacyc_pathway_code2desc_tsvgz)

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


    #####
    #####
    @staticmethod
    def pad_0_vecs(tsv_flpth, amp_mat):
        '''
        PICRUSt2 drops ids/samples that are all 0s
        Restore them in original order here
        By 'id' I mean 'amplicon'
        '''
        flnm2index =  {
            'path_abun_unstrat.tsv': 'sample', # func x sample
            'EC_pred_metagenome_unstrat.tsv': 'sample', # func x sample
            'KO_pred_metagenome_unstrat.tsv': 'sample', # func x sample
            'path_abun_predictions.tsv': 'amplicon', # amplicon x func
            'EC_predicted.tsv': 'amplicon', # amplicon x func
            'KO_predicted.tsv': 'amplicon', # amplicon x func
        }

        index = flnm2index[os.path.basename(tsv_flpth)]
        df_partial = pd.read_csv(tsv_flpth, sep='\t', index_col=0)

        # orient PICRUSt2 output matrix as something vs. func
        if index == 'sample':
            df_partial = df_partial.T
            id_l_full = amp_mat.obj['data']['col_ids']
        elif index == 'amplicon':
            id_l_full = amp_mat.obj['data']['row_ids']
        else:
            raise Exception()

        if df_partial.shape[0] == len(id_l_full):
            return

        df_full = pd.DataFrame(
            np.zeros((len(id_l_full), df_partial.shape[1])), 
            index=id_l_full,
            columns=df_partial.columns
        )

        df_full.loc[df_partial.index, df_partial.columns] = df_partial.values

        # undo orient
        if index == 'sample':
            df_full = df_full.T

        df_full.to_csv(tsv_flpth, sep='\t')
        

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
        # return variables are: output
        #BEGIN run_picrust2_pipeline
####################################################################################################
####################################################################################################

        logging.info('params: %s' % str(params))
        

        reset_var() # clear all fields but `debug`

        var.update(
            params=params,
            dfu=DataFileUtil(self.callback_url),
            kbr=KBaseReport(self.callback_url),
            fpu=FunctionalProfileUtil(self.callback_url, service_ver='dev'),
            gapi=GenericsAPI(self.callback_url, service_ver='dev'),
            shared_folder=self.shared_folder,
            run_dir=os.path.join(self.shared_folder, 'run_dir_picrust2_' + str(uuid.uuid4())),
            warnings= [],
        )

        os.mkdir(var.run_dir) # for this API-method run

        var.update(
            return_dir=os.path.join(var.run_dir, 'return'),
        )

        os.mkdir(var.return_dir) # for return input/output/logs etc.
    
        # TODO document `run_dir` structure

        #
        ##
        ### obj
        ####
        #####


        # instantiate

        amp_mat = AmpliconMatrix(params['amplicon_matrix_upa']) 
        if amp_mat.row_attrmap_upa:
            row_attrmap = AttributeMapping(amp_mat.row_attrmap_upa)
        else:
            msg = (
                "Input AmpliconMatrix "
                "does not have a row AttributeMapping to assign PICRUSt2 traits to. "
                "To create one, ??"
            )
            var.warnings.append(msg)



        # generate input files
        
        seq_flpth = os.path.join(var.return_dir, 'study_seqs.fna')
        seq_abundance_table_flpth = os.path.join(var.return_dir, 'study_seqs.tsv') 


        #amp_set.to_fasta(seq_flpth)
        amp_mat.to_fasta(seq_flpth)
        amp_mat.to_seq_abundance_table(seq_abundance_table_flpth)




        #
        ##
        ### args
        ####
        #####
        
        
        var.out_dir = os.path.join(var.return_dir, 'PICRUSt2_output')
        log_flpth = os.path.join(var.return_dir, 'cmd_log.txt')


        cmd_pipeline = ' '.join([
            'set -o pipefail &&',
            'source activate picrust2 &&',
            'picrust2_pipeline.py',
            '-s', seq_flpth,
            '-i', seq_abundance_table_flpth,
            '-o', var.out_dir,
            '--per_sequence_contrib',
            '-p 4',
            '|& tee', log_flpth,
        ])
        
    

        cmd_description = ' \\\n'.join([
            'cd %s &&' % var.out_dir,
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
        ### save Amplicon workflow objects 
        ####
        #####


        path_abun_predictions_tsv_gz_flpth = os.path.join(
            var.out_dir, 'pathways_out/path_abun_predictions.tsv.gz') 

        attribute = 'PICRUSt2 MetaCyc Pathway Predictions'
        source = 'kb_PICRUSt2/run_picrust2_pipeline'

        var.objects_created = []

        # if row AttributeMapping, 
        # update that and referencing objs
        if amp_mat.row_attrmap_upa is not None: 

            # update row AttributeMapping with traits
            id2traits_d = OutfileWrangler.parse_picrust2_traits(path_abun_predictions_tsv_gz_flpth)
            ind = row_attrmap.get_attribute_slot_warn(attribute, source)
            row_attrmap.update_attribute(ind, id2traits_d)
            row_attrmap_upa_new = row_attrmap.save()

            # update AmpliconMatrix which references row AttributeMapping
            amp_mat.update_row_attributemapping_ref(row_attrmap_upa_new)
            amp_mat_upa_new = amp_mat.save(name=params.get('output_name'))         

            var.objects_created.extend([
                {'ref': row_attrmap_upa_new, 'description': 'Added or updated attribute `%s`' % attribute}, 
                {'ref': amp_mat_upa_new, 'description': 'Updated row AttributeMapping reference'},
            ])

            FP_amp_mat_ref = amp_mat_upa_new
        
        else:

            FP_amp_mat_ref = params['amplicon_matrix_upa']
    
       

        #
        ##
        ### save FunctionalProfile objects
        ####
        #####

        def gunzip(read, write):
            '''
            Gunzip from `read` to `write`
            '''
            with gzip.open(read, 'rb') as fh_read:
                with open(write, 'wb') as fh_write:
                    shutil.copyfileobj(fh_read, fh_write)


        tsv_dir = os.path.join(var.shared_folder, 'tsv_dir_kbpicrust2' + str(uuid.uuid4()))
        os.mkdir(tsv_dir)

        # These are PICRUSt2 output tsvgz needed for FPs TODO make this list  and axisLabls mappings app globals
        tsvgz_relflpth_l = [
            'pathways_out/path_abun_unstrat.tsv.gz', # func x sample
            'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz', # func x sample
            'KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', # func x sample
            'pathways_out/path_abun_predictions.tsv.gz', # id x func
            'EC_predicted.tsv.gz', # id x func
            'KO_predicted.tsv.gz', # id x func
        ]

        # These are destination tsvs to give to FPU
        tsv_flpth_l = [os.path.join(tsv_dir, os.path.basename(tsvgz_relflpth[:-3])) for tsvgz_relflpth in tsvgz_relflpth_l]
        tsv_flpth_l[1] = os.path.dirname(tsv_flpth_l[1]) + '/EC_' + os.path.basename(tsv_flpth_l[1])
        tsv_flpth_l[2] = os.path.dirname(tsv_flpth_l[2]) + '/KO_' + os.path.basename(tsv_flpth_l[2])

        # Gunzip/copy
        for tsvgz_relflpth, tsv_flpth in zip(tsvgz_relflpth_l, tsv_flpth_l):
            gunzip(
                os.path.join(var.out_dir, tsvgz_relflpth),
                tsv_flpth
            )

        # Pad 0 ids/samples
        for tsv_flpth in tsv_flpth_l:
            OutfileWrangler.pad_0_vecs(tsv_flpth, amp_mat)

        dprint('tsv_flpth_l', run=locals())




        ## Community FPs ##

        if 'sample_set_ref' not in amp_mat.obj:
            msg = (
                'Sorry, input AmpliconMatrix %s does not have a SampleSet reference '
                'and so creating community FunctionalProfiles is prohibited. '
                'Please see importers to link a SampleSet'
                % amp_mat.name
            )
            logging.warning(msg)
            var.warnings.append(msg)
        
        else:

            var.objects_created.append(dict(
                ref=var.fpu.import_func_profile(dict(
                    workspace_id=var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_path_abun_unstrat' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[0],
                    profile_type='mg',
                    profile_category='community',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='MetaCyc vs. Sample',
                ))['func_profile_ref'],
                description='MetaCyc vs. Sample',
            ))

            var.objects_created.append(dict(
                ref=var.fpu.import_func_profile(dict(
                    workspace_id=var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_EC_pred_metagenome_unstrat' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[1],
                    profile_type='mg',
                    profile_category='community',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='EC vs. sample',
                ))['func_profile_ref'],
                description='EC vs. Sample',
            ))

            var.objects_created.append(dict(
                ref=var.fpu.import_func_profile(dict(
                    workspace_id=var.params['workspace_id'],
                    func_profile_obj_name='%s.PICRUSt2_KO_pred_metagenome_unstrat' % amp_mat.name,
                    original_matrix_ref=FP_amp_mat_ref,
                    profile_file_path=tsv_flpth_l[2],
                    profile_type='mg',
                    profile_category='community',
                    data_epistemology='predicted',
                    epistemology_method='PICRUSt2',
                    description='KO vs. Sample',
                ))['func_profile_ref'],
                description='KO vs. Sample',
            ))

        ## Amplicon FPs ##

        var.objects_created.append(dict(
            ref=var.fpu.import_func_profile(dict(
                workspace_id=var.params['workspace_id'],
                func_profile_obj_name='%s.PICRUSt2_path_abun_predictions' % amp_mat.name,
                original_matrix_ref=FP_amp_mat_ref,
                profile_file_path=tsv_flpth_l[3],
                profile_type='amplicon',
                profile_category='organism',
                data_epistemology='predicted',
                epistemology_method='PICRUSt2',
                description='Amplicon vs. MetaCyc',
            ))['func_profile_ref'],
            description='Amplicon vs. MetaCyc',
        ))

        var.objects_created.append(dict(
            ref=var.fpu.import_func_profile(dict(
                workspace_id=var.params['workspace_id'],
                func_profile_obj_name='%s.PICRUSt2_EC_predicted' % amp_mat.name,
                original_matrix_ref=FP_amp_mat_ref,
                profile_file_path=tsv_flpth_l[4],
                profile_type='amplicon',
                profile_category='organism',
                data_epistemology='predicted',
                epistemology_method='PICRUSt2',
                description='Amplicon vs. EC',
            ))['func_profile_ref'],
            description='Amplicon vs. EC',
        ))
 
        var.objects_created.append(dict(
            ref=var.fpu.import_func_profile(dict(
                workspace_id=var.params['workspace_id'],
                func_profile_obj_name='%s.PICRUSt2_KO_predicted' % amp_mat.name,
                original_matrix_ref=FP_amp_mat_ref,
                profile_file_path=tsv_flpth_l[5],
                profile_type='amplicon',
                profile_category='organism',
                data_epistemology='predicted',
                epistemology_method='PICRUSt2',
                description='Amplicon vs. KO',
            ))['func_profile_ref'],
            description='Amplicon vs. KO',
        ))


        #
        ##
        ### html report w heatmap(s)
        ####
        #####

        ##
        ## input tsv(gz)

        ''' Plotting input tsv tricky because may need log scale
        
        seq_abundance_table_tsvgz_flpth = os.path.join(var.run_dir, os.path.basename(seq_abundance_table_flpth) + '.gz')

        subprocess.run(
                'gzip -c %s > %s' % (seq_abundance_table_flpth, seq_abundance_table_tsvgz_flpth),
                stdout=sys.stdout,
                stderr=sys.stderr,
                shell=True,
                executable='/bin/bash'
        )
        '''

        ##
        ## PICRUSt2 output tsvgzs

        '''
        TSVs are:

        'pathways_out/path_abun_unstrat.tsv', #* [most important] 1.2M
        'pathways_out/path_abun_unstrat_per_seq.tsv',
        'pathways_out/path_abun_predictions.tsv',
        'EC_predicted.tsv', #* [nice-to-have] 100M TODO
        'KO_predicted.tsv', #* [nice-to-have] 358M TODO
        'EC_metagenome/pred_metagenome_unstrat.tsv',
        'KO_metagenome/pred_metagenome_unstrat.tsv',
        '''

        tsv_flnm_l = [
            'pathways_out/path_abun_unstrat.tsv', 
            'EC_predicted.tsv',
            'KO_predicted.tsv',
        ]

        tsvgz_flpth_l = [os.path.join(var.out_dir, tsv_flnm + '.gz') for tsv_flnm in tsv_flnm_l]
        
        ##
        ## report

        var.report_dir = os.path.join(var.run_dir, 'report')

        t0 = time.time()
        html_flpth = HTMLReportWriter(
                [cmd_pipeline, cmd_description], 
                tsvgz_flpth_l, 
                var.report_dir
        ).write()
        t = time.time() - t0

        dprint('Done with heatmaps. Took %.1f min' % (t/60))

        html_links = [{
            'path': var.report_dir,
            'name': os.path.basename(html_flpth),
        }]



        #
        ##
        ### return files
        ####
        #####


        file_links = [{
                'path': var.return_dir, 
                'name': 'PICRUSt2_results.zip', 
                'description': 'Input, output, cmd, intermediate files, log'
        }]
       
        params_report = {
            'warnings': var.warnings,
            'objects_created': var.objects_created,
            'file_links': file_links,
            'html_links': html_links,
            'direct_html_link_index': 0,
            'report_object_name': 'kb_PICRUSt2_report',
            'workspace_name': params['workspace_name'],
        }

        report = var.kbr.create_extended_report(params_report)

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
