import unittest
from unittest.mock import patch
import uuid
import pandas as pd
import numpy as np
import itertools
import shutil

from kb_PICRUSt2.kb_PICRUSt2Server import MethodContext
from kb_PICRUSt2.authclient import KBaseAuth as _KBaseAuth
from installed_clients.WorkspaceClient import Workspace

from kb_PICRUSt2.kb_PICRUSt2Impl import kb_PICRUSt2 
from kb_PICRUSt2.impl.kbase_obj import AmpliconMatrix, AttributeMapping
from kb_PICRUSt2.impl import appfile
from kb_PICRUSt2.impl.report import do_heatmap, HTMLReportWriter
from kb_PICRUSt2.impl.config import Var
from kb_PICRUSt2.impl.error import * # Exceptions
from kb_PICRUSt2.util.debug import dprint
from kb_PICRUSt2.util.validate import ValidationException
from mock import *
import config


class ReportTest(config.BaseTest):
####################################################################################################
####################################################################################################
    def test_large_heatmap(self):
        '''
        '''
        logging.info('Testing with test_large_heatmap')

        ##
        def write_random_tsv(flpth, dim, max):
            values = (np.random.random((dim, dim)) * max).round(decimals=2)
            df = pd.DataFrame(
                    values, 
                    index=['dummy_ind_%d' % i for i in range(dim)], 
                    columns=['dummy_col_%d' % i for i in range(dim)]
            )
            df.to_csv(flpth, sep='\t', compression='gzip')

        ##
        def has_n_htmls(dir, n):
            return len([flnm for flnm in os.listdir(report_dir) if flnm.endswith('.html')]) == n
    
        ## make `run_dir`
        run_dir = os.path.join(config.shared_folder, 'test_report_' + str(uuid.uuid4()))
        os.mkdir(run_dir)

        ###
        ### heatmap random large
        
        report_dir = os.path.join(run_dir, 'report_randomLargeHeatmap')
        fig_dir = os.path.join(report_dir, 'fig')
        os.makedirs(fig_dir)

        tsvgz_flpth = os.path.join(report_dir, 'random.tsv.gz')
        html_flpth = os.path.join(report_dir, 'heatmap_random.html')

        write_random_tsv(tsvgz_flpth, dim=7000, max=1500)

        do_heatmap(tsvgz_flpth, html_flpth)



####################################################################################################
####################################################################################################
    def test_report(self):
        '''
        Check them (`cd test_local/workdir/tmp && firefox test_report_*/report_dir_*/*.html &)
        '''

        logging.info('Testing with test_report')


        ##
        def write_random_tsv(flpth, dim, max):
            values = (np.random.random((dim, dim)) * max).round(decimals=2)
            df = pd.DataFrame(
                    values, 
                    index=['dummy_ind_%d' % i for i in range(dim)], 
                    columns=['dummy_col_%d' % i for i in range(dim)]
            )
            df.to_csv(flpth, sep='\t', compression='gzip')

        ##
        def has_n_htmls(dir, n):
            return len([flnm for flnm in os.listdir(report_dir) if flnm.endswith('.html')]) == n
    
        ## make `run_dir`
        run_dir = os.path.join(config.shared_folder, 'test_report_' + str(uuid.uuid4()))
        os.mkdir(run_dir)



        ###
        ### heatmap enigma17770by511
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_enigma17770by511')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_enigma17770by511 = '/kb/module/test/data/by_dataset_input/enigma17770by511/return/PICRUSt2_output/pathways_out/path_abun_unstrat.tsv.gz'
            html_flpth = os.path.join(report_dir, 'heatmap_enigma17770by511.html')

            do_heatmap(tsvgz_flpth_enigma17770by511, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### heatmap enigma50by30
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_enigma50by30')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_enigma50by30 = '/kb/module/test/data/by_dataset_input/enigma50by30/return/PICRUSt2_output/pathways_out/path_abun_unstrat.tsv.gz'
            html_flpth = os.path.join(report_dir, 'heatmap_enigma50by30.html')

            do_heatmap(tsvgz_flpth_enigma50by30, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### heatmap random
        with self.subTest():
            
            report_dir = os.path.join(run_dir, 'report_random_2k')
            fig_dir = os.path.join(report_dir, 'fig')
            os.makedirs(fig_dir)

            tsvgz_flpth_random = os.path.join(report_dir, 'random.tsv.gz')
            html_flpth = os.path.join(report_dir, 'heatmap_random.html')

            write_random_tsv(tsvgz_flpth_random, dim=2000, max=100)

            do_heatmap(tsvgz_flpth_random, html_flpth)

            self.assertTrue(has_n_htmls(report_dir, 1))


        ###
        ### HTMLReportWriter with 1 tsvgz
        with self.subTest():

            cmd_l = ['wingardium leviosa']

            tsvgz_flpth_l = [
                tsvgz_flpth_enigma17770by511,
            ]

            report_dir = os.path.join(run_dir, 'report_HTMLReportWriter_1heatmap')

            HTMLReportWriter(cmd_l, tsvgz_flpth_l, report_dir).write()

            self.assertTrue(has_n_htmls(report_dir, 7))


        ###
        ### HTMLReportWriter with 3 tsvgz
        with self.subTest():

            cmd_l = ['expecto patronum', 'accio heatmap', 'luminos']

            report_dir = os.path.join(run_dir, 'report_HTMLReportWriter_3heatmaps')
            os.mkdir(report_dir)

            # copy tsvgzs into `report_dir` with unique names
            # since some have same filenames

            tsvgz_flpth_l = [
                    tsvgz_flpth_enigma17770by511,
                    tsvgz_flpth_enigma50by30,
                    tsvgz_flpth_random,
            ]

            HTMLReportWriter(cmd_l, tsvgz_flpth_l, report_dir).write()

            self.assertTrue(has_n_htmls(report_dir, 7))


