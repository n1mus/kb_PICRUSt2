import uuid
import logging
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import os
import sys
import time
import plotly
import plotly.graph_objects as go
import traceback
import retrying
import itertools

from .config import Var
from ..util.debug import dprint

t0 = None
REPORT_HEIGHT = 800 # px
MAX_DYN_LEN = 500
#MAX_DYN_SIZE = MAX_DYN_LEN ** 2

'''
Max matrix dim lengths, (sometimes assuming squarishness as upper bound):
* Clustering: >3000
* To use kaleido: ~3000. 3000^2 is largest total size supported by kaleido's chromium JSON stdin.
* To render in narrative: as little as possible. 600^2?
'''

# TODO a element all: inhert; ?
# TODO test things like empty df, large ..
# TODO buttons cut off when wrap 2nd+ row
# TODO Cmd - wrap text after a certain width
# TODO long y-labels -> label and title cut off

####################################################################################################
####################################################################################################
def do_heatmap(tsv_flpth, html_flpth): # TODO log coloring for func x sample?
    '''
    tsv_flpth: data to heatmap. it is a TSV GZ in PICRUSt2's Var.out_dir
    html_flpth: where to write plotly interactive html
    cluster: scipy clustering
    '''
    global t0

    df = pd.read_csv(tsv_flpth, sep='\t', index_col=0) # default infer compression from file name extension
    tsv_flnm = os.path.basename(tsv_flpth)

    dprint('df.shape # original', run=locals())

    ###
    ### subset
    original_shape = df.shape
    subset = False
    if df.shape[0] > MAX_DYN_LEN or df.shape[1] > MAX_DYN_LEN:
        subset = True

        row_ordering = df.sum(axis=1).values.argsort()[::-1]
        col_ordering = df.sum(axis=0).values.argsort()[::-1]
    
        df = df.iloc[row_ordering, col_ordering]
        df = df.iloc[:MAX_DYN_LEN,:MAX_DYN_LEN]


    dprint('df.shape # subset', run=locals())

    ###
    ###
    logging.info('Clustering heatmap for %s' % os.path.basename(tsv_flpth))

    t0 = time.time()
    row_ordering = leaves_list(linkage(df))
    t_cluster_row = time.time() - t0
    t0 = time.time()
    col_ordering = leaves_list(linkage(df.T))
    t_cluster_col = time.time() - t0

    dprint('t_cluster_row', 't_cluster_col', run=locals())

    df = df.iloc[row_ordering, col_ordering]

    ###
    ###
    logging.info('Generating plotly interactive heatmap for %s' % os.path.basename(tsv_flpth))

    t0 = time.time()
    fig = go.Figure(go.Heatmap(
        z=df.values,
        y=df.index.tolist(),
        x=df.columns.tolist(),
        #colorbar={'len': 0.3},
    ))
    t_go_heatmap = time.time() - t0

    dprint('t_go_heatmap', run=locals())

    fig.update_layout(
        #height=1000,
        #width=1500,
        #title_text=os.path.basename(tsv_flpth) + '<br>shape=' + str(df.shape),
        title=dict(
            text=(
                os.path.basename(tsv_flpth) + '<br>' + 
                'shape%s=%s' % (
                    ('<sub>unsubset</sub>' if subset else ''),
                    original_shape,
                )
            ),
            x=0.5,
        ),
        xaxis_title=Var.tsvTsvgzFlnm2AxisLabels.get(tsv_flnm, (('test - no axis labels for this flnm',)*2))[1],
        yaxis_title=Var.tsvTsvgzFlnm2AxisLabels.get(tsv_flnm, (('test - no axis labels for this flnm',)*2))[0],
        xaxis_tickangle=45,
        #yaxis_tickangle=45,
    )

    ###
    ###
    logging.info('Writing plotly interactive heatmap at %s' % html_flpth)

    t0 = time.time()
    fig.write_html(
        html_flpth,
        #default_height='100%',
        #default_width='100%',
    )
    t_write_html = time.time() - t0

    dprint('t_write_html', run=locals())
    
    

        

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
class HTMLReportWriter:

    # incoming TSVs should match this order so can iteratively generate tabcontent divs with right ids
    # don't go by TSV file names from PCRSt2 because those can be duplicates
    fig_id_l = [ 
        'amplicon_ec',
        'amplicon_ko',
        'amplicon_metacyc',
        'community_ec',
        'community_ko',
        'community_metacyc',
    ]

####################################################################################################
####################################################################################################
    def __init__(self, cmd_l, tsv_flpth_l, report_dir): 
        '''
        Input:
        * cmd_l - list of app CLI commands
        * tsv_flpth_l
        * report_dir - not by app-global for some reason

        Sort of self contained, feed it these ^
        '''
        self.replacement_d = {}

        # cycle tsv_flpth_l for testing
        n = len(self.fig_id_l)
        if len(tsv_flpth_l) < n:
            c = itertools.cycle(tsv_flpth_l)
            tsv_flpth_l = [
                next(c) for i in range(n)
            ]

        #
        self.tsv_flpth_l = tsv_flpth_l
        self.cmd_l = cmd_l

        #
        if not os.path.exists(report_dir):
            os.mkdir(report_dir)

        self.report_dir = report_dir


####################################################################################################
####################################################################################################
    def _compile_cmd(self):
        
        txt = ''

        for cmd in self.cmd_l:
            txt += (
                '<p>\n'
                '<code>' + cmd.replace('\n','<br>') + '</code>\n'
                '</p>\n'
            )

        self.replacement_d['CMD_TAG'] = txt


####################################################################################################
####################################################################################################
    def _compile_figures(self):

        logging.info('Compiling report figures')

        #
        html_flpth_l = [
            os.path.join(self.report_dir, fig_id + '.html') 
            for fig_id in self.fig_id_l
        ]
       
        #
        for flpth_tup in zip(self.tsv_flpth_l, html_flpth_l):
            do_heatmap(*flpth_tup)


        # tabcontent divs with iframes
        txt = ''
        for fig_id, html_flpth in zip(self.fig_id_l, html_flpth_l):
            txt += (
                '<div id="%s" class="tabcontent" %s>\n'
                '<iframe src="%s" '
                'scrolling="no" seamless="seamless"'
                '></iframe>\n'
                '</div>\n\n'
                % (
                    fig_id,
                    ('style="display:inline-flex;"' if fig_id == 'community_metacyc' else ''),
                    os.path.basename(html_flpth),
                )
            )


        self.replacement_d['HEATMAPS_TAG'] = txt


####################################################################################################
####################################################################################################
    def write(self):
        self._compile_cmd()
        self._compile_figures() # TODO stress test heatmaps

        
        REPORT_HTML_TEMPLATE_FLPTH = '/kb/module/lib/kb_PICRUSt2/template/report.html'
        html_flpth = os.path.join(self.report_dir, 'report.html')

        with open(REPORT_HTML_TEMPLATE_FLPTH, 'r') as src_fp:
            with open(html_flpth, 'w') as dst_fp:
                for line in src_fp:
                    if line.strip() in self.replacement_d:
                        dst_fp.write(self.replacement_d[line.strip()].strip() + '\n')
                    else:
                        dst_fp.write(line)
        

        return html_flpth







