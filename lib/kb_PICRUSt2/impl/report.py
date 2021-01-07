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

REPORT_HEIGHT = 800 # px
MAX_DYN_LEN = 500
#MAX_DYN_SIZE = MAX_DYN_LEN ** 2

'''
Max matrix dim lengths, (sometimes assuming squarishness as upper bound):
* Clustering: >3000
* To use kaleido: ~3000. 3000^2 is largest total size supported by kaleido's chromium JSON stdin.
* To render in narrative: as little as possible. 600^2?
'''

# TODO test things like empty df, large ..
# TODO full function in hover? too much data ...

####################################################################################################
####################################################################################################
def do_heatmap(tsv_flpth, html_flpth, axis_labels): # TODO log coloring for func x sample?
    '''
    tsv_flpth: data to heatmap. it is a TSV GZ in PICRUSt2's Var.out_dir
    html_flpth: where to write plotly html
    '''

    df = pd.read_csv(tsv_flpth, sep='\t', index_col=0) # default infer compression from file name extension
    tsv_flnm = os.path.basename(tsv_flpth)

    ###
    ### subset
    original_shape = df.shape
    subset = df.shape[0] > MAX_DYN_LEN or df.shape[1] > MAX_DYN_LEN
    if subset is True:
        row_ordering = df.sum(axis=1).values.argsort()[::-1]
        col_ordering = df.sum(axis=0).values.argsort()[::-1]
    
        df = df.iloc[row_ordering, col_ordering]
        df = df.iloc[:MAX_DYN_LEN,:MAX_DYN_LEN]

    ###
    ###
    row_ordering = leaves_list(linkage(df))
    col_ordering = leaves_list(linkage(df.T))

    df = df.iloc[row_ordering, col_ordering]

    ###
    ###
    fig = go.Figure(go.Heatmap(
        z=df.values,
        y=df.index.tolist(),
        x=df.columns.tolist(),
    ))

    fig.update_layout(
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
        xaxis_title=axis_labels[1],
        yaxis_title=axis_labels[0],
        xaxis_tickangle=45,
    )

    ###
    ###
    fig.write_html(
        html_flpth,
    )

    
    

        

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
class HTMLReportWriter:


####################################################################################################
####################################################################################################
    def __init__(self, cmd_l, tsv_flpth_l, report_dir): 
        '''
        tsv_flpth_l should have TSVs corresponding to Var.id_l etc.
        '''

        self.replacement_d = {}

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
                '<p class="fixwhitespace">\n'
                '<code>' + cmd + '</code>\n'
                '</p>\n'
            )

        self.replacement_d['CMD_TAG'] = txt


####################################################################################################
####################################################################################################
    def _compile_figures(self):

        #
        html_flpth_l = [
            os.path.join(self.report_dir, fig_id + '.html') 
            for fig_id in Var.id_l
        ]
       
        #
        for flpth_tup in zip(self.tsv_flpth_l, html_flpth_l, Var.axis_labels):
            do_heatmap(*flpth_tup)


        # tabcontent divs with iframes
        txt = ''
        for fig_id, html_flpth in zip(Var.id_l, html_flpth_l):
            txt += (
                '<div id="%s" class="tabcontent" %s>\n'
                '<iframe src="%s" '
                'scrolling="no" seamless="seamless"'
                '></iframe>\n'
                '</div>\n\n'
                % (
                    fig_id,
                    ('style="display:inline-flex;"' if fig_id == 'metagenome_metacyc' else ''),
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







