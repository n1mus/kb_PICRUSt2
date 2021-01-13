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

'FAPROTAX Functions <taxonomy="RDP Clsf Taxonomy, <gene="ssu", minWords="default", conf=0.8>">'

####################################################################################################
####################################################################################################
def do_heatmap(tsv_fp, html_fp, axis_labels): # TODO log coloring for func x sample?
    '''
    tsv_fp: data to heatmap. it is a TSV GZ in PICRUSt2's Var.out_dir
    html_fp: where to write plotly html
    '''

    df = pd.read_csv(tsv_fp, sep='\t', index_col=0) # default infer compression from file name extension
    tsv_flnm = os.path.basename(tsv_fp)

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
        colorbar=dict(
            title=dict(
                text='Abundance',
                side='right',
            ),
        ),
    ))

    fig.update_layout(
        title=dict(
            text=(
                os.path.basename(tsv_fp) + '<br>' + 
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
        html_fp,
    )

    
    

        

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
class HTMLReportWriter:
    '''
    Needs to know Var.report_dir
    '''


####################################################################################################
####################################################################################################
    def __init__(self, cmd_l): 
        '''
        '''

        self.replacement_d = {}

        #
        self.cmd_l = cmd_l

        if not os.path.exists(Var.report_dir):
            os.mkdir(Var.report_dir)



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



        button_l = []
        content_l = []
        for per in ['amplicon', 'metagenome']:
            for func in Var.func_l:
                if not Var.params.getd(func):
                    continue


                fig_id = per + '_' + func
                func_name = Var.func_2_cfg[func]['name']
                fig_title = per.title() + ' ' + func_name

                ind = 0 if per=='amplicon' else 1
                tsv_fp = os.path.join(Var.out_dir, Var.func_2_cfg[func]['relfp'][ind])
                html_fp = os.path.join(Var.report_dir, fig_id + '.html')

                axis_labels = (
                    ('amplicon', func_name) if per=='amplicon' else
                    (func_name, 'sample')
                )

                do_heatmap(tsv_fp, html_fp, axis_labels)

                button_l.append(
                    '''<button class="tablinks %s" onclick="openTab(event, '%s')">%s</button>'''  
                    % (
                        'active' if fig_id == 'metagenome_metacyc' else '',
                        fig_id, 
                        fig_title,
                    ) 
                )

                content_l.append(
                    '<div id="%s" class="tabcontent" %s>\n' % (
                        fig_id,
                        ('style="display:inline-flex;"' if fig_id == 'metagenome_metacyc' else ''),
                    ) +
                    '<iframe src="%s" scrolling="no" seamless="seamless"></iframe>\n' % os.path.basename(html_fp) +
                    '</div>\n'
                )

        self.replacement_d['HEATMAP_BUTTON_TAG'] = '\n'.join(button_l)
        self.replacement_d['HEATMAP_CONTENT_TAG'] = '\n'.join(content_l)


####################################################################################################
####################################################################################################
    def write(self):
        self._compile_cmd()
        self._compile_figures() # TODO stress test heatmaps

        
        REPORT_HTML_TEMPLATE_FLPTH = '/kb/module/lib/kb_PICRUSt2/template/report.html'
        html_fp = os.path.join(Var.report_dir, 'report.html')

        with open(REPORT_HTML_TEMPLATE_FLPTH, 'r') as src_fh:
            with open(html_fp, 'w') as dst_fh:
                for line in src_fh:
                    s = line.strip()
                    if s in self.replacement_d:
                        dst_fh.write(self.replacement_d[s].strip() + '\n')
                    else:
                        dst_fh.write(line)
        

        return html_fp







