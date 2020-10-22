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
import seaborn as sns
import matplotlib.pyplot as plt
import cProfile
import traceback
import retrying

from .config import var
from .dprint import dprint

t0 = None
MAX_LEN = 2800 # 3000^2 is largest size supported by kaleido's chromium JSON stdin TODO adjust to < 9M? <8M?



####################################################################################################
####################################################################################################
def do_heatmap(tsv_flpth, png_flpth, html_flpth): # TODO log coloring?
    '''
    tsv_flpth: data to heatmap. it is a TSV GZ in PICRUSt2's var.out_dir
    png_flpth: where to write png heatmap
    html_flpth: where to write plotly interactive html
    cluster: scipy clustering
    png_from: which package to use to write to png

    plotly can generate the interactive html in good time, but for some reason can take
    all night to write to static image. this reason is that orca is a server contacted through
    a port. using kaleido instead limits the size of the df to write though since it sends
    a json through stdin to a chromium server and that limits the JSON
    '''
    global t0

    df = pd.read_csv(tsv_flpth, sep='\t', index_col=0)
    tsv_flnm = os.path.basename(tsv_flpth)

    dprint('df.shape', run=locals())

    ###
    ###
    if df.shape[0] > MAX_LEN or df.shape[1] > MAX_LEN:
        row_ordering = df.sum(axis=1).values.argsort()[::-1]
        col_ordering = df.sum(axis=0).values.argsort()[::-1]
    
        df = df.iloc[row_ordering, col_ordering]
        df = df.iloc[:MAX_LEN,:MAX_LEN]


    dprint('df.shape', run=locals())

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
        colorbar={'len': 0.3},
    ))
    t_go_heatmap = time.time() - t0

    dprint('t_go_heatmap', run=locals())

    fig.update_layout(
        #height=1000,
        #width=1500,
        title_text=os.path.basename(tsv_flpth) + '<br>shape=' + str(df.shape),
        title_x=0.5,
        xaxis_title=var.tsvTsvgzFlnm2AxisLabels.get(tsv_flnm, (('no axis labels for this flnm',)*2))[1],
        yaxis_title=var.tsvTsvgzFlnm2AxisLabels.get(tsv_flnm, (('no axis labels for this flnm',)*2))[0],
        xaxis_tickangle=45,
        #yaxis_tickangle=45,
    )

    ###
    ###
    logging.info('Writing plotly interactive heatmap at %s' % html_flpth)

    t0 = time.time()
    fig.write_html(
        html_flpth,
        default_height='100%',
        default_width='100%',
    )
    t_write_html = time.time() - t0

    dprint('t_write_html', run=locals())
    
    ###
    ###
    logging.info('Writing plotly static heatmap at %s' % png_flpth)
    t0 = time.time()

    fig.write_image(
        png_flpth, 
        height=700,
        width=1100,
        engine='kaleido',
    )
    
    t_write_image = time.time() - t0

    dprint('t_write_image', run=locals())
    

        

####################################################################################################
####################################################################################################
class HTMLReportWriter:

    def __init__(self, cmd_l, tsv_flpth_l, report_dir): 
        '''
        Input:
        * cmd_l - list of app CLI commands
        * tsv_flpth_l
        * report_dir - not by app-global for some reason
        '''
        self.replacement_d = {}

        #
        self.tsv_flpth_l = tsv_flpth_l
        self.cmd_l = cmd_l

        #
        if not os.path.exists(report_dir):
            os.mkdir(report_dir)

        self.report_dir = report_dir


    def _compile_cmd(self):
        
        txt = ''

        for cmd in self.cmd_l:
            txt += '<pre><code>' + cmd + '</code></pre>\n'

        self.replacement_d['CMD_TAG'] = txt


    def _compile_figures(self):

        logging.info('Compiling all report figures')

        #
        self.fig_dir = os.path.join(self.report_dir, 'fig')
        os.mkdir(self.fig_dir)

        #
        png_flpth_l = [os.path.join(self.fig_dir, os.path.basename(tsv_flpth)[:-4] + '.png') 
            for tsv_flpth in self.tsv_flpth_l]
        html_flpth_l = [os.path.join(self.report_dir, os.path.basename(tsv_flpth)[:-4] + '.html') 
            for tsv_flpth in self.tsv_flpth_l]
        
        for flpth_tup in zip(self.tsv_flpth_l, png_flpth_l, html_flpth_l):
            try:
                do_heatmap(*flpth_tup)
            except Exception as e:
                traceback.print_exc()
                logging.info(
                    'Error occurred when generating heatmap for `%s`. '
                    'Error is:\n`%s`\n' 
                    'Time from last checkpoint is %.2fs\n'
                    'Aborting heatmapping' 
                    % (flpth_tup[0], traceback.format_exc(), time.time() - t0)
                )
                var.warnings.append('Error occurred heatmapping %s. Aborting heatmaps' % flpth_tup[0])
                self.replacement_d['FIGURES_TAG'] = (
                    '<p><i>'
                    'Sorry, an error occurred generating a heatmap for %s. '
                    'See returned files to access the TSVs'
                    '</i></p>' 
                    % os.path.basename(flpth_tup[0])
                )
                return

        def get_relative_fig_path(flpth):
            '''fig path relative to report.html'''
            return '/'.join(flpth.split('/')[-2:])

        # build replacement string
        txt = '<p><i>Heatmaps are restricted by top %d vector sums per dimension. See returned files for full TSVs.</i></p>' % MAX_LEN 
        txt += '<div id="imgLink">\n'
        for png_flpth, html_flpth in zip(png_flpth_l, html_flpth_l):
            txt += '<p><a href="%s" target="_blank"><img alt="%s" src="%s" title="Open to interact"></a></p>\n' % (
                os.path.basename(html_flpth),
                os.path.basename(png_flpth),
                get_relative_fig_path(png_flpth)
            )

        txt += '</div>\n'

        self.replacement_d['FIGURES_TAG'] = txt


    def write(self):
        self._compile_cmd()
        self._compile_figures() # TODO stress test heatmaps

        
        REPORT_HTML_TEMPLATE_FLPTH = '/kb/module/ui/output/report.html'
        html_flpth = os.path.join(self.report_dir, 'report.html')

        with open(REPORT_HTML_TEMPLATE_FLPTH, 'r') as src_fp:
            with open(html_flpth, 'w') as dst_fp:
                for line in src_fp:
                    if line.strip() in self.replacement_d:
                        dst_fp.write(self.replacement_d[line.strip()])
                    else:
                        dst_fp.write(line)
        

        return html_flpth







