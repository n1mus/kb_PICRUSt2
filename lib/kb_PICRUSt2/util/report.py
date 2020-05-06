import uuid
import logging
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import os
import sys
import time
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt

from .config import _globals
from .dprint import dprint



def do_heatmap(tsvgz_flpth, png_flpth, html_flpth, cluster=False):
    df = pd.read_csv(tsvgz_flpth, sep='\t', index_col=0, compression='gzip')

    if cluster:
        t0 = time.time()
        row_ordering = leaves_list(linkage(df))
        col_ordering = leaves_list(linkage(df.T))
        t_cluster = time.time() - t0

        dprint('t_cluster', run=locals())

        df = df.iloc[row_ordering, col_ordering]

    logging.info('Generating interactive heatmap for %s' % os.path.basename(tsvgz_flpth))

    t0 = time.time()
    fig = go.Figure(go.Heatmap(
        z=df.values,
        y=df.index.tolist(),
        x=df.columns.tolist(),
        colorbar={'len': 0.3}
        ))
    t_go_heatmap = time.time() - t0

    dprint('t_go_heatmap', run=locals())

    fig.update_layout(
        height=1000,
        width=1500,
        title_text=os.path.basename(tsvgz_flpth)[:-3] + ' ' + str(df.shape),
        title_x=0.5
        )

    logging.info('Writing heatmap at %s' % html_flpth)
    t0 = time.time()
    fig.write_html(html_flpth)
    t_write_html = time.time() - t0

    dprint('t_write_html', run=locals())
    
    """ # orca server process reconnection error
    logging.info('Writing heatmap at %s' % png_flpth)
    t0 = time.time()
    fig.write_image(png_flpth)
    t_write_image = time.time() - t0

    dprint('t_write_image', run=locals())
    """
    
    logging.info('Generating static heatmap image')

    t0 = time.time()
    ax = sns.heatmap(df)
    t_sns_heatmap = time.time() - t0

    dprint('t_sns_heatmap', run=locals())

    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
    fig = ax.get_figure()
    fig.set_size_inches(20, 10)
    
    logging.info('Writing heatmap at %s' % png_flpth)

    t0 = time.time()
    fig.savefig(png_flpth)
    t_savefig = time.time() - t0

    dprint('t_savefig', run=locals())

    plt.clf()
        

class HTMLReportWriter:

    def __init__(self, cmd_l, tsvgz_flpth_l):
        '''
        out_files - [fixRank, filterByConf]
        params_prose - all str
        '''
        self.replacement_d = {}

        #
        self.tsvgz_flpth_l = tsvgz_flpth_l
        self.cmd_l = cmd_l

        #
        self.report_dir = os.path.join(_globals.shared_folder, str(uuid.uuid4()))
        os.mkdir(self.report_dir)


    def _compile_cmd(self):
        
        txt = ''

        for cmd in self.cmd_l:
            txt += '<pre><code>' + cmd + '</code></pre>\n'

        self.replacement_d['CMD_TAG'] = txt


    def _compile_figures(self):
        #
        self.fig_dir = os.path.join(self.report_dir, 'fig')
        os.mkdir(self.fig_dir)

        #
        png_flpth_l = [os.path.join(self.fig_dir, os.path.basename(tsvgz_flpth)[:-7] + '.png') 
            for tsvgz_flpth in self.tsvgz_flpth_l]
        html_flpth_l = [os.path.join(self.report_dir, os.path.basename(tsvgz_flpth)[:-7] + '.html') 
            for tsvgz_flpth in self.tsvgz_flpth_l]
        
        for tsvgz_flpth, png_flpth, html_flpth in zip(self.tsvgz_flpth_l, png_flpth_l, html_flpth_l):
            do_heatmap(tsvgz_flpth, png_flpth, html_flpth)


        def get_relative_fig_path(flpth):
            return '/'.join(flpth.split('/')[-2:])

        # build replacement string
        txt = '<div id="imgLink">\n'
        for png_flpth, html_flpth in zip(png_flpth_l, html_flpth_l):
            txt += '<p><a href="%s" target="_blank"><img alt="%s" src="%s" title="Open to interact"></a></p>\n' % (
                os.path.basename(html_flpth),
                os.path.basename(png_flpth),
                get_relative_fig_path(png_flpth))

        txt += '</div>\n'

        self.replacement_d['FIGURES_TAG'] = txt


    def write(self):
        self._compile_cmd()
        #self._compile_figures() # TODO stress test heatmaps

        
        REPORT_HTML_TEMPLATE_FLPTH = '/kb/module/ui/output/report.html'
        html_flpth = os.path.join(self.report_dir, 'report.html')

        with open(REPORT_HTML_TEMPLATE_FLPTH, 'r') as src_fp:
            with open(html_flpth, 'w') as dst_fp:
                for line in src_fp:
                    if line.strip() in self.replacement_d:
                        dst_fp.write(self.replacement_d[line.strip()])
                    else:
                        dst_fp.write(line)
        

        return self.report_dir, html_flpth







