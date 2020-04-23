import logging

from installed_clients.kb_GenericsReportClient import kb_GenericsReport

from .config import _globals




def get_html_dir_l(tsv_flpth_l):


    kbgr = kb_GenericsReport(_globals.callback_url)

    html_dir_l = []

    for tsv_flpth in tsv_flpth_l:
        logging.info('Building heatmap for %s with API-call to kb_GenericsReport' % tsv_flpth)
        html_dir = kbgr.build_heatmap_html({
            'tsv_file_path': tsv_flpth,
            'cluster_data': True,
            })['html_dir']

        html_dir_l.append(html_dir)

    return html_dir_l










