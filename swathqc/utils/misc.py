import base64
from io import BytesIO
import numpy as np


def check_file_type(key):
    if key.lower().endswith('osw'):
        return 'osw'
    if key.lower().endswith('tsv'):
        return 'tsv'
    if key.lower().endswith('featureXML'):
        return 'featureXML'


def subset_only_targets(df, onlyTargets):
    if onlyTargets==True:
        if check_decoy_column(df) is None:
            ret_df = df
        else:
            decoy = check_decoy_column(df)
            ret_df = df[df[decoy] == 0]
    else:
        ret_df = df
    return ret_df


def img_to_html(figure):
    """
      @brief
      @param
      @return
      """
    buffer = BytesIO()
    figure.savefig(buffer, format='png')
    buffer.seek(0)

    data_uri = base64.b64encode(buffer.read()).decode('ascii')
    html = '<img class="myImages" src="data:image/png;base64,{0}" width=100%>'.format(data_uri)
    return html


def check_decoy_column(df):
    if 'decoy' in df.columns:
        decoy = 'decoy'
    else:
        if 'DECOY' in df.columns:
            decoy = 'DECOY'
        else:
            decoy = None
    return decoy


def set_uni_ylimits(fig):
    """
    set same y axis for all subplots in figure
    get all ylimits and set to tallest ylimit
    :param fig:
    :return:
    """
    axes = fig.axes
    ylimits = []
    for n, subplot in np.ndenumerate(axes):
        ylimits.append(subplot.get_ylim())
    top = max(ylimits, key=lambda item: item[1])[1]
    for n, subplot in np.ndenumerate(axes):
        subplot.set_ylim(0, top)


def set_uni_xlimits(fig):
    """
    set same x axis for all subplots in figure
    :param fig:
    :return:
    """
    axes = fig.axes
    xlimits = []
    for n, subplot in np.ndenumerate(axes):
        xlimits.append(subplot.get_xlim())
    top = max(xlimits, key=lambda item: item[1])[1]
    for n, subplot in np.ndenumerate(axes):
        subplot.set_xlim(0, top)