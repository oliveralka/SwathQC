import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
from swathqc.utils.misc import *


def id_over_rt_subplot(ax, df, color, min, title):
    if 'leftWidth' in df.columns:
        RTcolumn = 'RT'
    if 'LEFT_WIDTH' in df.columns:
        RTcolumn = 'EXP_RT'
    df[RTcolumn] = df.loc[:, RTcolumn].apply(lambda x: x/60)
    # bin RT range

    minRT = int(df.loc[:, RTcolumn].min())
    maxRT = int(df.loc[:, RTcolumn].max())
    if (maxRT - minRT) <= min:
        RTrange = np.arange(minRT, maxRT)
        ax.set_xlabel('RT [min]')
    else:
        RTrange = np.arange(minRT, maxRT, min)
        ax.set_xlabel('RT [in {} min intervals]'.format(str(min)))
    #RTrange = np.arange(int(df.loc[:, RTcolumn].min()), int(df.loc[:, RTcolumn].max()), min)

    xlabs = [min * i for i in range(len(RTrange))]

    # group by time interval and count
    # library patch
    ax.hist(df.loc[:, RTcolumn], bins=RTrange, edgecolor="k", color=color)
    plt.title(title)
    #plt.xlabel('RT[min]')
    plt.ylabel('# of IDs')
    plt.xticks(RTrange)

    # style x and y labels:
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(12)
    for label, i in zip(ax.xaxis.get_ticklabels(), xlabs):

         #label is a Text instance
        label.set_rotation(45)
        label.set_fontsize(12)
        #label.set_text(i)


def id_over_rt(df, min=5, color='coral', title='ID over RT', figsize=(10,5)):
    """

    :param df: pandas dataframe
    :param min: optional bin size, default 5 min
    :param title: optional title, default 'ID over RT'
    :param args:
    :param kwargs:
    :return: figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    # copy df before plotting!
    temp_df = df.copy()
    id_over_rt_subplot(ax, temp_df, color, min, title)
    fig.tight_layout()
    return fig


def report_id_over_rt(dfdict, color='coral', min=5, cols=2):
    # check of dict empty, this makes calling the fuction from run report easier.
    if not any(dfdict):
        return ''

    else:
        plt.rcParams.update({'font.size': 14})
        # get dataframe length and keys
        keys = dfdict.keys()
        N = len(dfdict)
        keyindex = zip(range(N), keys)
        # set up grid
        rows = int(math.ceil(N / cols))
        gs = gridspec.GridSpec(rows, cols)
        fig = plt.figure(figsize=(20, 7*rows))

        # plot each dataframe in dict

        for n, key in keyindex:
            if key.lower().endswith('osw'):
                subdict=dfdict[key]
                df = subdict['feature']
            else:
                df = dfdict[key]

            temp_df = df.copy()
            ax = fig.add_subplot(gs[n])
            id_over_rt_subplot(ax, temp_df, color, min, title=key)

        set_uni_ylimits(fig)
        set_uni_xlimits(fig)
        plt.suptitle('ID over RT')
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        return img_to_html(fig)


def describe(keys, min):
    """return html string"""
    html = "Plots the number of features identified in each retention time interval. Each bar represents a 5 minute RT interval. " \
           "Subplots are generated from any given swath (colored: red) and pyprophet or TRIC (both colored: blue) input files. " \
           "Please compare runs and make sure, the number of Feature IDs are stable over RT dimension and between runs. " \
           " <br>Files used to create plots: <ul><li>".format(min)
    html += "</li><li>".join(k for k in list(keys))
    html += "</ul>"
    return html


