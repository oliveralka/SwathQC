import math
import matplotlib.pyplot as plt
from swathqc.utils.misc import *
from matplotlib import gridspec
from scipy import stats
import numpy as np


plt.rcParams.update({'font.size': 14})
# Todo: adjust this metric
# Todo: only plotting from pp tsv export?

def describe():
    html = "Scatter plot of run RT against RT difference between library and run." \
           "<br>can be plotted from pyprophet or tric output"
    return html


def add_linear_regress(ax, df, xcol, ycol):
    x = df[xcol]
    y = df[ycol]
    # zero line
    ax.plot([min(x), max(x)], [0, 0], 'g-')

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line_x = np.arange(x.min(), x.max())
    line_y = slope * line_x + intercept
    ax.plot(line_x, line_y, 'r-', label='$%.2fx + %.2f$, $R^2=%.2f$' % (slope, intercept, r_value ** 2))
    ax.legend()


def rt_correlation_subplot(ax, df, title,  x_rt='RT', y_rt ='delta_rt', decoy='decoy', group='transition_group_id', mode='mean'):
    # mode can be either mean, byrun, or all
    # if iRT, use iRT, else use RT
    if 'iRT' in df.columns:
        x_rt = 'iRT'
        y_rt = 'delta_iRT'
        xlab = 'iRT [min]'
        ylab = 'delta iRT [min]'
        if 'LIBRARY_RT' in df.columns:
            df['delta_iRT'] = df['iRT'] - df['LIBRARY_RT']
            group = 'precursor'
    else:
        xlab = 'RT [s]'
        ylab = 'delta RT [s]'
    # convert rt to minutes
    #df[y_rt] = df['delta_rt'].apply(lambda x: x / 60)
    #df['RT'] = df['RT'].apply(lambda x: x / 60)
    if mode == 'mean':
        df[df[decoy] == 0].groupby(group).mean().plot.scatter(x=x_rt, y=y_rt, ax=ax)
    if mode == 'byrun':
        pass
        # use different color for each run
    if mode == 'all':
        df[df[decoy] == 0].plot.scatter(x=x_rt, y=y_rt, ax=ax)



    # plot mean for all features

    # plot delta for every feature
    #df[df[decoy] == 0].plot.scatter(x=x_rt, y=y_rt, ax=ax)

    plt.ylabel(ylab)
    plt.xlabel(xlab)
    add_linear_regress(ax, df, x_rt, y_rt)
    plt.title(title)


def rt_correlation(df, title='RT correlation', xcol='RT', ycol ='delta_rt', decoycol='decoy', groupcol='transition_group_id', figsize=(10, 7), mode='mean'):
    fig, ax = plt.subplots(figsize=figsize)
    temp_df = df.copy()
    # only plot from features identified as targets
    rt_correlation_subplot(ax, temp_df, title, xcol, ycol, decoycol, groupcol, mode)
    # check and set column names
    fig.tight_layout()
    return fig


def report_rt_correlation(dfdict, suptitle='RT correlation with library'):
    if not any(dfdict):
        return ''
    else:

        cols = 1
        keys = dfdict.keys()
        N = len(dfdict)
        if N < cols:
            cols = N
        keyindex = zip(range(N), keys)
        rows = int(math.ceil(N / cols))

        gs = gridspec.GridSpec(rows, cols)
        fig = plt.figure(figsize=(12, 7 * rows))

        for n, key in keyindex:
            if key.lower().endswith('osw'):
                subdict = dfdict[key]
                df = subdict['rt']
            if key.lower().endswith('tsv'):
                df = dfdict[key]

            ax = fig.add_subplot(gs[n])
            rt_correlation_subplot(ax, df, key)
            plt.subplots_adjust(hspace=.2, wspace=.001)
        fig.tight_layout()
        fig.subplots_adjust(top=(1-0.12/rows))
        fig.suptitle(suptitle)

        return img_to_html(fig)

