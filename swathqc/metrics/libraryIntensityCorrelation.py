import math

import matplotlib.pyplot as plt
from swathqc.utils.misc import *
from matplotlib import gridspec


def lib_intensity_corr_subplot(ax, df, colname, title):
    """ boxplot of column 'var_library_corr' or 'VAR_LIBRARY_CORR (in all except pp exprot tsv)"""
    lib_corr = df[colname]
    # get mean, median and std error to place in textbox
    mu = np.mean(lib_corr)
    median = np.median(lib_corr)
    sigma = np.std(lib_corr)
    textstr = '\n'.join((
        r'$\mu=%.2f$' % (mu,),
        r'$\mathrm{median}=%.2f$' % (median,),
        r'$\sigma=%.2f$' % (sigma,)))
    ax.boxplot(lib_corr)
    #plt.title('Intensity Correlation with Library')
    ax.set_xlabel(title)
    ax.set_ylabel(colname)
    # vertical line for mean
    # plt.axvline(mean,color='b', linewidth=1)

    # these are matplotlib.patch.Patch properties
    props = dict(facecolor='w', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)


def lib_intensity_corr(df, colname='VAR_INTENSITY_CORR', title='Library intensity correlation', figsize=(10,10)):
    fig, ax = plt.subplots(figsize=figsize)
    # copy df before plotting!

    temp_df = df.copy()
    if colname not in df.columns:
        if colname.lower() in df.columns:
            colname = colname.lower()
    lib_intensity_corr_subplot(ax, temp_df, colname, title=title)
    fig.tight_layout()
    return fig


def report_lib_intensity_corr(dfdict, suptitle='Library Intensity Correlation', cols=3):
    # osw: FEATURE_MS2
    if not any(dfdict):
        return ''
    else:
        keys = dfdict.keys()

        N = len(dfdict)
        if N <= cols:
            cols = N

        keyindex = zip(range(N), keys)
        rows = int(math.ceil(N / cols))
        gs = gridspec.GridSpec(rows, cols)
        fig = plt.figure(figsize=(10, 5 * rows))

        for n, key in keyindex:
            if key.lower().endswith('tsv'):
                colname = 'var_library_corr'
                df = dfdict[key]
                if colname not in df.columns:
                    return 'plot not available from pyprophet tsv export'
            if key.lower().endswith('osw'):
                subdict = dfdict[key]
                df = subdict['feature']
                colname = 'VAR_LIBRARY_CORR'

            # set ax
            ax = fig.add_subplot(gs[n])
            lib_intensity_corr_subplot(ax, df, colname, key)

        fig.subplots_adjust(hspace=.2, wspace=.001)
        fig.tight_layout()
        fig.suptitle(suptitle)
        fig.subplots_adjust(top=0.9)
        return img_to_html(fig)