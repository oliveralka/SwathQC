'''
from: swath tsv only
'''
import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
from swathqc.utils.misc import *


def describe(data):
    data = data.to_html
    html = "For each sequence, OpenSwathWorkflow does an in situ digest (assuming Trypsin) and places the number of missed " \
           "cleavages for each peptide in a column 'MC'. This columns can currently only be found in the tsv output from " \
           "OpenSwathWorkflow, so this plot will only be available if swath tsv file(s) were provided.<br>" \
           "Processing data: all peptides tagged as decoys are dropped before plotting and transitions groups are grouped" \
           "so that each peptide is counted only once.<br> QC target: sample preparation"
    html += data
    return html


def missed_cleavages_data(df, onlyTarget, transitionCol='transition_group_id', mcCol='MC'):
    df = subset_only_targets(df, onlyTarget)
    # subset data
    temp_df = df[[transitionCol, mcCol]].drop_duplicates(transitionCol).groupby(
        mcCol).count().reset_index()
    temp_df = temp_df.rename(columns={'transition_group_id': 'peptides'})
    return temp_df


def missed_cleavages_subplot(ax, df, color, title, onlyTarget, transitionCol='transition_group_id', mcCol='MC'):
    props = dict(boxstyle='round', facecolor='w', alpha=0.5)
    # place a text box in upper left in axes coords
    # onyl plot if missed cleavages not null
    temp_df = missed_cleavages_data(df, onlyTarget, transitionCol, mcCol)
    # peptides without missed cleavages
    nomisses = temp_df[temp_df[mcCol] == 0].iloc[0, 1]
    textstr = "Peptides without \nmissed cleavage: \n{}".format(nomisses)
    temp_df[temp_df[mcCol] != 0]['peptides'].plot(kind='bar', ax=ax, color=color, legend=None)
    ax.set_alpha(0.8)
    ax.set_title(title)  # , fontsize=14)
    ax.set_xlabel('# missed cleavages')
    ax.set_ylabel('count (transitions groups)')
    ax.xaxis.set_tick_params(rotation=0)
    ax.text(0.6, 0.9, textstr, transform=ax.transAxes, verticalalignment='top', bbox=props)
    b, t = ax.get_ylim()
    for i in ax.patches:
        ax.text(i.get_x(), i.get_height() + t / 100, str(i.get_height()))


def missed_cleavages(df, color='coral', title='PeakWidth over RT', onlyTarget=True, figsize=(10.5)):
    fig, ax = plt.subplots(figsize=figsize)
    temp_df = df.copy()
    # only plot from features identified as targets

    # check and set column names
    missed_cleavages_subplot(ax, temp_df, color, title, onlyTarget)
    fig.tight_layout()
    return fig


def report_missed_cleavages(dfdict, color='coral', suptitle='Missed Cleavages in peptides', onlyTarget=True, cols=3):

    #if not any(dfdict):
     #   return ('<p class="missing">can only be plotted from swath TSV </p>', describe(), 'masserror')

    keys = dfdict.keys()
    N = len(dfdict)
    if N < cols:
        cols = N
    keyindex = zip(range(N), keys)
    rows = int(math.ceil(N / cols))

    gs = gridspec.GridSpec(rows, cols)
    fig = plt.figure(figsize=(12, 7 * rows))

    for n, key in keyindex:
        ax = fig.add_subplot(gs[n])
        df = dfdict[key]
        missed_cleavages_subplot(ax, df, color, title=key, onlyTarget=onlyTarget)

    set_uni_ylimits(fig)
    fig.suptitle(suptitle)
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)

    return img_to_html(fig)