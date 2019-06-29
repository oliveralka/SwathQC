import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
from swathqc.utils.misc import *


def num_peaks(entrystr):
    temp_list = [i for i in entrystr.split(';') if float(i) != 0]
    return len(temp_list)


def num_transition_data(df, fileType='osw'):
    if fileType == 'tsv':
        df['peak_count'] = df['aggr_Peak_Area'].apply(lambda x: num_peaks(x))
        to_plot = df[df['decoy'] == 0].groupby('peak_count').count().reset_index().rename(columns={'transition_group_id':'number of available peaks'})
    if fileType == 'osw':
        # from FEATURE_TRANSITION table
        # either APEX_INTENSITY or AREA_INTENSITY
        to_plot = df[df.APEX_INTENSITY != 0.0].groupby('FEATURE_ID').count().groupby('TRANSITION_ID').count().reset_index().rename(
            columns={'TRANSITION_ID': 'peak_count', 'AREA_INTENSITY': 'number of available peaks'})

    return to_plot


def num_transitions_subplot(ax, df, title, color):
    """number of transitions for each identified peptide
    datasource: pptsv, swathtsv, stwathosw:feature_transition table
    :param pptsv dataframe
    :return figure as html
    """
    # from tsv:
    df.plot.bar(x='peak_count', y='number of available peaks', color=color, edgecolor='k', rot=0, ax=ax, legend=None)

    ax.set_title(title)

    ax.set_ylabel("number of features")
    ax.set_xlabel("number of transitions")
    # tooltips
    #for i, bar in enumerate(ax.get_children()):
    #    tooltip = mpld3.plugins.LineLabelTooltip(bar, label=str(i))
    #    mpld3.plugins.connect(plt.gcf(), tooltip)
    b, t = ax.get_ylim()
    for i in ax.patches:
        ax.text(i.get_x(), i.get_height() + t/100, str(i.get_height()))


def num_transitions(df, title='Num of Transitions/Feature', color='coral', figsize=(10,10)):
    fig, ax = plt.subplots(figsize=figsize)
    # copy df before plotting!

    temp_df = df.copy()
    if 'aggr_Peak_Area' in temp_df.columns:
        temp_df['peak_count'] = temp_df['aggr_Peak_Area'].apply(lambda x: num_peaks(x))
        to_plot = temp_df[temp_df['decoy'] == 0].groupby('peak_count').count().reset_index().rename(
            columns={'transition_group_id': 'number of available peaks'})
    if "FEATURE_ID" in temp_df.columns:
        # from FEATURE_TRANSITION table
        # either APEX_INTENSITY or AREA_INTENSITY
        to_plot = temp_df[temp_df.APEX_INTENSITY != 0.0].groupby('FEATURE_ID').count().groupby(
            'TRANSITION_ID').count().reset_index().rename(
            columns={'TRANSITION_ID': 'peak_count', 'AREA_INTENSITY': 'number of available peaks'})

    num_transitions_subplot(ax, to_plot, title, color)
    fig.tight_layout()
    return fig


def report_num_transitions(dfdict, color='coral', cols=3):

    if not any(dfdict):
        return ''
    else:
        keys = dfdict.keys()
        N = len(dfdict)
        if N < cols:
            cols = N
        keyindex = zip(range(N), keys)
        rows = int(math.ceil(N / cols))
        gs = gridspec.GridSpec(rows, cols)
        fig = plt.figure(figsize=(10, 7 * rows))

        # collect yaxis limits
        ylimits = []
        # add subplot for each passed dataframe
        for n, key in keyindex:
            fileType= check_file_type(key)
            if fileType =='tsv':
                df = dfdict[key]
                df_temp = df.copy()
                to_plot = num_transition_data(df_temp, 'tsv')
            if fileType == 'osw':
                subdict = dfdict[key]
                df = subdict['featureTransition']
                df_temp = df.copy()
                to_plot = num_transition_data(df_temp, 'osw')

            ax = fig.add_subplot(gs[n])
            # plot on ax
            num_transitions_subplot(ax, to_plot, key, color)

        set_uni_ylimits(fig)

        plt.suptitle('Number of Transitions per Feature', size=14)
        fig.tight_layout()
        fig.subplots_adjust(top=(1-0.12/rows))
        return img_to_html(fig)


def describe():
    html = "The library contains (or at least  it should) 6 transitions for each feature. " \
           "For each feature, the number of recovered transitions are counted from the column " \
           " 'aggr_peak_area' in tsv files and and 'AREA_INTENSITY' from the FEATURE_TRANSITION_TABLE in osw files." \
           "Only non zero area intensities will be counted. " \
           "Files from pyprophet or TRIC contains all features and transitions from merged swath files, so the total " \
           "number will be much higher. But the number of features with only 1 - 3 transitions should go down considerably, " \
           "as these low quality features get filtered out."


    return html