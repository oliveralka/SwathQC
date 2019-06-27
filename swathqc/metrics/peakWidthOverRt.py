''' create plot and description for distribution of peak width over RT
Peak Width is calculated from leftWidth and rightWidth
(or maybe alternatively from fwhm?)'''
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from swathqc.utils.misc import *


def pw_over_rt_subplot(ax, df, title, onlyTarget=False, min=5):
    if onlyTarget:
        if check_decoy_column(df) is None:
            temp_df = df
        else:
            decoy = check_decoy_column(df)
            temp_df = df[df[decoy] == 0]
    else:
        temp_df = df

    if 'leftWidth' in df.columns:
        RTcolumn = 'RT'
        RWcolumn = 'rightWidth'
        LWcolumn = 'leftWidth'
    if 'LEFT_WIDTH' in df.columns:
        RTcolumn = 'EXP_RT'
        RWcolumn = 'RIGHT_WIDTH'
        LWcolumn = 'LEFT_WIDTH'
    df[RTcolumn] = df.loc[:, RTcolumn].apply(lambda x: x / 60)
    minRT = int(temp_df.loc[:, RTcolumn].min())
    maxRT = int(temp_df.loc[:, RTcolumn].max())
    if (maxRT - minRT) <= min:
        RTrange = np.arange(minRT, maxRT)
        ax.set_xlabel('RT [min]')
    else:
        RTrange = np.arange(minRT, maxRT, min)
        ax.set_xlabel('RT [in {} min intervals]'.format(str(min)))


    xlabs = [min * i for i in range(len(RTrange))]
    # subset dataframe
    df_to_plot = temp_df.loc[:, [RWcolumn, LWcolumn, RTcolumn]]
    # calculate width of peak
    df_to_plot['widthDiff'] = df_to_plot.loc[:, RWcolumn] - df_to_plot.loc[:, LWcolumn]
    df_to_plot['bin'] = pd.cut(df_to_plot.loc[:, RTcolumn], RTrange)
    df_to_plot.boxplot(column='widthDiff', by='bin',  rot=90, grid=False, ax=ax)

    ax.set_ylabel('PeakWidth')
    plt.title(title)
    plt.suptitle("")

    xtickNames = plt.setp(ax, xticklabels=xlabs)
    plt.setp(xtickNames, rotation=45, fontsize=12)



def pw_over_rt(df, title='PeakWidth over RT', onlyTarget=True, figsize=(10,5), min=5):

    # set ax
    fig, ax = plt.subplots(figsize=figsize)
    temp_df = df.copy()
    # only plot from features identified as targets

    # check and set column names
    pw_over_rt_subplot(ax, temp_df, title, onlyTarget, min)
    fig.tight_layout()
    return fig


def report_pw_over_rt(dfdict, suptitle='PeakWidth over RT', min=5, cols=1):
    """
    boxplots of Peak Width for given dict of dataframesw, arrange plots deppending on size of dataframe
    determine osw or tsv source from dataframe keys
    :param dfdict: dict of dataframes, either from tsv or osw FEATURE table
    :param min: time interval on RT, default 5 min
    :param cols: number of subplots per column, default: 1
    :return: tuple: (fig as html, description as html, sectionid as string)
    """
    if not any(dfdict):
        return ''
    else:

        keys = dfdict.keys()
        N = len(dfdict)
        keyindex = zip(range(N), keys)

        rows = int(math.ceil(N / cols))
        gs = gridspec.GridSpec(rows, cols)
        fig = plt.figure(figsize=(10, 4 * rows))


        for n, key in keyindex:
            if key.lower().endswith('osw'):
                subdict = dfdict[key]
                df = subdict['feature']
            if key.lower().endswith('tsv'):
                df = dfdict[key]
            # set ax
            temp_df = df.copy()
            ax = fig.add_subplot(gs[n])
            pw_over_rt_subplot(ax, temp_df, key)
        set_uni_ylimits(fig)
        plt.suptitle(suptitle)
        fig.tight_layout()
        fig.subplots_adjust(top=(1-0.14/rows))
        return img_to_html(fig)


def describePWoverRT(keys, minutes):

    html = "This plot shows the variation of peak width over the RT dimension. <br> The data for these plots was taken from " \
           "OpenSwathWorkflow output: <br><ul><li>"
    html+="</li><li>".join(k for k in list(keys))
    html = html + "</ul>Each boxplot corresponds to a {} minute time interval and shows mean peak wifth and variation." \
           "Peak Width is calculated from the leftWidth and rightWidth columns. Alternativley, full widh at half" \
           "maximum (fwhm) could also be used to represent this metric".format(minutes)
    return html