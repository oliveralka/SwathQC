import matplotlib.pyplot as plt
from matplotlib import gridspec
from swathqc.utils.misc import *


def dscore_osw_subplot(df, title, gs_inner, fig):
    contexts = set(df.CONTEXT)
    for i, context in zip(range(len(contexts)), contexts):
        sub_df = df[df.CONTEXT == context]
        ax = plt.Subplot(fig, gs_inner[i])
        sub_df[sub_df.decoy == 0].d_score.plot.density(label='Targets', title=title, ax=ax)
        sub_df[sub_df.decoy == 1].d_score.plot.density(label='Decoys',  ax=ax)
        #sub_df[sub_df.decoy == 0].d_score.plot.hist(color='blue', bins=40, alpha=0.5, label="Targets", ax=ax)
        #sub_df[sub_df.decoy == 1].d_score.plot.hist(color='green', bins=40, alpha=0.5,   label='Decoys', ax=ax)
        ax.set_xlabel('d-score')
        fig.add_subplot(ax)


def dscore_tsv_subplot(fig, df, key, gs_inner):
    ax1 = plt.Subplot(fig, gs_inner[0])
    #ax2 = plt.Subplot(fig, gs_inner[1])
    df[df.decoy == 0].d_score.plot.density(label='Targets', ax=ax1, title=key)
    df[df.decoy == 1].d_score.plot.density(label='Decoys', ax=ax1)

    #df[df.decoy == 0].d_score.plot.hist(color='blue', bins=40, alpha=0.5, label="Targets", ax=ax2)
    #df[df.decoy == 1].d_score.plot.hist(color='green', bins=40, alpha=0.5,   label='Decoys', ax=ax2)


def pyprophet_dscore(df, title='', figsize=(10, 5)):

    #Todo: for package, read in data from file location insetad of dataframe?

    fig, ax = plt.subplots(figsize=figsize)
    temp_df = df.copy()

    fig.tight_layout()
    return fig



def report_dscore(dfdict, suptitle='d-score density plots'):
    keys = dfdict.keys()
    N = len(dfdict)
    keyindex = zip(range(N), keys)
    rows = N

    gs_outer = gridspec.GridSpec(rows, 1)
    fig = plt.figure(figsize=(12, 7 * rows))
    for n, key in keyindex:
        df = dfdict[key]
        if 'osw' in key:
            # unpack
            subdf_peptide= df['peptideScores']
            #subdf_protein = df['proteinScores']
            #TODO: fix plotting from osw
            subcols = 3
            gs_inner = gridspec.GridSpecFromSubplotSpec(1, subcols, subplot_spec=gs_outer[n], wspace=0.2, hspace=0.2)
            dscore_osw_subplot(subdf_peptide, key, gs_inner, fig)

        if 'tsv' in key:
            ax = fig.add_subplot(gs_outer[n])
            df[df.decoy == 0].d_score.plot.density(label='Targets', ax=ax, title=key, linewidth=3, legend=True)
            df[df.decoy == 1].d_score.plot.density(label='Decoys', ax=ax, linewidth=3, legend=True)
            ax.set_xlabel('d-score')
            #df[df.decoy == 0].d_score.plot.hist(color='blue', bins=40, alpha=0.5, label="Targets", ax=ax)
            #df[df.decoy == 1].d_score.plot.hist(color='green', bins=40, alpha=0.5,   label='Decoys', ax=ax)

    plt.suptitle(suptitle)
    set_uni_ylimits(fig)
    fig.tight_layout()
    fig.subplots_adjust(top=(1-0.14/rows))


    return img_to_html(fig)



#Todo: implement plot mscores/Qvalues

def mscore_plot_osw(df, key, gs_inner, fig):
    contexts = set(df.CONTEXT)
    for i, context in zip(range(len(contexts)), contexts):
        sub_df = df[df.CONTEXT == context]
        ax = plt.Subplot(fig, gs_inner[i])
        sub_df[sub_df.decoy == 0].d_score.plot.density(label='Targets', title=(context+key[-8:]), ax=ax)
        sub_df[sub_df.decoy == 1].d_score.plot.density(label='Decoys',  ax=ax)
        sub_df[sub_df.decoy == 0].d_score.plot.hist(color='blue', bins=40, alpha=0.5, label="Targets", ax=ax)
        sub_df[sub_df.decoy == 1].d_score.plot.hist(color='green', bins=40, alpha=0.5,   label='Decoys', ax=ax)
        ax.set_xlabel('d-score')
        fig.add_subplot(ax)


def mscore_plot_tsv(fig, df, key, gs_inner):

    ax1 = plt.Subplot(fig, gs_inner[0])
    ax2 = plt.Subplot(fig, gs_inner[1])
    df[df.decoy == 0].d_score.plot.density(label='Targets', ax=ax1, title=key)
    df[df.decoy == 1].d_score.plot.density(label='Decoys', ax=ax1)

    df[df.decoy == 0].d_score.plot.hist(color='blue', bins=40, alpha=0.5, label="Targets", ax=ax2)
    df[df.decoy == 1].d_score.plot.hist(color='green', bins=40, alpha=0.5,   label='Decoys', ax=ax2)


def qvalue_plot(dfdict):
    keys = dfdict.keys()
    N = len(dfdict)
    keyindex = zip(range(N), keys)
    rows = N

    gs_outer = gridspec.GridSpec(rows, 1)
    fig = plt.figure(figsize=(20, 7 * rows))
    for n, key in keyindex:
        df = dfdict[key]
        contexts = set(df.CONTEXT)
        gs_inner = gridspec.GridSpecFromSubplotSpec(1, subplot_spec=gs_outer[n], wspace=0.2, hspace=0.2)
        dscore_osw_subplot(df, key, gs_inner, fig)

    plt.suptitle('Q Value plot')


    html = "d-scores <br>" \
               "density plot or histogram?".format()
    return (img_to_html(fig), html, 'dscore_density')
