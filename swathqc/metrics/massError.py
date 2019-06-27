import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
from swathqc.utils.misc import *


def mean_error(entrystr):
    """for given string, that contains Apek m/z and masserror separated by semicolon
    split string into list, take every second element (the masserror)
    :return list of mean masserror """
    temp_list = entrystr.split(';')
    tmp_list = list(map(float, temp_list))
    return np.mean(list(tmp_list[i] for i in range(1, len(tmp_list), 2)))


def flatten_mass_errors(masserrorlist):
    """ for list of strings, that contain mass_errors, flatten into one list.
    :return: list of masserrors"""
    return_list = []
    for errors in masserrorlist:
        temp_list = errors.split(';')
        temp_list = list(map(float, temp_list))
        tmp_list = [temp_list[i] for i in range(1, len(temp_list),2)]
        return_list.append(tmp_list)
    return_list = [y for x in return_list for y in x]
    return return_list


def mass_error_subplot(ax, df, color, title, onlyTarget):

    df = subset_only_targets(df, onlyTarget)

    # masserror column into list of string
    masserrors = list(df.loc[:, 'masserror_ppm'].dropna())
    # lsit of mean mass errors
    #meanerrors = [mean_error(i) for i in masserrors]

    errors = flatten_mass_errors(masserrors)


    bins = np.arange(min(errors), max(errors), 1)
    mu = np.mean(errors)
    median = np.median(errors)
    sigma = np.std(errors)
    textstr = '\n'.join((
        r'$\mu=%.2f$' % (mu,),
        r'$\mathrm{median}=%.2f$' % (median,),
        r'$\sigma=%.2f$' % (sigma,)))

    ax.hist(errors, bins=bins, alpha=0.8, color=color)
    plt.title(title)
    plt.xlabel('masserror [ppm]')
    plt.ylabel('count')
    # vertical line for mean
    # plt.axvline(mean,color='b', linewidth=1)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='w', alpha=0.5)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)


def massdev_score_subplot(ax, df, color, title, onlyTarget):
    # TODO: don't really need that one, remember to delete
    df = subset_only_targets(df, onlyTarget)
    """Plot va_massdev_score from osw FEATURE_MS2"""
    masserrors = list(df.loc[:, 'VAR_MASSDEV_SCORE'].dropna())
    # lsit of mean mass errors

    bins = np.arange(min(masserrors), max(masserrors), 1)
    ax.hist(masserrors, bins=bins, alpha=0.8, color=color)
    plt.title(title)
    plt.xlabel('mass_dev_score')
    plt.ylabel('count')


def mass_error(df, color='coral', title='Masserror ppm', onlyTarget=False, figsize=(5,5)):
    fig, ax = plt.subplots(figsize=figsize)
    temp_df = df.copy()
    # only plot from features identified as targets

    # check and set column names
    mass_error_subplot(ax, temp_df, color, title, onlyTarget)
    fig.tight_layout()
    return fig


def describe():

    html = "From TSV: histogram of mean mass errors of each feature. The mean is calculated from the mass errors of all ions in a feature." \
           "<br>Fom OSW: plot massdev_score, a subscore calculated by OpenSwathWorkflow. A massdev_score close to 0 is desired"
    return html


def report_mass_error(dfdict, color='coral', suptitle='Masserror ppm', cols=3):
    #plt.rcParams.update({'font.size': 14})

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
        if key.lower().endswith('osw'):
            ax.text('can only be plotted from swath TSV files', ha="center", va="center")

        mass_error_subplot(ax, df, color, title=key, onlyTarget=False)


    #set_uni_ylimits(fig)
    #set_uni_xlimits(fig)
    plt.suptitle(suptitle)
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)

    return img_to_html(fig)
# --------------------------------------
# two dimensional histogram, not needed

def two_dim_hist_subplot(ax, df, key):

    df_to_plot = df[['masserror_ppm', 'm/z']].dropna()
    df_to_plot['meanerror'] = df_to_plot['masserror_ppm'].apply(lambda x: mean_error(x))
    xdata = df_to_plot['m/z']
    ydata = df_to_plot['meanerror']
    bins = np.arange(min(xdata), max(xdata), 1)

    ax.hist2d(xdata, ydata, bins=[100, 100], alpha=0.5)
    plt.ylabel('mean masserror of transition')
    plt.ylabel('m/z')
    plt.title(key)


def plt_two_dim_hist(dfdict, cols=3):
    keys = dfdict.keys()
    N = len(dfdict)
    if N < cols:
        cols = N
    keyindex = zip(range(N), keys)
    rows = int(math.ceil(N / cols))

    gs = gridspec.GridSpec(rows, cols)
    fig = plt.figure(figsize=(20, 10 * rows))

    for n, key in keyindex:
        ax = fig.add_subplot(gs[n])
        df = dfdict[key]
        two_dim_hist_subplot(ax, df, key)
    plt.suptitle('2D Histogram masserrors')
    #plt.colorbar()

    return img_to_html(fig)
