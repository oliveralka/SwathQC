'''
plot a venn diagram of coverage of he assay library, based on peptides
'''
from __future__ import absolute_import
import math
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib_venn import venn3, venn2
from swathqc.utils.misc import *


# TODO: library coverage on protein, peptide and ?? level



def checkFilyType(key):
    if key.lower().endswith('osw'):
        return 'osw'
    if key.lower().endswith('tsv'):
        return 'tsv'
    if key.lower().endswith('featureXML'):
        return 'featureXML'

def description():
    description= "Venn diagram of library coverage. This can only be plotted if the library assay list was provided in tsv format. " \
       "Libraries in pqp format are currently not supported. <br>" \
       "Library coverage will be plotted from tsv files, if they were provided, otherwise it can also be plotted from osw files. <br>" \
       ""
    return description

def library_coverage_subplot(ax, library, file, title, libraryCol, dataCol, decoyCol):
    """
    add venn diagram to represent library coverage and decoys on given axis,
    call this function from plot
    :param ax: axes to plot on
    :param library: dataframe of library from tsv
    :param file: dataframe of swath output, osw or tsv
    :param source: either tsv or osw, tet the right column names
    :param title: string, title of subplot
    :return:
    """

    # process dataframes into sets for venn diagram
    # for DIA: drop (UniMod:4) from peptide string
    dia_peptides = set(file[file[decoyCol] == 0].loc[:, dataCol])
    lib_peptides = set(library[library[decoyCol] == 0].loc[:, libraryCol])


    # plot
    v = venn2([lib_peptides, dia_peptides], set_labels=('Library', 'DIA', ' '), ax=ax)
    # library patch
    v.get_patch_by_id('10').set_alpha(0.8)
    v.get_patch_by_id('10').set_color('coral')

    # DIA patch
    v.get_patch_by_id('11').set_alpha(0.8)
    v.get_patch_by_id('11').set_color('#3A78A4')
    # decoy patch
    v.get_patch_by_id('01').set_alpha(0.8)
    v.get_patch_by_id('01').set_color('#81e448')

    # adjust subset label positions
    v.get_label_by_id("A").set_x(-0.4)

    b = v.get_label_by_id("B")
    x, y = b.get_position()
    b.set_position((x + 0.1, y + 0.2))

    # d = v.get_label_by_id("C")
    # x, y = d.get_position()
    # d.set_position((x + 0.1, y))
    ax.set_title(title)



def library_cov_subplot(ax, library, file, title, libraryCol, dataCol, decoyCol):
    """
    add venn diagram to represent library coverage and decoys on given axis,
    call this function from plot
    :param ax: axes to plot on
    :param library: dataframe of library from tsv
    :param file: dataframe of swath output, osw or tsv
    :param source: either tsv or osw, tet the right column names
    :param title: string, title of subplot
    :return:
    """

    # process dataframes into sets for venn diagram
    # for DIA: drop (UniMod:4) from peptide string
    dia_peptides = set(file[file[decoyCol] == 0].loc[:, dataCol])
    lib_peptides = set(library[library[decoyCol] == 0].loc[:, libraryCol])
    decoys = set([0])

    # plot
    v = venn3([lib_peptides, dia_peptides, decoys], set_labels=('Library', 'DIA', ' '), ax=ax)
    # library patch
    v.get_patch_by_id('100').set_alpha(0.8)
    v.get_patch_by_id('100').set_color('coral')

    # DIA patch
    v.get_patch_by_id('110').set_alpha(0.8)
    v.get_patch_by_id('110').set_color('#3A78A4')
    # decoy patch
    v.get_patch_by_id('001').set_alpha(0.8)
    v.get_patch_by_id('001').set_color('#81e448')

    # adjust subset label positions
    v.get_label_by_id("A").set_x(-0.4)

    b = v.get_label_by_id("B")
    x, y = b.get_position()
    b.set_position((x + 0.1, y + 0.2))

    # d = v.get_label_by_id("C")
    # x, y = d.get_position()
    # d.set_position((x + 0.1, y))
    ax.set_title(title)

def library_coverage(lib, df, figsize=(5,5), title='Library Coverage (peptides)', libraryCol='UNMODIFIED_SEQUENCE',
                     dataCol='UNMODIFIED_SEQUENCE', decoyCol='DECOY'):
    """
    :param: library as dataframe
    :param dataframe with the columns
    :return: figure
    """
    temp_lib = lib.copy()
    temp_df = df.copy()
    # Todo, check if library has decoys!!
    if 'PeptideSequence' in temp_lib.columns:
        temp_lib = temp_lib.rename(columns={'PeptideSequence':'UNMODIFIED_SEQUENCE', 'decoy': 'DECOY'})
    if 'FullPeptideName' in temp_df.columns:
        temp_df = temp_df.rename(columns={'FullPeptideName': 'UNMODIFIED_SEQUENCE', 'decoy': 'DECOY'})


    fig, ax = plt.subplots(figsize=figsize)
    library_coverage_subplot(ax, temp_lib, temp_df, title, libraryCol, dataCol, decoyCol)

    return fig


def report_libray_coverage(lib, dfdict, suptitle='Library Coverage (peptides)', libraryCol='UNMODIFIED_SEQUENCE',
                     dataCol='UNMODIFIED_SEQUENCE', decoyCol='DECOY', cols=3):
    """
    make subplot for each given dataframe in dictionary
    :param lib: library as dataframe
    :param dfdict: dictionary of dataframes to be plotted
    :return: figure
    """
    if not any(dfdict):
        return ''

    else:
        if isinstance(lib, dict):
            lib = lib['peptide']
        if 'PeptideSequence' in lib.columns:
            lib = lib.rename(columns={'PeptideSequence': 'UNMODIFIED_SEQUENCE', 'decoy':'DECOY'})
        N = len(dfdict)
        rows = int(math.ceil(N / cols))

        gs = gridspec.GridSpec(rows, cols)
        fig = plt.figure(figsize=(20, 7 * rows))

        keys = dfdict.keys()
        keyindex = zip(range(N), keys)
        for n, key in keyindex:
            filetype = checkFilyType(key)
            df = dfdict[key]
            if filetype == 'tsv':
                df = df.rename(columns={'decoy': 'DECOY', 'FullPeptideName': 'UNMODIFIED_SEQUENCE'})
            # set ax
            ax = fig.add_subplot(gs[n])
            # plot on ax
            library_coverage_subplot(ax, lib, df, key, libraryCol, dataCol, decoyCol)
        fig.subplots_adjust(hspace=.2, wspace=.001)
        plt.suptitle(suptitle)

        return img_to_html(fig)






