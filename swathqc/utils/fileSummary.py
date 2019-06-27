""" defualt column names are for libraries
- get number of peptide, proteins, transitions, pretien groups and precursors
-
"""

import pandas as pd
# TODO: transitions are now only from osw, maybe include from tsv again
# constants: columns
PYPROPHET_TSV_COLUMNS={'peptideColumn' : 'FullPeptideName',
           'proteinColumn':'ProteinName',
           'transitionColumn':'aggr_Peak_Area',
           'precursorColumn': None}

LIBRARY_COLUMNS={'peptideColumn' : 'PeptideSequence',
           'proteinColumn':'ProteinName',
           'transitionColumn':'transition_name',
           'uniprotColumn': 'UniprotID',
           'precursorColumn': 'PrecursorMZ'}

OSW_COLUMNS = {'peptideColumn' : 'FID',
           'proteinColumn':'ID',
           'transitionColumn':'ID',
           'precursorColumn': None}


def num_of_peptides(df, peptideColumn='PeptideSequence'):
    # only for data where on row represents a transtition,
    # not for pyprophet export
    return df.groupby([peptideColumn]).count().shape[0]


def num_of_proteins_group(df, proteinColumn = 'ProteinName'):
    return df.groupby([proteinColumn]).count().shape[0]


def num_of_proteins(df, proteinColumn = 'ProteinName'):
    temp_list = list(set(df[proteinColumn]))
    totallist = []
    for i in temp_list:
        if type(i) is int:
            continue
        else:
            sublist = i.split('/')[1:]
            for x in sublist:
                subsublist = x.split('|')
                if len(subsublist) < 2:
                    continue
                else:
                    totallist.append(subsublist[1])
    return len(set(totallist))


def num_of_transitions(df, transitionColumn='transition_name'):
    return df.groupby([transitionColumn]).count().shape[0]


def count_transitions_from_aggr(entrystr):
    temp_list = [i for i in entrystr.split(';') if float(i) != 0.0]
    return len(temp_list)


def num_of_transitions_from_aggr(df, aggrColumn='aggr_Peak_Area'):
    # for pyprophet tsv export
    return df[aggrColumn].apply(lambda x: count_transitions_from_aggr(x)).sum()


def num_of_precursors(df, precursorColumn='PrecursorMz'):
    return df.groupby([precursorColumn]).count().shape[0]


def library_summary_tsv(df):
    """
    the difference between library summary anf fiel_summary is, that library_sumary also counts the number of precursors
    :param df: library as pandas df
    :return: dictionary
    """
    summaryDict = {'ProteinGroups': num_of_proteins_group(df),
               'Proteins': num_of_proteins(df),
               'Peptides': num_of_peptides(df),
               'Precursors': num_of_precursors(df),
               'Transitions': num_of_transitions(df)}
    return summaryDict


def library_summary_pqp(df):
    transition = df['transition']
    peptide = df['peptide']
    protein = df['protein']
    #TODO make this better
    summaryDict = {'Proteins': len(list(set(protein[protein['DECOY'] == 0]['ID']))),
                'Peptides': len(list(set(peptide[peptide['DECOY'] == 0]['ID']))),
               'Transitions':  transition[transition['DECOY'] == 0]['ID'].count()}

    return summaryDict

def library_summary(df):
    if isinstance(df, dict):
        temp_dict =  library_summary_pqp(df)
    else:
        temp_dict = library_summary_tsv(df)

    lib_sum_df = pd.DataFrame(temp_dict, index=[0])
    #lib_sum_df = lib_sum_df.reindex(['Species', 'Proteins', 'ProteinGroups', 'Peptides', 'Precursors', 'Transitions'], axis=1)

    return lib_sum_df


def file_summary(df, columnNames=PYPROPHET_TSV_COLUMNS, onlyTargets=True, decoyColumn='decoy'):
    # file summary tsv (swath or pyprophet)
    if decoyColumn in df.columns or 'DECOY' in df.columns:
        if decoyColumn not in df.columns:
            decoyColumn = 'DECOY'
        target_df = df[df[decoyColumn] == 0]
        decoy_df = df[df[decoyColumn] == 1]
        if onlyTargets:
            summaryDict = {'ProteinGroups': num_of_proteins_group(target_df, columnNames['proteinColumn']),
                           'Proteins': num_of_proteins(target_df, columnNames['proteinColumn']),
                           'Peptides': num_of_peptides(target_df, columnNames['peptideColumn']),
                           'Transitions': num_of_transitions_from_aggr(target_df)}
        else:
            summaryDict = {'ProteinGroups': num_of_proteins_group(target_df, columnNames['proteinColumn']),
                           'Proteins': num_of_proteins(target_df, columnNames['proteinColumn']),
                           'Peptides': num_of_peptides(target_df, columnNames['peptideColumn']),
                           'Transitions': num_of_transitions(target_df, columnNames['transitionColumn']),
                           'DecoyProteinGroups': num_of_proteins_group(decoy_df, columnNames['proteinColumn']),
                           'DecoyProteins': num_of_proteins(decoy_df, columnNames['proteinColumn']),
                           'DecoyPeptides': num_of_peptides(decoy_df, columnNames['peptideColumn']),
                           'DecoyTransitions': num_of_transitions_from_aggr(decoy_df, columnNames['transitionColumn'])
                           }

    else:
        summaryDict = {'ProteinGroups': num_of_proteins_group(df, columnNames['proteinColumn']),
                       'Proteins': num_of_proteins(df, columnNames['proteinColumn']),
                       'Peptides': num_of_peptides(df, columnNames['peptideColumn']),
                       'Transitions': num_of_transitions_from_aggr(df, columnNames['transitionColumn'])}

    return summaryDict


def file_summary_osw(dfpeptide, dfprotein, dftransition, onlyTargets,  columnNames=OSW_COLUMNS, decoyColumn='DECOY'):

    if decoyColumn in dfpeptide.columns or 'DECOY' in dfpeptide.columns:
        if decoyColumn not in dfpeptide.columns:
            decoyColumn = 'decoy'
        target_df_prot = dfprotein[dfprotein[decoyColumn] == 0]
        decoy_df_prot = dfprotein[dfprotein[decoyColumn] == 1]
        target_df_peptide = dfpeptide[dfpeptide[decoyColumn] == 0]
        decoy_df_peptide = dfpeptide[dfpeptide[decoyColumn] == 1]
        #target_df_transition = dftransition[dftransition[decoyColumn] == 0]
        #decoy_df_transition = dftransition[dftransition[decoyColumn] == 1]

        if onlyTargets:
            summaryDict = {'ProteinGroups': num_of_proteins_group(target_df_prot, columnNames['proteinColumn']),
                           'Proteins': num_of_proteins(target_df_prot, columnNames['proteinColumn']),
                           'Peptides': num_of_peptides(target_df_peptide, columnNames['peptideColumn']),
                           'Transitions': num_of_transitions(dftransition, columnNames['transitionColumn'])}
        else:
            summaryDict = {'ProteinGroups': num_of_proteins_group(target_df_prot, columnNames['proteinColumn']),
                           'Proteins': num_of_proteins(target_df_prot, columnNames['proteinColumn']),
                           'Peptides': num_of_peptides(target_df_peptide, columnNames['peptideColumn']),
                           'Transitions': num_of_transitions(target_df_peptide, columnNames['transitionColumn']),
                           'DecoyProteinGroups': num_of_proteins_group(decoy_df_prot, columnNames['proteinColumn']),
                           'DecoyProteins': num_of_proteins(decoy_df_prot, columnNames['proteinColumn']),
                           'DecoyPeptides': num_of_peptides(decoy_df_peptide, columnNames['peptideColumn']),
                           'DecoyTransitions': num_of_transitions(dftransition, columnNames['transitionColumn'])
                           }

    else:
        summaryDict = {'ProteinGroups': num_of_proteins_group(dfprotein, columnNames['proteinColumn']),
                       'Proteins': num_of_proteins(dfprotein, columnNames['proteinColumn']),
                       'Peptides': num_of_peptides(dfpeptide, columnNames['peptideColumn']),
                       'Transitions': num_of_transitions(dfpeptide, columnNames['transitionColumn'])}

    return summaryDict



def library_summary_by_species(df, proteinColumn='ProteinName', speciesList=['HUMAN', 'YEAS', 'ECOLI']):
    """

    :param df:
    :param proteinColumn:
    :param speciesList: list of species as they are called in proteinColumns
    :return: dataframe
    """
    listOfDicts = []
    totalSummaryDict = library_summary_tsv(df)
    totalSummaryDict['Species'] = 'TOTAL'
    listOfDicts.append(totalSummaryDict)

    for i in speciesList:
        # subset dataframe:
        df[i] = df[proteinColumn].str.contains(i).astype(int)
        temp_df = df[df[i] == 1]
        summaryDict = library_summary_tsv(temp_df)
        summaryDict['Species'] = i
        listOfDicts.append(summaryDict)

    df = pd.DataFrame(listOfDicts)
    df = df.reindex(['Species', 'Proteins', 'ProteinGroups', 'Peptides', 'Precursors', 'Transitions'], axis=1)
    return df


def file_summary_by_species(df, columnNames, speciesList=['HUMAN', 'YEAS', 'ECOLI']):
    """

    :param df:
    :param columnNames:
    :param speciesList:
    :return: df
    """
    proteinColumn=columnNames['proteinColumn']
    listOfDicts = []
    totalSummaryDict = file_summary(df, columnNames)
    totalSummaryDict['Species'] = 'TOTAL'
    listOfDicts.append(totalSummaryDict)

    for i in speciesList:
        # subset dataframe:
        df[i] = df[proteinColumn].str.contains(i).astype(int)
        temp_df = df[df[i] == 1]
        summaryDict = file_summary(temp_df, columnNames)
        summaryDict['Species'] = i
        listOfDicts.append(summaryDict)

    df = pd.DataFrame(listOfDicts)
    df = df.reindex(['Species', 'Proteins', 'ProteinGroups', 'Peptides',  'Transitions'], axis=1)
    return df


def input_files_summary_table(dfdict, onlyTargets=True):

    listOfDicts = []
    for key, dfs in dfdict.items():
        if key.lower().endswith('tsv'):
            temp_dict = file_summary(dfs, PYPROPHET_TSV_COLUMNS, onlyTargets)
            temp_dict['File'] = key
            listOfDicts.append(temp_dict)
        if key.lower().endswith(('osw', 'pqp')):
            dfs = dfdict[key]
            transition = dfs['featureTransition']
            peptide = dfs['peptide']
            protein = dfs['protein']
            temp_dict = file_summary_osw(protein, peptide, transition, OSW_COLUMNS, onlyTargets)
            temp_dict['File'] = key
            listOfDicts.append(temp_dict)


    columnsTitlesDecoy = ['File', 'Proteins', 'ProteinGroups', 'Peptides', 'Transitions',
                     'DecoyProteinGroups', 'DecoyProteins', 'DecoyPeptides', 'DecoyTransitions']
    columnsTitles= ['File', 'Proteins', 'ProteinGroups', 'Peptides', 'Transitions']

    table = pd.DataFrame(listOfDicts)
    if 'DecoyProteins' in table.columns:
        table = table.reindex(columns=columnsTitlesDecoy)
    else:
        table = table.reindex(columns=columnsTitles)

    return table.to_html()




