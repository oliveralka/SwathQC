'''reading sqlite based data fromats,
including:
    - .osw
    - .pqp
read in only the necessary information for plotting and summaries
'''
import sqlite3
import pandas as pd


def conn_osw(osw_file):
    """connect to sqlite based database (osw or pqp) and return connection
    :param osw file
    :return conn, database connection object"""
    conn = sqlite3.connect(osw_file)
    return conn


def list_osw_tables(conn):
    """list all available tables in database
    :param conn: connection to databse use conn_osw)
    :return tables: list of table names as strings"""
    conn_cursor = conn.cursor()
    conn.text_factory = str
    res = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = [name[0] for name in res]
    return tables


def list_osw_columns(conn, tables):
    """ get set of column names in database
    :param conn: databse connection
    :param tables: list of tables in database
    :return: sorted list of all column names
    """
    names = []
    for t in tables:
        res = conn.execute("select * from " + t)
        for description in res.description:
            names.append(description[0])
    return sorted(list(set(names)))


def column_in_osw(osw_file, col):
    """ check if a column is present in any table of the dataframe
    :param osw_file:
    :param col:
    :return:
    """
    conn = conn_osw(osw_file)
    tables = list_osw_tables(conn)
    cols = list_osw_columns(conn, tables)
    return col in cols


def osw_table_to_df(osw_file, table_name):
    """
    read a table into pandas dataframe
    :param osw_file:
    :param table_name:
    :return:
    """
    conn = conn_osw(osw_file)
    df = pd.read_sql_query("SELECT * FROM " + table_name, conn)
    conn.close()
    return df


def query_osw_to_df(file, query):
    """

    :param file: osw file
    :param query: string of sql query
    :return: dataframe
    """
    conn = conn_osw(file)
    df = pd.read_sql_query(query, conn)
    return df


def rt_from_osw(file):
    """
    from pyprophet osw file
    :param file:
    :return:
    """
    conn = conn_osw(file)
    df = pd.read_sql_query("""SELECT FEATURE.NORM_RT as iRT, 
                                     FEATURE.EXP_RT as RT,
                                     FEATURE.DELTA_RT as delta_rt,
                                     FEATURE.PRECURSOR_ID as precursor, 
                                     PRECURSOR.LIBRARY_RT, 
                                     PRECURSOR.DECOY as decoy
                                     FROM PRECURSOR
                                     INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID """, conn)
    return df


def pyprophet_peptide_score_from_osw(file):
    """read in peptide information from two tables"""
    conn = conn_osw(file)
    df = pd.read_sql_query("""SELECT SCORE_PEPTIDE.PEPTIDE_ID as ID,
                            SCORE_PEPTIDE.SCORE as d_score,
                            SCORE_PEPTIDE.QVALUE,
                            SCORE_PEPTIDE.PVALUE,
                            SCORE_PEPTIDE.PEP,
                            SCORE_PEPTIDE.CONTEXT,
                            SCORE_PEPTIDE.RUN_ID,
                            PEPTIDE.ID as precursor, 
                            PEPTIDE.DECOY as decoy
                            FROM PEPTIDE
                            INNER JOIN SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PEPTIDE.ID """, conn)
    return df


def pyprophet_protein_score_from_osw(file):
    """read in protein information from two tables"""
    conn = conn_osw(file)
    df = pd.read_sql_query("""SELECT SCORE_PROTEIN.PROTEIN_ID,
                                SCORE_PROTEIN.SCORE as d_score,
                                SCORE_PROTEIN.QVALUE,
                                SCORE_PROTEIN.PVALUE,
                                SCORE_PROTEIN.PEP,
                                SCORE_PROTEIN.CONTEXT,
                                SCORE_PROTEIN.RUN_ID,
                                PROTEIN.ID, 
                                PROTEIN.DECOY as decoy
                            FROM PROTEIN
                            INNER JOIN SCORE_PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = PROTEIN.ID  """, conn)
    return df


def feature_tables_from_osw(file):
    """Information from FEATURE and FEATURE_MS2 for plotting"""
    conn = conn_osw(file)
    df = pd.read_sql_query("""SELECT FEATURE.ID, 
                            FEATURE.NORM_RT as normRT, 
                            FEATURE.EXP_RT as RT, 
                            FEATURE.DELTA_RT as delta_RT, 
                            FEATURE.LEFT_WIDTH leftWidth, 
                            FEATURE.RIGHT_WIDTH as rightWidth, 
                            FEATURE_MS2.FEATURE_ID, 
                            FEATURE_MS2.VAR_MASSDEV_SCORE,
                            FEATURE_MS2.VAR_MASSDEV_SCORE_WEIGHTED, 
                            FEATURE_MS2.VAR_LIBRARY_CORR
                        FROM FEATURE 
                        INNER JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID """, conn)
    return df


def feature_transition_table_from_osw(file):
    """Information from FEATURE_TRANSITION for plotting"""
    conn = conn_osw(file)
    df = pd.read_sql_query("""SELECT TRANSITION_ID, 
                            FEATURE_ID, 
                            AREA_INTENSITY,
                            APEX_INTENSITY
                        FROM FEATURE_TRANSITION""", conn)
    return df


def swath_osw_to_dfdict(file):
    subdict = {'feature': feature_tables_from_osw(file),
               'featureTransition': feature_transition_table_from_osw(file),
               'peptide': osw_table_to_df(file, 'PEPTIDE'),
               'protein': osw_table_to_df(file, 'PROTEIN')}
    return subdict


def pp_osw_to_dfdict(file):
    subdict = {'feature': feature_tables_from_osw(file),
               'featureTransition': feature_transition_table_from_osw(file),
               'peptide': osw_table_to_df(file, 'PEPTIDE'),
               'protein': osw_table_to_df(file, 'PROTEIN'),
               'rt': rt_from_osw(file),
               'peptideScores': pyprophet_peptide_score_from_osw(file),
               'proteinScores': pyprophet_protein_score_from_osw(file)}
    return subdict



