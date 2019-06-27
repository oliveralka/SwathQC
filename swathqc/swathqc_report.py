""" run swathqc_report.py to generate formatted HTML report """

from __future__ import print_function
import argparse
import sys
import os
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from swathqc import *

PATH = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_ENVIRONMENT = Environment(
    autoescape=False,
    loader=FileSystemLoader(os.path.join(PATH, 'templates')),
    trim_blocks=False)


def render_template(template_filename, content, title):
    # build navigation
    navigation = [(i['name'], i['header']) for i in content]
    return TEMPLATE_ENVIRONMENT.get_template(template_filename).render(sections=content, title=title,
                                                                       navigation=navigation)


def prepare_section(section_id, header=None, paragraph=None, images=None, table=None):
    section_dict = {
        'name': section_id,
        'header': header,
        'images': images,
        'table': table,
        'paragraph': paragraph
    }

    return section_dict


def main():
    def valid_file(choices, fname):
        """check file extension and raise parser error"""
        ext = os.path.splitext(fname)[1][1:]
        if ext.lower() not in choices:
            parser.error("file doesn't end with one of {}".format(choices))
        return fname

    # add command line arguments
    parser = argparse.ArgumentParser(description="Create HTML report for openSwathWorkflow, pyprophet and tric output. "
                                                 "This is intended for (but not limited to) the sqlite based workflow.")
    parser.add_argument('-s', help="swath output, valid formats are osw or tsv,", dest='swath_files', nargs='+',
                        type=lambda s: valid_file(("tsv", "osw"), s))  # eiuther featureXML, tsv or osw
    parser.add_argument('-p', help="output from pyprophet, recommended is tsv, valid are tsv and osw",
                        dest='pp_file', nargs='+', type=lambda s: valid_file(("tsv", "osw"), s))
    parser.add_argument('-ppdscores', help="pp osw files for plotting dscores",
                        dest='ppscores', nargs='+', type=lambda s: valid_file(("osw"), s))
    parser.add_argument('-t', help="tric tsv files", dest='tric_files', nargs='+',
                        type=lambda s: valid_file(("tsv"), s), required=False)
    # parser.add_argument('-j', "--qc_json", help="input in json format", dest='input_json', type=lambda s: valid_file(("json"), s))
    parser.add_argument('-l', help='library file in tsv or pqp format', dest='lib',
                        type=lambda s: valid_file(("tsv", "pqp"), s), required=False)
    parser.add_argument('-o', help='outfile name, default: "report.html"', default='report.html', dest='outfile',
                        type=str, required=False)
    parser.add_argument('-title', help='add Title', dest='title',
                        type=str, required=False)
    # parse arguments
    args = parser.parse_args()
    # print help and exit if no argument gets passed
    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()

    # get the system path separator, (flexibility for windows and linux)
    separator = os.path.sep
    # dfdict = {'}
    # keep track of read in files
    filedict = {}
    # read in files
    dictlist = []

    swathDict = {}
    if args.swath_files is not None:
        swathFiles = []
        for i in args.swath_files:
            f = str.split(i, separator)[-1]
            swathFiles.append(f)
            # read available files into dataframes
            if f.endswith("featureXML"):
                continue
            if f.endswith("tsv"):
                swathDict[f] = pd.read_csv(i, sep='\t')
            if f.endswith("osw"):
                # read all necessary tables into dictionary of dataframes:
                # {'file1.osw': {'feature': 'featuretable1', 'peptide': 'peptidetable1'},
                # 'file2.osw': {'feature': 'featuretable2', 'peptide': 'peptidetable2'}}
                subdict = {'feature': feature_tables_from_osw(i),
                           'featureTransition': feature_transition_table_from_osw(i),
                           'peptide': osw_table_to_df(i, 'PEPTIDE'),
                           'protein': osw_table_to_df(i, 'PROTEIN')}
                swathDict[f] = subdict
        filedict['swath files'] = swathFiles
        dictlist.append(swathDict)

    ppDict = {}
    if args.pp_file is not None:
        ppFiles = []
        for ppFile in args.pp_file:
            ppfile = str.split(ppFile, separator)[-1]
            ppFiles.append(ppfile)
            if ppfile.endswith('tsv'):
                pp = pd.read_csv(ppFile, sep='\t')
                ppDict[ppfile] = pp
            if ppfile.endswith('osw'):
                subdict = {'feature': feature_tables_from_osw(ppFile),
                           'featureTransition': feature_transition_table_from_osw(ppFile),
                           'peptide': osw_table_to_df(ppFile, 'PEPTIDE'),
                           'protein': osw_table_to_df(ppFile, 'PROTEIN'),
                           'rt': rt_from_osw(ppFile)}
                ppDict[ppfile] = subdict
        filedict['PyProphet files'] = ppFiles
        dictlist.append(ppDict)

    tricDict = {}
    if args.tric_files is not None:
        tricFiles = []
        for tricFile in args.tric_files:
            tricF = str.split(tricFile, separator)[-1]
            tricFiles.append(tricF)
            tricDict[tricF] = pd.read_csv(tricFile, sep="\t")
        filedict['TRIC files'] = tricFiles
        dictlist.append(tricDict)

    if args.lib is not None:
        libfile = str.split(args.lib, separator)[-1]
        filedict['library'] = [libfile]
        if libfile.lower().endswith('tsv'):
            lib = pd.read_csv(args.lib, sep='\t')
        if libfile.lower().endswith('pqp'):
            subdict = {'peptide': osw_table_to_df(args.lib, 'PEPTIDE'),
                       'protein': osw_table_to_df(args.lib, 'PROTEIN'),
                       'transition': osw_table_to_df(args.lib, 'TRANSITION')}
            lib = subdict
    # make sections with plots to pass to html template
    # list of dicts

    if args.swath_files is not None:
        min_length = len(args.swath_files)
    else:
        if args.pp_file is not None:
            min_length = len(args.pp_file)
        else:
            min_length = 1

    if args.lib is not None:
        lib_summary = ('Summary of assay Library', library_summary(lib).to_html())
    else:
        lib_summary = ('Library Summary', 'no assay library provided')
    dfdict = dict(j for i in dictlist for j in i.items())
    filetable = pd.DataFrame({k: pd.Series(v[:min_length]) for k, v in filedict.items()}).fillna('-')

    sections = [prepare_section('top', header="Info", paragraph='',
                                table=[('Input files', filetable.to_html())
                                       #('Summary of identified target peptides, proteins and transitions in DIA input files',
                                       #input_files_summary_table(dfdict))
                                       ]),
                prepare_section('idoverrt', header="ID over RT",
                                images=[report_id_over_rt(swathDict), report_id_over_rt(ppDict),
                                        report_id_over_rt(tricDict, color="#FFB350")]),
                prepare_section('numtrans', header='#Transitions',
                                images=[report_num_transitions(swathDict, color="#3A78A4"),  #
                                        report_num_transitions(ppDict),
                                        report_num_transitions(tricDict, color="#FFB350")]),
                prepare_section('libintcorr', header='Library Intensity', images=[report_lib_intensity_corr(swathDict),
                                                                                  report_lib_intensity_corr(ppDict)]),
                prepare_section('pwoverrt', header="PW over RT",
                                images=[report_pw_over_rt(swathDict), report_pw_over_rt(ppDict)]),

                prepare_section('rtcorr', header='RT correlation', images=report_rt_correlation(ppDict))

                ]
    if args.swath_files is not None:
        if args.swath_files[0].endswith('tsv'):
            sections.append(
                prepare_section('masserror', header='Mass Error', images=report_mass_error(swathDict, color="#3A78A4")))
            sections.append(prepare_section('mc', header='Missed Cleavages',
                                            images=report_missed_cleavages(swathDict, color="#3A78A4")))
    if args.lib is not None:
        print('Library was provided')
        sections.append(
            prepare_section('libcov', header='Library Coverage', images=[report_libray_coverage(lib, ppDict)]))

    if args.ppscores is not None:
        scoreDict = {}
        for scoreFile in args.ppscores:
            scorefile = str.split(scoreFile, separator)[-1]
            ppFiles.append(ppfile)
            peptidescores = pyprophet_peptide_score_from_osw(scoreFile)
            proteinscores = pyprophet_protein_score_from_osw(scoreFile)
            subdict = {'peptideScores': peptidescores, 'proteinScores': proteinscores}
            scoreDict[scorefile] = subdict

        sections.append(prepare_section('dscores', header='d-Scores', images=report_dscore(scoreDict)))

    #   if args.ppscores is None:
    #        sections.append(prepare_section('dscores', header='d-Scores', images=report_dscore(ppDict)))

    # env = Environment(loader=FileSystemLoader('templates/'))  # or parh to template location
    # template_vars = temp_vars

    # write out html string into report html file:
    html_out = render_template('base.j2', sections, args.title)
    # weasy print rewires some more packages, GTK+
    # HTML(html_out).write_pdf(args.outfile, stlysheets=["templates/style.css"])
    with open(args.outfile, "w") as file:
        file.write(html_out)


if __name__ == "__main__":
    main()

# try:
#     import matplotlib
#     matplotlib.use('Agg')
#     from matplotlib.backends.backend_pdf import PdfPages
#     import matplotlib.pyplot as plt
# except ImportError:
#     plt = None
