# SWATH QC 

SwathQC is a python tool which can be used for extraction and visualization of 
QC data für SWATH analysis with OpenSWATH. 


It can either be used as a command line tool to generate a html report, or
used imported as an external library python scripts. The second option oers
more exibility. 

An example jupyter notebook is provided on github.

To generate a report run swathqc reprot.py with the following options:
-s SWATH analysis data either osw or tsv format
-p FDR filtered export from PyProphet in tsv format
-ppdscores merged osw from PyProphet for plotting dscores
-t TRIC aligned output in tsv format
-l library assay library as tsv or pqp format
-title title on top o the HTML report
-o output filename, default: report.html
