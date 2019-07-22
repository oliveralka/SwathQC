# SWATH QC 

SwathQC is a python tool which can be used for extraction and visualization of 
QC data für SWATH analysis with OpenSWATH. 

It can either be used as a command line tool to generate a html report, or
used imported as an external library python scripts.
An example jupyter notebook is provided on github.

To generate a report run swathqc reprot.py with the following options: <br/>
-s SWATH analysis data either osw or tsv format <br/>
-p FDR filtered export from PyProphet in tsv format <br/>
-ppdscores merged osw from PyProphet for plotting dscores <br/>
-t TRIC aligned output in tsv format <br/>
-l library assay library as tsv or pqp format <br/>
-title title on top o the HTML report <br/>
-o output filename, default: report.html <br/>
