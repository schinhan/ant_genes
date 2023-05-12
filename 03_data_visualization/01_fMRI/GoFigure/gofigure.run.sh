#!/bin/bash

# GO-Figure clustering: summary visualisation of list of GO terms
# purpose: summarize results of gene set enrichment analysis
# for detailed explanations and instructions see [https://gitlab.com/evogenlab/GO-Figure]
# use standard input consisting of a tabular format file with two columns: GO term, P-value

# requirements: python 3 with the following packages:
# numpy, matplotlib, seaborn, scikit-learn, adjustText 

# set directory
cd ../../03_data_visualization/01_fMRI/GoFigure

#clone gitlab project
git clone https://gitlab.com/evogenlab/GO-Figure.git

# run programm for all ANT conditions (input: all FDR-adjusted GO terms)
# settings: max. 20 clusters displayed, similarity cutoff = 0.5

# alerting
python gofigure.py -i gofigure.input.fdr.alert.txt -o alerting_results -a 20 -si 0.5  -c 'log10-pval' -e 100 -p Blues_r  -k 'xx-large' -t "Alerting" 

# control
python gofigure.py -i gofigure.input.fdr.control.txt -o control_results -a 20 -si 0.5  -c 'log10-pval' -e 100 -p Blues_r -k 'xx-large' -t "Control"

# orienting  
python gofigure.py -i gofigure.input.fdr.orient.txt -o orient_results -a 20 -si 0.5  -c 'log10-pval' -e 100 -p Blues_r -k 'xx-large' -t "Orienting"

# use these visualisations or create your own figures with the output files