# create gene expression atlas based on Lausanne Atlas
# for detailed instructions see:
# [https://abagen.readthedocs.io/en/stable/index.html]

# requirements: python 3.6+ with the abagen python package
# you can also clone the corresponding gitlab project via:
# git clone https://github.com/rmarkello/abagen.git

# set directory to where your brain atlas is
import os
os.chdir('../../02_data/02_parcellation')  

# fetch AHBA microarray data
import abagen
files = abagen.fetch_microarray(donors='all', verbose=0)

# define a parcellation 
# make sure atlas and information file are added to path
atlas = {'image': ['Lausanne.250.L.10k_fsavg_5.label.gii', 'Lausanne.250.R.10k_fsavg_5.label.gii'],
        'info': 'atlas_laus250.csv'
}

# parcellate expression data
expression = abagen.get_expression_data(atlas['image'], atlas['info'])
print (expression)

#save expression data to csv
import pandas as pd
expression.to_csv('lausanne_parc_atlas.csv',  sep=',')
