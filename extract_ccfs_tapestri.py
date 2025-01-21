# Import mosaic libraries
import missionbio.mosaic as ms

# Import these to display entire dataframes
from IPython.display import display, HTML

# Import graph_objects from the plotly package to display figures when saving the notebook as an HTML
import plotly as px
import plotly.graph_objects as go

# Import additional packages for specific visuals
import matplotlib.pyplot as plt
import plotly.offline as pyo
# pyo.init_notebook_mode()
import numpy as np
import seaborn as sns
from missionbio.plotting.multimap import MultiMap
from missionbio.plotting.heatmap import Heatmap

# Import COMPASS for imputation
from missionbio.mosaic.algorithms.compass import COMPASS
import pandas as pd
import os
import glob

file_xlxs = "/orfeo/LTS/CDSLab/LT_storage/CLL/scDNA/h5_file/Tapestry CLL Targeted Panel.xlsx"
dfs = pd.read_excel(file_xlxs, sheet_name=None)
wgs_muts = dfs["X_Tapestry"]

# h5 files to be read 
h5_path = "/orfeo/LTS/CDSLab/LT_storage/CLL/scDNA/h5_file/"
files = glob.glob(h5_path+"*.h5")


# maxi function to extract data and get the ccfs of specific mutations
def get_portion_cells(h5path, # path of the h5 files
                      patient, # patient
                      sample, # sample
                      ids_file, # table with the mutations and ccf from bulk 
                      whitelist = True, # should we use a whitelist? 
                      filt = False # should we filter the variants?
                     ): 
    patient_id = patient
    sample_id = sample
    cff_name = "CCF_" + sample_id
    muts = ids_file[(ids_file["Case"]==patient_id) & (ids_file["TAPESTRI"]=="YES") & (ids_file[cff_name]!=0)]
    muts["Tapestri_id"] = muts[['Chr','Start','Ref','Var']].apply(lambda x : '{}:{}:{}/{}'.format(x[0],x[1],x[2],x[3]), axis=1)

    if whitelist == True:
        whitelist = list(set(list(muts[['Chr','Start','Ref','Var']].apply(lambda x : '{}:{}:{}/{}'.format(x[0],x[1],x[2],x[3]), axis=1))))
    else:
        whitelist = []
    # Specify the h5 file to be used in this analysis: h5path = '/path/to/h5/file/test.h5'
    # load one sample for testing

    # Load the data
    sample = ms.load(h5path, 
                     raw=False, 
                     filter_variants=filt, 
                     single=True, 
                     whitelist = whitelist, 
                     filter_cells=False)
    
    vars_kept = list(sample.dna.filter_variants())
    vars_kept = list(set(vars_kept + whitelist))

    # filter the variants
    sample_dna_filt = sample.dna[sample.dna.barcodes(), vars_kept]

    col_attributes = pd.DataFrame(sample_dna_filt.col_attrs)
    # col_attributes

    NGT = sample_dna_filt[:,vars_kept].layers["NGT"]

    # create the table
    # Iterate over columns of the array 
    
    result = {}

    # get the total number of cells 
    tot_cells = NGT.shape[0]
    
    # get the mutations 
    mutations = sample_dna_filt[:,vars_kept].ids()
    
    for id_indx, mut in enumerate(mutations):
        col = NGT[:, id_indx]  # Extract column
        unique, counts = np.unique(col, return_counts=True)
    
        props = dict(zip(unique, counts/tot_cells))
        # print(props)
        props.setdefault(2, 0)
        props.setdefault(1, 0)
           
        result[mut] = props[1]+props[2]
        
    CCF_table = pd.DataFrame.from_dict(result, orient="index", columns=["CCF_" + patient_id + "_" + sample_id])
    return CCF_table

# example run
df = get_portion_cells("/orfeo/LTS/CDSLab/LT_storage/CLL/scDNA/h5_file/CT339_T2_POS_Test_11.dna.h5", 
                  patient = "CT339", 
                  sample = "T2_pos", # sample
                  ids_file = wgs_muts, # table with the mutations and ccf from bulk 
                  whitelist = True, # should we use a whitelist? 
                  filt = False)