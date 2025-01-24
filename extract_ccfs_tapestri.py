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

samples = ["CT287","TS187","CT339","CT344","CT48","CT525","CZ40"] # ,"RM238"] removed bc not sure which time point the h5 file refers to 

# prepare the input for the extraction of data
metadata = {}

for i in samples:
    file = glob.glob(h5_path+"*"+i+"*")
    for n in range(len(file)):
        time_point = re.findall("_T[1-3]", file[n], re.IGNORECASE)
        sample_type = re.findall("neg|pos", file[n], re.IGNORECASE)
        
        if((len(time_point) > 0) and (len(sample_type) > 0)):
            time_point = re.sub("_", "", time_point[0])
            sample_type = sample_type[0]
        
            res = {"path" : file[n], 
                "patient" : i, 
                "sample" : time_point + "_" + sample_type.lower()}

            metadata[i + "_" + time_point + "_" + sample_type.lower()] = res



# maxi function to extract data and get the ccfs of specific mutations
def get_portion_cells(metadata_dict, # dictionary with information on patient, sample and path of the h5 file
                      ids_file, # table with the mutations and ccf from bulk 
                      whitelist = True, # should we use a whitelist? 
                      filt = False # should we filter the variants?
                     ): 
    patient_id = metadata_dict["patient"]
    sample_id = metadata_dict["sample"]
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
    sample = ms.load(metadata_dict["path"], 
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

    # add VAF? report also them in the final table
        
    CCF_table = pd.DataFrame.from_dict(result, orient="index", columns=["Tapestri_CCF_" + patient_id + "_" + sample_id])
    CCF_table['mutation'] = CCF_table.index

    bulk_reduced = muts[["Tapestri_id", cff_name]]
    bulk_reduced.rename(columns = {cff_name : "Bulk_" + cff_name})

    CCF_table_final = pd.merge(left=CCF_table, right=bulk_reduced, left_on='mutation', right_on='Tapestri_id')
    #CCF_table_final = CCF_table_final[['mutation', "Tapestri_CCF_" + patient_id + "_" + sample_id, "Bulk_" + cff_name]]
    return CCF_table_final

# produce the data
os.mkdir("excel_positions")
results_ccf = {}

for i in metadata.keys():
    results_ccf[i] = get_portion_cells(metadata_dict = metadata[i],
                       ids_file = wgs_muts, # table with the mutations and ccf from bulk 
                       whitelist = True, # should we use a whitelist? 
                       filt = False)
    results_ccf[i].to_csv("excel_positions/" + i + "_comparison_tapestri_vs_bulk_excel_positions.csv")