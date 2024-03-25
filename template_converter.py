# Converting Michael Brown's template format to that used by EAZY
# ------------------------------------------------------------

import pandas as pd
import glob
import os
import numpy as np

# Set the directory
in_dir = 'G:/AGN/hlsp_agnsedatlas_multi_multi_all_multi_v1_collection/templates_observed/'
out_dir = 'C:/Users/eddie/PycharmProjects/SEDTemplate_conda/templates/hlsp_agnsedatlas_observed/'

# Get the txt files  in the directory
txt_files = glob.glob(in_dir + '*.txt')

for txt_files in txt_files:

    # Read the txt file
    data = pd.read_csv(txt_files, sep=" ", comment='#', header=None, skipinitialspace=True)

    #scale_fix = np.multiply(data[2], data[3])
    output_data = data[[0,2]]

    # Get the base name of the txt file
    base_name = os.path.basename(os.path.splitext(txt_files)[0])

    # Write the selected data to a dat file
    output_data.to_csv(out_dir + f'{base_name}.dat', sep=' ', index=False, header=False)
