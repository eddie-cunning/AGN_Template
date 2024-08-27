# will run each template individually
# each test needs a new test name, otherwise ti will append to the old name directory and main file.
# change the key, field, zstep for each test
# for each EAZY run, the following parameters will be returned:
"""
    - object ID
    - template
    - NMAD
    - outlier count
    The following will be in a new file as they are per object:
    - flux density model
    - phot redshift
    - chi2
"""
# ----------------------------------------------------------------------------------------------------------------------
import os
import glob
import matplotlib.pyplot as plt
import warnings
import numpy as np
import pandas as pd
from astropy.utils.exceptions import AstropyWarning
import time
import importlib
import sys
import eazy
import astropy.stats
import multiprocessing as mp

np.seterr(all='ignore')
warnings.simplefilter('ignore', category=AstropyWarning)
# print('EAZYCODE = '+ str(os.getenv('EAZYCODE')) + '\n')

# print(time.ctime() + '\n')

# print(sys.version + '\n')

# for module in ['numpy', 'scipy', 'matplotlib','astropy','eazy']:
#    mod = importlib.import_module(module)
#    print('{0:>20} : {1}'.format(module, mod.__version__))
# ----------------------------------------------------------------------------------------------------------------------
os.getcwd()

# Load ZFOURGE catalogue from local drive
test_title = 'individual_bin0.3to0.4_0.05'  # title of the test, eg. 1,2, A, B, Initial.
field = 'uds'  # 'cdfs', 'cosmos', or 'uds'

# Choose ID key for the catalogue
# A key designates what object you want included
id_key = 'fraction0.3to0.4'

# Directories for key, name keys anything, just to keep track of any complex object choices made in catalogue_prepare.ipynb
id_key_dict = {
    'normal': f'inputs/alternate_catalogues/{field}.normal.cat',
    'fraction0.0to0.1': f'inputs/alternate_catalogues/{field}.fraction.bin0.0to0.1.cat',
    'fraction0.1to0.2': f'inputs/alternate_catalogues/{field}.fraction.bin0.1to0.2.cat',
    'fraction0.2to0.3': f'inputs/alternate_catalogues/{field}.fraction.bin0.2to0.3.cat',
    'fraction0.3to0.4': f'inputs/alternate_catalogues/{field}.fraction.bin0.3to0.4.cat',
    'fraction0.4to0.5': f'inputs/alternate_catalogues/{field}.fraction.bin0.4to0.5.cat',
    'ir_agn': f'inputs/alternate_catalogues/{field}.ir_agn.cat',
    'radio_agn': f'inputs/alternate_catalogues/{field}.radio_agn.cat',
    'xray_agn': f'inputs/alternate_catalogues/{field}.xray_agn.cat',
    'only_agn_0.4': f'inputs/alternate_catalogues/{field}.only_agn_above_0.4.cat',
    'only_agn_0.5': f'inputs/alternate_catalogues/{field}.only_agn_above_0.5.cat',
    'lacy': f'inputs/alternate_catalogues/{field}.lacy_wedge.cat',
    'donley': f'inputs/alternate_catalogues/{field}.donley_wedge.cat',
    'useflag': f'inputs/alternate_catalogues/{field}.useflag.cat'
}

# AGN templates allocation
loop_number = 1  # what loop you are on

use_galaxy_templates = True  # set to True to use galaxy templates as well

# ----------------------------------------------------------------------------------------------------------------------
# Directories
output_location = 'G:/honours/outputs'
if not os.path.isdir(f'{output_location}/{field}/{test_title}'):
    os.makedirs(f'{output_location}/{field}/{test_title}')

# Where to Save number data

key_data_file = f'{output_location}/{field}/{test_title}/{test_title}_data.csv'
headings = ['id_key', 'zstep', 'loop_number', 'agn_templates', 'galaxy templates', 'total_obj', 'mean_agn_frac',
            'spec_count', 'outlier_count', 'nmad_val']

key_data = pd.DataFrame(columns=headings)
if not os.path.isfile(key_data_file):
    key_data.to_csv(key_data_file, index=False)

# ----------------------------------------------------------------------------------------------------------------------
# Read the Catalogue

main_cat = pd.read_csv(id_key_dict[id_key])  # get the catalogue for the id_key
main_cat.to_csv('inputs/eazy_auto.cat',
                index=False)  # create a new catalogue, allows for change to be made in this cell

# Setting up the main catalogue
main = pd.read_csv('inputs/eazy_auto.cat', sep=" ", comment="#", header=None,
                   skipinitialspace=True)  # opening cut cat, and adjusting it
headers = pd.read_csv('inputs/eazy_auto.cat', sep=" ", header=None, nrows=1).iloc[0]
headers = headers[1:]
main.columns = headers

total_count = len(main)  # all objects in the range

# ----------------------------------------------------------------------------------------------------------------------
# Load any templates from the AGN template library

temp_param = 'templates/eazy_v1.3.spectra.param'  # basic parameter file, no agn templates
last_id = 9  # last id in the parameter file
empty_param = 'templates/eazy_v1.3_empty.param'  # empty parameter file
agn_param = 'templates/eazy_v1.3_AGN_auto.param'  # parameter file with agn templates

# opening the parameter files, and reading the contents
with open(temp_param) as f:
    original_galaxy = f.read()

with open(empty_param) as f:
    original_empty = f.read()

agn_dir = 'templates/hlsp_agnsedatlas_rest/'  # dir with all agn templates
agn_temp_all = os.listdir(agn_dir)


def agn_template_loader(templates, use_galaxy_templates=False):

    """
    Function to load AGN templates to the parameter file
    templates: list of templates to be added
    use_galaxy_templates: set to True to use galaxy templates as well
    """
    if use_galaxy_templates:
        copy_original_galaxy = original_galaxy
        no_of_templates_added = len(templates)
        if no_of_templates_added == 0:
            open(agn_param, 'w').write(copy_original_galaxy)
            print('No AGN templates added, just using EAZY galaxy templates')
            return
        for i in range(no_of_templates_added):
            current_id = last_id + 1 + i
            copy_original_galaxy = copy_original_galaxy + f'\n{current_id}   {agn_dir}{agn_temp_all[templates[i]]}   1.0 0 1.0    '
        open(agn_param, 'w').write(copy_original_galaxy)
        print(f'AGN templates added to the parameter file, {templates}, {agn_param}, {last_id} galaxy templates used')
        return
    else:
        copy_original_galaxy = original_empty
        no_of_templates_added = len(templates)
        if no_of_templates_added == 0:
            open(agn_param, 'w').write(copy_original_galaxy)
            print('No AGN templates added, no Galaxy templates used')
            return
        for i in range(no_of_templates_added):
            current_id = 0 + i
            copy_original_galaxy = copy_original_galaxy + f'\n{current_id}   {agn_dir}{agn_temp_all[templates[i]]}   1.0 0 1.0    '
            open(agn_param, 'w').write(copy_original_galaxy)
            print(f'AGN templates added to the parameter file {agn_param}, no galaxy templates used')
            return


# ----------------------------------------------------------------------------------------------------------------------
# AGN Info
agn_per_dir = f'inputs/{field}_agn_frac.txt'  # file with AGN fractions for each object, prepared in catalogue_prepare.ipynb
all_bayes = pd.read_csv(agn_per_dir, sep="\s+", comment='#')
main = pd.merge(main, all_bayes, on='id', how='left')  # AGN fraction for each object

has_fraction = 'bayes.agn.fracAGN' in main.columns
if not has_fraction:
    main = pd.merge(main, all_bayes, on='id', how='left') # AGN fraction for each object

mean_frac = np.mean(main['bayes.agn.fracAGN'])
positive_agn = main[main['bayes.agn.fracAGN'] > 0]
# ----------------------------------------------------------------------------------------------------------------------

# EAZY parameters

# following files should be in the same dir as the project
param_file = 'base.param'  # base parameter file, does not include all information
translate_file = glob.glob(f'zfourge/{field}/eazy/{field}.*.translate')

params = {}  # setting field specific parameters
params['Z_STEP'] = 0.05  # redshift step, defines the precision of each fit, 0.005 default

# inputs
params['TEMPLATES_FILE'] = agn_param  # parameter file containing which templates will be used
params['CACHE_FILE'] = f'{output_location}/{field}/{test_title}/tempfilt_{field}_{id_key}_{params["Z_STEP"]}.dat' # template cache file, not used
params['CATALOG_FILE'] = f'inputs/eazy_auto.cat'  # for cut catalogue created in the earlier cell


# ----------------------------------------------------------------------------------------------------------------------
def eazy_single_template(template):
    """
    Function to run EAZY for a single template, given by the number identifier
    """

    # Setup for output
    agn_sed = [template]  # AGN templates to be added, comma separated list
    output_directory = f'{output_location}/{field}/{test_title}/{field}_{test_title}_{id_key}_{agn_sed}_{use_galaxy_templates}'  # output directory for images
    params['MAIN_OUTPUT_FILE'] = output_directory

    agn_template_loader([template], use_galaxy_templates)

    # ------------------------------------------------------------------------------------------------------------------

    self = eazy.photoz.PhotoZ(param_file=param_file, translate_file=translate_file[0], zeropoint_file=None,
                              params=params, load_prior=True, load_products=False)

    # ------------------------------------------------------------------------------------------------------------------
    # Iterative Zero-point corrections

    NITER = 3  # no. of iterations
    NBIN = np.minimum(self.NOBJ // 100, 180)  # no. of bins

    for iter in range(NITER):
        sn = self.fnu / self.efnu
        clip = (sn > 1).sum(axis=1) > 4  # Generally make this higher to ensure reasonable fits
        self.iterate_zp_templates(idx=self.idx[clip], update_templates=False,
                                  update_zeropoints=True, iter=iter, n_proc=8,
                                  save_templates=False, error_residuals=False,
                                  NBIN=NBIN, get_spatial_offset=False)
    # ------------------------------------------------------------------------------------------------------------------
    # Turn off error corrections derived above
    self.set_sys_err(positive=True)

    # Full catalog
    sample = np.isfinite(self.ZSPEC)

    # fit_parallel renamed to fit_catalog 14 May 2021
    self.fit_catalog(self.idx[sample], n_proc=8)

    # Show zspec-zphot comparison
    zmax = 6
    fig = self.zphot_zspec(zmax=zmax)
    fig.savefig(
        f'{output_location}/{field}/{test_title}/zphot_zspec_{field}_{id_key}_{params["Z_STEP"]}_{agn_sed}_{use_galaxy_templates}.png')

    # ------------------------------------------------------------------------------------------------------------------
    main['ZSPEC'] = self.ZSPEC
    main['ZPHOT'] = self.zbest  # adding zspec and zphot to the agn_frac df
    main['chi2'] = self.chi2_best / self.NFILT  # adding chi2/N to the main df

    # sort main df by zspec
    main_red = main[main['ZSPEC'] > 0.005]  # filter
    main_red = main_red[main_red['ZPHOT'] > 0.02]  # filter
    main_red = main_red[main_red['ZSPEC'] <= zmax]  # filter
    main_red = main_red.sort_values(by='ZSPEC')  # sort

    # total NMAD
    zspec_nmad = np.array(main_red['ZSPEC'])
    zphot_nmad = np.array(main_red['ZPHOT'])
    dz = (zphot_nmad - zspec_nmad) / (1 + zspec_nmad)
    total_nmad = astropy.stats.mad_std(dz)

    # outliers
    spec_count = len(dz)
    catastrophic_limit = 0.15  # catastrophic_limit in eazy code
    outliers = np.abs(dz) >= catastrophic_limit
    outliers_count = sum(outliers)

    # ------------------------------------------------------------------------------------------------------------------

    key_data = pd.DataFrame(columns=headings)
    key_data.loc[0] = [id_key, params['Z_STEP'], loop_number, agn_sed, use_galaxy_templates, total_count, mean_frac,
                       spec_count, outliers_count, total_nmad]
    key_data.to_csv(key_data_file, mode='a', index=False, header=False)

    induvidual_data = pd.DataFrame(columns=['id', 'phot_redshift', 'chi2'])
    induvidual_data['id'] = self.idx
    induvidual_data['phot_redshift'] = self.zbest
    induvidual_data['chi2'] = self.chi2_best
    for i in range(self.fmodel.shape[1]):
        induvidual_data[f'band_{i}'] = self.fmodel[:, i]
    induvidual_data.to_csv(
        f'{output_location}/{field}/{test_title}/induvidual_data_{field}_{id_key}_{params["Z_STEP"]}_{agn_sed}_{use_galaxy_templates}.csv',
        index=False)
    # ------------------------------------------------------------------------------------------------------------------


# Looping
no_of_templates = len(agn_temp_all)
all_time = []

if __name__ == '__main__':

    for j in range(no_of_templates):

        pool = mp.Pool(processes=4, maxtasksperchild=500)  # EAZY runs mutliprocessing, starting the pool here to avoid memory issues

        i = j

        if i >= no_of_templates:
            print('All templates done')
            break

        all_time.append(time.ctime())
        eazy_single_template(i)
        pool.close()

        print(f'template {i} done')
        print(all_time)
        print('--------------------------------------------------------------')