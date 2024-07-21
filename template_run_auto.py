# will run each template individually
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

np.seterr(all='ignore')
warnings.simplefilter('ignore', category=AstropyWarning)
print('EAZYCODE = '+ str(os.getenv('EAZYCODE')) + '\n')

print(time.ctime() + '\n')

print(sys.version + '\n')

for module in ['numpy', 'scipy', 'matplotlib','astropy','eazy']:
    mod = importlib.import_module(module)
    print('{0:>20} : {1}'.format(module, mod.__version__))
# ----------------------------------------------------------------------------------------------------------------------
os.getcwd()

# Load ZFOURGE catalogue from local drive
test_title = 'individual' # title of the test, eg. 1,2, A, B, Initial.
field = 'cdfs' #'cdfs', 'cosmos', or 'uds'

# Choose ID key for the catalogue
# A key designates what object you want included
id_key = 'normal'

# Directories for key, name keys anything, just to keep track of any complex object choices made in catalogue_prepare.ipynb
id_key_dict = {'normal': 'inputs/alternate_catalogues/cdfs.range.(0, -1).cat',
      'fraction_bin_10': f'inputs/alternate_catalogues/{field}.fraction.bin10.0.cat',
   'luminosity_bin_10' : f'inputs/alternate_catalogues/{field}.luminosity.bin10.0.cat'}

# AGN templates allocation
loop_number = 1 # what loop you are on

use_galaxy_templates = True # set to True to use galaxy templates as well

# ----------------------------------------------------------------------------------------------------------------------

# Directories
if not os.path.isdir(f'outputs/{field}/{test_title}'):
    os.makedirs(f'outputs/{field}/{test_title}')

# Where to Save number data

key_data_file = f'outputs/{field}/{test_title}/{test_title}_data.csv'
headings = ['id_key', 'zstep', 'loop_number', 'agn_templates', 'galaxy templates', 'total_obj', 'mean_agn_frac', 'spec_count', 'outlier_count', 'nmad_val']

key_data = pd.DataFrame(columns=headings)
if not os.path.isfile(key_data_file):
    key_data.to_csv(key_data_file, index=False)

# ----------------------------------------------------------------------------------------------------------------------
#Read the Catalogue

main_cat = pd.read_csv(id_key_dict[id_key]) # get the catalogue for the id_key
main_cat.to_csv('inputs/eazy_test.cat', index=False) # create a new catalogue, allows for change to be made in this cell

#Setting up the main catalogue
main = pd.read_csv('inputs/eazy_test.cat', sep=" ", comment="#", header=None, skipinitialspace=True) # opening cut cat, and adjusting it
headers = pd.read_csv('inputs/eazy_test.cat', sep=" ", header=None, nrows=1).iloc[0]
headers = headers[1:]
main.columns = headers

total_count = len(main) # all objects in the range

# ----------------------------------------------------------------------------------------------------------------------
# Load any templates from the AGN template library

temp_param = 'templates/eazy_v1.3.spectra.param' # basic parameter file, no agn templates
last_id = 9 # last id in the parameter file
empty_param = 'templates/eazy_v1.3_empty.param' # empty parameter file
agn_param = 'templates/eazy_v1.3_AGN.param' # parameter file with agn templates

# opening the parameter files, and reading the contents
with open(temp_param) as f:
    original_galaxy = f.read()

with open(empty_param) as f:
    original_empty = f.read()

agn_dir = 'templates/hlsp_agnsedatlas_observed/' # dir with all agn templates
agn_temp_all = os.listdir(agn_dir)
def agn_template_loader(templates, use_galaxy_templates=False):
    """
    Function to load AGN templates to the parameter file
    templates: list of templates to be added
    use_galaxy_templates: set to True to use galaxy templates as well
    """
    if use_galaxy_templates:
      copy = original_galaxy
      no_of_templates = len(templates)
      if no_of_templates == 0:
        open(agn_param, 'w').write(copy)
        print('No AGN templates added, just using EAZY galaxy templates')
        return
      for i in range(no_of_templates):
        id = last_id + i
        copy = copy + f'\n{id}   {agn_dir}{agn_temp_all[templates[i]]}   1.0 0 1.0    '
      open(agn_param, 'w').write(copy)
      print(f'AGN templates added to the parameter file, {templates}, {agn_param}, {last_id} galaxy templates used')
      return
    else:
      copy = original_empty
      no_of_templates = len(templates)
      if no_of_templates == 0:
        open(agn_param, 'w').write(copy)
        print('No AGN templates added, no Galaxy templates used')
        return
      for i in range(no_of_templates):
        id = 0 + i
        copy = copy + f'\n{id}   {agn_dir}{agn_temp_all[templates[i]]}   1.0 0 1.0    '
      open(agn_param, 'w').write(copy)
      print(f'AGN templates added to the parameter file {agn_param}, no galaxy templates used')
      return

# ----------------------------------------------------------------------------------------------------------------------
#AGN Info
agn_per_dir = f'inputs/{field}_agn_frac.txt' # file with AGN fractions for each object, prepared in catalogue_prepare.ipynb
all_bayes = pd.read_csv(agn_per_dir, sep="\s+", comment='#')
main = pd.merge(main, all_bayes, on='id', how='left') # AGN fraction for each object
mean_frac = np.mean(main['bayes.agn.fracAGN'])
positive_agn = main[main['bayes.agn.fracAGN'] > 0]
# ----------------------------------------------------------------------------------------------------------------------

# EAZY parameters

# following files should be in the same dir as the project
param_file = 'base.param' #base parameter file, does not include all information
translate_file = glob.glob(f'zfourge/{field}/eazy/{field}.*.translate')

params = {} # setting field specific parameters
params['Z_STEP'] = 0.005 # redshift step, defines the precision of each fit, 0.005 default

#inputs
params['TEMPLATES_FILE'] = 'templates/eazy_v1.3_AGN.param' # parameter file containing which templates will be used
params['CACHE_FILE'] = f'zfourge/{field}/{field}.tempfilt'
params['CATALOG_FILE'] = f'inputs/eazy_test.cat' # for cut catalogue created in the earlier cell



# ----------------------------------------------------------------------------------------------------------------------
def eazy_single_template(template):
    """
    Function to run EAZy for a single template, given by the number identifier
    """

    #Setup for output
    agn_sed = [template]  # AGN templates to be added, comma separated list
    output_directory = f'outputs/{field}/{test_title}/{field}_{test_title}_{id_key}_{agn_sed}_{use_galaxy_templates}'  # output directory for images
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
    fig = self.zphot_zspec(zmax=8)
    fig.savefig(
        f'outputs/{field}/{test_title}/zphot_zspec_{field}_{id_key}_{params["Z_STEP"]}_{agn_sed}_{use_galaxy_templates}.png')

    # ------------------------------------------------------------------------------------------------------------------
    # Derived parameters (z params, RF colors, masses, SFR, etc.)
    warnings.simplefilter('ignore', category=RuntimeWarning)
    zout, hdu = self.standard_output(simple=False,
                                    rf_pad_width=0.5, rf_max_err=2,
                                    prior=True, beta_prior=True,
                                    absmag_filters=[],
                                    extra_rf_filters=[])

    # 'zout' also saved to [MAIN_OUTPUT_FILE].zout.fits

    # ------------------------------------------------------------------------------------------------------------------
    main['ZSPEC'] = self.ZSPEC
    main['ZPHOT'] = self.zbest  # adding zspec and zphot to the agn_frac df
    main['chi2'] = self.chi2_best / self.NFILT  # adding chi2/N to the main df

    # sort main df by zspec
    main_red = main[main['ZSPEC'] > 0.005]  # filter
    main_red = main_red[main_red['ZPHOT'] > 0.02]  # filter
    main_red = main_red.sort_values(by='ZSPEC')  # sort

    # sort main df by agn fraction
    main_agn = main_red.dropna(subset=['bayes.agn.fracAGN'])  # filter
    main_agn = main_agn.sort_values(by='bayes.agn.fracAGN')  # sort

    # total NMAD
    dz = (np.array(main_red['ZPHOT']) - np.array(main_red['ZSPEC'])) / (1 + np.array(main_red['ZSPEC']))
    total_nmad = astropy.stats.mad_std(dz)

    # outliers
    spec_count = len(dz)
    catastrophic_limit = 0.15  # catastrophic_limit in eazy code
    outliers = np.abs(dz) >= catastrophic_limit
    outliers_count = sum(outliers)

    no_of_bins = 15  # no. of bins

    # ------------------------------------------------------------------------------------------------------------------
    logbins_start = np.log10(min(main_red['ZSPEC']))  # using logarithmic bins
    logbins_end = np.log10(max(main_red['ZSPEC'])) + 1e-10  # adding a small value to include the last value
    logbins = np.logspace(logbins_start, logbins_end, num=no_of_bins + 1)
    counts_red, bins_red = np.histogram(main_red['ZSPEC'], bins=logbins)
    print(f'Counts: {counts_red}')
    print(f'Total: {len(main_red["ZSPEC"])}')
    bin_centers_red = (np.array(bins_red[:-1]) + np.array(bins_red[1:])) / 2

    last_val = 0  # last val is the first value within a bin (say object 2334), while new val is the last
    nmad_red_val = []
    outlier_frac_red = []
    chi2_red = []
    for bin_loop in range(no_of_bins):
        new_val = last_val + counts_red[bin_loop] - 1
        zspec_val = main_red['ZSPEC'][last_val:new_val]
        zphot_val = main_red['ZPHOT'][last_val:new_val]
        dz = (np.array(zphot_val) - np.array(zspec_val)) / (1 + np.array(zspec_val))
        nmad_red_val.append(astropy.stats.mad_std(dz))  # NMAD
        outliers_bin = np.abs(dz) >= catastrophic_limit
        bin_fraction = sum(outliers_bin) / counts_red[bin_loop]
        outlier_frac_red.append(bin_fraction)
        chi2_med = np.median(main_red['chi2'][last_val:new_val])  # chi2
        chi2_red.append(chi2_med)
        last_val = new_val + 1  # adding 1 to skip the last value of the previous bin

    print(f'NMAD: {nmad_red_val}')
    print(f'Outliers: {outlier_frac_red}')
    print(f'Fits: {chi2_red}')

    plt.clf()
    fig, ax = plt.subplots(3, 1, sharex=True)

    # NMAD
    ax[0].plot(bin_centers_red, nmad_red_val, 'r--')
    ax[0].set_ylabel('NMAD')
    ax[0].text(0.25, 3.2, f'Total NMAD: {total_nmad:.4f}', fontsize=8, horizontalalignment='center',
               verticalalignment='center', transform=plt.gca().transAxes)

    # Outliers
    ax[1].plot(bin_centers_red, outlier_frac_red, 'b--')
    ax[1].set_ylabel('Outlier Fraction')

    # Chi2
    ax[2].plot(bin_centers_red, chi2_red, 'm--')
    ax[2].set_xlabel('Redshift')
    ax[2].set_ylabel('Chi2')

    fig.tight_layout()

    plt.savefig(
        f'outputs/{field}/{test_title}/RED_{field}_{id_key}_{params["Z_STEP"]}_{agn_sed}_{use_galaxy_templates}.png')

    # ------------------------------------------------------------------------------------------------------------------

    # NMAD against AGN fraction

    counts_agn, bins_agn = np.histogram(main_agn['bayes.agn.fracAGN'], bins=no_of_bins)
    print(f'Counts: {counts_agn}')
    bin_centers_agn = (np.array(bins_agn[:-1]) + np.array(bins_agn[1:])) / 2

    last_val = 0  # last val is the first value within a bin (say object 2334), while new val is the last
    nmad_agn_val = []
    outlier_frac_agn = []
    chi2_agn = []
    for bin_loop in range(no_of_bins):
        new_val = last_val + counts_agn[bin_loop] - 1
        zspec_val = main_agn['ZSPEC'][last_val:new_val]
        zphot_val = main_agn['ZPHOT'][last_val:new_val]
        dz = (np.array(zphot_val) - np.array(zspec_val)) / (1 + np.array(zspec_val))
        nmad_agn_val.append(astropy.stats.mad_std(dz))  # NMAD
        outliers_bin = np.abs(dz) >= 0.015  # outlier
        bin_fraction = sum(outliers_bin) / counts_agn[bin_loop]
        outlier_frac_agn.append(bin_fraction)
        chi2_med = np.median(main_agn['chi2'][last_val:new_val])  # chi2
        chi2_agn.append(chi2_med)
        last_val = new_val + 1  # adding 1 to skip the last value of the previous bin

    print(f'NMAD: {nmad_agn_val}')
    print(f'Outliers: {outlier_frac_agn}')
    print(f'Fits: {chi2_agn}')

    plt.clf()
    fig, ax = plt.subplots(3, 1, sharex=True)

    # NMAD
    ax[0].plot(bin_centers_agn, nmad_agn_val, 'r--')
    ax[0].set_ylabel('NMAD')
    ax[0].text(0.25, 3.2, f'Total NMAD: {total_nmad:.4f}', fontsize=8, horizontalalignment='center',
               verticalalignment='center', transform=plt.gca().transAxes)

    # Outliers
    ax[1].plot(bin_centers_agn, outlier_frac_agn, 'b--')
    ax[1].set_ylabel('Outlier Fraction')

    # Chi2
    ax[2].plot(bin_centers_agn, chi2_agn, 'm--')
    ax[2].set_xlabel('AGN Fraction')
    ax[2].set_ylabel('Chi2')

    fig.tight_layout()

    plt.savefig(
        f'outputs/{field}/{test_title}/AGN_{field}_{id_key}_{params["Z_STEP"]}_{agn_sed}_{use_galaxy_templates}.png')

    # ------------------------------------------------------------------------------------------------------------------

    # limit the no of objects to those that are present in all bands, and have a positive redshift
    flux_observed = []
    flux_residual_objects = []
    # any(self.fnu[i] == -99) or
    for i in range(len(self.fnu)):
        if self.zbest[i] < 0.001:
            continue
        else:
            flux_observed.append(self.fnu[i])
            flux_residual_objects.append(i)
    flux_observed = np.array(flux_observed)

    # find the residuals of the objects
    flux_model = self.fmodel
    flux_residual = np.zeros_like(flux_observed)
    residual_uncertainties = []
    for row in range(len(flux_residual_objects)):
        im = flux_residual_objects[row]
        residual_uncertainties.append(self.efnu[im] / flux_observed[row])
        for col in range(flux_observed.shape[1]):  # no. of bands
            if flux_observed[row, col] == -99:  # if the flux is -99, it is not recorded and should not be shown as such
                flux_residual[row, col] = -99
                continue
            else:
                flux_residual[row, col] = (flux_observed[row, col] - flux_model[im, col]) / flux_observed[row, col]
    residual_uncertainties = np.array(residual_uncertainties)

    # redshift the objects
    object_redshifted = []
    pivots = self.pivot
    for i in range(len(flux_residual_objects)):
        im = flux_residual_objects[i]
        redshifted = pivots / (1 + self.zbest[im])
        object_redshifted.append(redshifted)
    object_redshifted = np.array(object_redshifted)

    # ------------------------------------------------------------------------------------------------------------------

    # plot the residuals
    wavelength_flat = object_redshifted.flatten()
    residual_flat = flux_residual.flatten()
    uncertainties_flat = residual_uncertainties.flatten()

    # binning
    obj_per_bin = 2000  # EAZY used 2000
    no_of_bins_res = np.ceil(len(residual_flat) // obj_per_bin)
    res_sorted = pd.DataFrame(
        {'wavelength': wavelength_flat, 'residual': residual_flat, 'uncertainties': uncertainties_flat})
    res_sorted = res_sorted.sort_values(by='wavelength')
    res_sorted = res_sorted.loc[res_sorted['residual'] != -99]
    res_sorted = res_sorted.reset_index(drop=True)
    res_sorted_copy = res_sorted  # not abs
    res_sorted = abs(res_sorted)
    res_bin = pd.DataFrame(columns=["Median Residuals", "Median Wavelength", "Median Uncertainties", "Confidence"])
    res_bin_copy = pd.DataFrame(columns=["Median Residuals", "Median Wavelength", "Median Uncertainties", 'Confidence'])
    for i in range(int(no_of_bins_res)):
        min_loop = i * obj_per_bin
        max_loop = ((i + 1) * obj_per_bin) - 1
        med_res = np.median(res_sorted['residual'][min_loop:max_loop])
        med_res_copy = np.median(
            res_sorted_copy['residual'][min_loop:max_loop])  # only the residuals should be affected by abs
        med_wave = np.median(res_sorted['wavelength'][min_loop:max_loop])
        med_unc = np.median(res_sorted['uncertainties'][min_loop:max_loop])
        confidence = (50 / 68) * np.std(
            res_sorted['residual'][min_loop:max_loop])  # want 50% confidence interval, not 1sigma
        res_bin.loc[i] = [med_res, med_wave, med_unc, confidence]
        res_bin_copy.loc[i] = [med_res_copy, med_wave, med_unc, confidence]

    # ------------------------------------------------------------------------------------------------------------------

    plt.clf()
    fig, ax = plt.subplots(2, 1, sharex=True)
    fig.gca().set_xscale('log')
    fig.gca().xaxis.set_major_formatter('{x:.0f}')

    ax[0].plot(wavelength_flat, residual_flat, 'k,', alpha=0.1)
    ax[0].plot(res_bin_copy["Median Wavelength"], res_bin_copy["Median Residuals"], 'r--')
    ax[0].axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax[0].set_ylim(-0.6, 0.6)
    ax[0].set_xlim(1e3, 6e4)
    ax[0].set_ylabel('F-T/F')
    ax[0].set_title('Residuals of the objects at Restframe')

    ax[1].plot(res_bin["Median Wavelength"], res_bin["Median Residuals"], 'r--')
    ax[1].plot(res_bin["Median Wavelength"], res_bin["Median Uncertainties"], 'b-.')
    # ax[1].plot(res_bin["Median Wavelength"], res_bin["Confidence"], 'k-.')
    ax[1].set_ylim(0, 0.6)
    ax[1].set_ylabel('Median Residuals')
    ax[1].set_xlabel('Wavelength (Angstrom)')

    fig.tight_layout()

    plt.savefig(
        f'outputs/{field}/{test_title}/residuals_{field}_{id_key}_{params["Z_STEP"]}_{agn_sed}_{use_galaxy_templates}.png')

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
        f'outputs/{field}/{test_title}/induvidual_data_{field}_{id_key}_{params["Z_STEP"]}_{agn_sed}_{use_galaxy_templates}.csv',
        index=False)
    # ------------------------------------------------------------------------------------------------------------------

#Looping
no_of_templates = len(agn_temp_all)
all_time = []

for i in range(no_of_templates):
    if __name__ == '__main__':
        all_time.append(time.ctime())
        eazy_single_template(i)
        print(f'template {i} done')
        print(all_time)
        print('--------------------------------------------------------------')
