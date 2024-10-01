# will run each template individually
# each test needs a new test name, otherwise ti will append to the old name directory and main file.
# change the key, field, zstep for each test
# for each EAZY run, saves the self object to a pickle file, which can be loaded and used for further analysis
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
import global_settings as gs
import pickle

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
# The naming convention suggests that 'individual' is at the start of the tests from this file
template_key_dict = gs.get_template_dict() # loading is the same for any test, so it can be done here

# ----------------------------------------------------------------------------------------------------------------------
# Directories
output_location = 'G:/honours/outputs'

# ----------------------------------------------------------------------------------------------------------------------
# EAZY parameters

# following files should be in the same dir as the project
param_file = 'base.param'  # base parameter file, does not include all information
agn_param = 'templates/eazy_v1.3_AGN_auto.param'  # parameter file with agn templates

# ----------------------------------------------------------------------------------------------------------------------


def eazy_single_template(test_title, template, field, id_key, template_key, use_galaxy_templates, use_prior, template_combos, z_step):
    """
    Function to run EAZY for a single template, given by the number identifier
    """
    if not os.path.isdir(f'{output_location}/{field}/{test_title}'):
        os.makedirs(f'{output_location}/{field}/{test_title}')

    # ------------------------------------------------------------------------------------------------------------------
    field = field
    id_key = id_key
    template_key = template_key
    id_key_dict = gs.get_id_dict(field)
    template_key = template_key

    # AGN templates allocation
    use_galaxy_templates = use_galaxy_templates  # set to True to use galaxy templates as well
    # Load any templates from the AGN template library

    agn_dir = template_key_dict[template_key]  # dir with all agn templates
    agn_temp_all = os.listdir(agn_dir)

    params = {}
    params['Z_STEP'] = z_step
    params['CATALOG_FILE'] = 'inputs/eazy_auto.cat'
    params['TEMPLATE_COMBOS'] = template_combos
    params['APPLY_PRIOR'] = use_prior
    params['TEMPLATES_FILE'] = agn_param
    params['CACHE_FILE'] = f'{output_location}/{field}/{test_title}/tempfilt_{field}_{id_key}_{z_step}.dat'
    translate_file = glob.glob(f'zfourge/{field}/eazy/{field}.*.translate')
    # ------------------------------------------------------------------------------------------------------------------
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
    # ------------------------------------------------------------------------------------------------------------------
    # AGN Info
    agn_per_dir = f'inputs/{field}_agn_frac.txt'  # file with AGN fractions for each object, prepared in catalogue_prepare.ipynb
    all_bayes = pd.read_csv(agn_per_dir, sep="\s+", comment='#')
    main = pd.merge(main, all_bayes, on='id', how='left')  # AGN fraction for each object

    has_fraction = 'bayes.agn.fracAGN' in main.columns
    if not has_fraction:
        main = pd.merge(main, all_bayes, on='id', how='left')  # AGN fraction for each object

    mean_frac = np.mean(main['bayes.agn.fracAGN'])
    # ------------------------------------------------------------------------------------------------------------------
    # Setup for output
    agn_sed = template  # AGN templates to be added, comma separated list
    templates_use = gs.check_template(agn_sed, agn_temp_all)  # check if all templates are wanted
    output_directory = gs.save_directory(output_location, field, test_title, id_key, template_key, agn_sed, use_galaxy_templates, z_step, template_combos) # output directory
    params['MAIN_OUTPUT_FILE'] = output_directory

    gs.agn_template_loader(templates_use, agn_param=agn_param, agn_dir=agn_dir, agn_temp_all=agn_temp_all, use_galaxy_templates=use_galaxy_templates)

    # ------------------------------------------------------------------------------------------------------------------
    print(params)
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
    self.crps_val = gs.crps(self.lnp, self.zgrid, self.ZSPEC)
    mean_CRPS = np.mean(self.crps_val[self.crps_val >= 0])
    total_nmad, outlier_nmad, outlier_count, outlier_fraction = gs.nmad_calc(main_red['ZPHOT'], main_red['ZSPEC'], outlier=True)
    outlier_scatter = gs.rms_calc(main_red['ZPHOT'], main_red['ZSPEC'], outlier=True)
    bias = np.median(main_red['ZPHOT'] - main_red['ZSPEC'])
    spec_count = len(main_red['ZSPEC'])

    # ------------------------------------------------------------------------------------------------------------------

    gs.save_key_data(output_location=output_location, field=field, test_title=test_title,
                     key_input=[id_key, z_step, template_key, templates_use, use_galaxy_templates,
                                template_combos, total_count, mean_frac, spec_count, outlier_count, total_nmad,
                                outlier_nmad, outlier_scatter, outlier_fraction, bias, mean_CRPS])

    with open(f'{output_directory}_individual_data.pkl', 'wb') as file:
        pickle.dump(self, file)
    # ------------------------------------------------------------------------------------------------------------------


# Looping
loop_style = 'nothing'  # loop style, just saved loadouts for loops
"""
main loops through all fields, id_keys, and template sets that are set
single_loop loops through all templates in the agn template directory for a single field and id_key
recommendation is to use to test a recommendation
"""
#

if __name__ == '__main__':
    if loop_style == 'main':

        print('Running main loop')

        test_title = 'main_template_sets_3'  # title of the test, eg. 1,2, A, B, Initial.


        all_time = []
        all_fields = ['cdfs2', 'cosmos2', 'uds']
        all_id_keys = ['normal', 'lacy_no_ovlp', 'xray_agn_no_ovlp', 'donley']
        all_template_keys = ['recommendation', 'EAZY', 'atlas_all', 'XMM']
        templates = {'atlas_rest': ['all'], 'atlas_all': ['all'], 'XMM': ['all']} # what templates for each key

        recommendation_df = pd.read_csv(f'{output_location}/other_data/recommendation_best-no-of-templates.csv')  # load the recommendation file

        for field in all_fields:
            for id_key in all_id_keys:
                for template_key in all_template_keys:
                    print('--------------------------------------------------------------')
                    print(f'Field: {field}, id_key: {id_key}, template_key: {template_key}')
                    print('--------------------------------------------------------------')

                    pool = mp.Pool(processes=4,
                                   maxtasksperchild=500)  # EAZY runs mutliprocessing, starting the pool here to avoid memory issues
                    all_time.append(time.ctime())

                    if template_key == 'EAZY': # Want to run only the eazy templates
                        agn_sed = []
                        template_set = 'atlas_rest' # doesn't matter here
                        use_galaxy_templates = True
                        z_step = 0.05
                        t_combos = 'a'
                        use_prior = 'y'
                        eazy_single_template(test_title, agn_sed, field, id_key, template_set, use_galaxy_templates, use_prior, t_combos, z_step)

                    elif template_key == 'recommendation': # Want to run the recommendation templates

                        if field == 'uds' and id_key == 'x_ray_agn_no_ovlp': # no recommendations for this field and id_key
                            continue

                        template_list = recommendation_df[(recommendation_df['field'] == field)
                                                          & (recommendation_df['id_key'] == id_key)]['templates'].iloc[0]
                        agn_sed = gs.stringlist_to_list(template_list)
                        print(f'Using the templates {agn_sed}')
                        template_set = 'atlas_rest'
                        use_galaxy_templates = True
                        z_step = 0.05
                        t_combos = 'a'
                        use_prior = 'y'
                        eazy_single_template(test_title, agn_sed, field, id_key, template_set, use_galaxy_templates, use_prior,
                                             t_combos, z_step)

                    else:
                        for template1 in templates[template_key]: # used for single template fits
                            template2 = template1 # if the loop crashes, change this to start at any template
                            use_galaxy_templates = False
                            z_step = 0.05
                            t_combos = 1
                            use_prior = 'n'
                            eazy_single_template(test_title, template2, field, id_key, template_key, use_galaxy_templates, use_prior,
                                                 t_combos, z_step)

                    pool.close()
                    print(all_time)
                    print('--------------------------------------------------------------')

    elif loop_style == 'single_loop':

        print('Running single loop')
        template_set = 'atlas_rest'
        use_galaxy_templates = True
        use_prior = 'y'
        t_combos = 'a'
        z_step = 0.05
        agn_dir = template_key_dict[template_set]  # dir with all agn templates
        agn_temp_all = os.listdir(agn_dir)
        field = 'uds'
        id_key = 'donley'

        test_title = f'individual_{field}_{id_key}_{z_step}_{use_galaxy_templates}'  # title of the test, eg. 1,2, A, B, Initial.

        no_of_templates = len(agn_temp_all)
        all_time = []

        for j in range(no_of_templates):

            pool = mp.Pool(processes=4, maxtasksperchild=500)  # EAZY runs multiprocessing, starting the pool here to avoid memory issues

            i = j
            if i >= no_of_templates:
                print('All templates done')
                break
            all_time.append(time.ctime())
            eazy_single_template(test_title, [i], field, id_key, template_set, use_galaxy_templates, use_prior, t_combos,
                                 z_step)
            pool.close()

            print(f'template {i} done')
            print(all_time)
            print('--------------------------------------------------------------')

    elif loop_style == 'recommendation':
        # loops through all the recommended templates through the given field and id_key

        print('Running recommendation loop')
        template_set = 'atlas_rest'
        agn_dir = template_key_dict[template_set]
        use_galaxy_templates = True
        use_prior = 'y'
        t_combos = 'a'
        z_step = 0.05
        field = 'cosmos2'
        id_key = 'normal'
        recommendation_list = [27, 9, 29, 22, 37, 31, 30, 25, 35, 24] # insert the recommended templates here

        test_title = f'recommendation_added_{field}_{id_key}'  # title of the test, eg. 1,2, A, B, Initial, change the second chunk to the method used to get the recommendations.

        for idi, template in enumerate(recommendation_list):

            templates_to_run = recommendation_list[:idi+1]
            print(templates_to_run)
            eazy_single_template(test_title, templates_to_run, field, id_key, template_set, use_galaxy_templates, use_prior, t_combos, z_step)

    print('######################################################################################')
    print('All done')