# functions with values used in multiple files

import os


def is_working():
    print('Working')


def get_id_dict(field):
    """
    Return a dictionary with the keys for the different catalogues
    """
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
    return id_key_dict


def get_template_dict():
    """
    Return a dictionary with the keys for the different template sets
    """
    template_key_dict = {
        'atlas_rest': 'templates/hlsp_agnsedatlas_rest/',
        'atlas_observed': 'templates/hlsp_agnsedatlas_observed/',
        'atlas_composite': 'templates/hlsp_agnsedatlas_rest_composite/',
        'atlas_2014': 'templates/hlsp_agnsedatlas_2014/',
        'atlas_all': 'templates/hlsp_agnsedatlas_all/',
        'XMM': 'templates/MARA23032010/'
    }
    return template_key_dict


def nmad_calc(phot, spec, outlier=False):
    """
    Calculate the NMAD of the difference between photometric and spectroscopic redshifts
    """
    import numpy as np
    catastrophic_limit = 0.15 # catastrophic_limit in eazy code
    if len(phot) != len(spec):
        print('Lengths of the photometric and spectroscopic redshifts do not match')
        print(f'Phot: {len(phot)}, Spec: {len(spec)}')
        return
    phot = np.array(phot)
    spec = np.array(spec)
    diff = (phot - spec) / (1 + spec)
    if outlier:
        diff_noout = diff[abs(diff) < catastrophic_limit]
        outlier_count = len(diff) - len(diff_noout)
        if len(diff_noout) == 0:
            # if there are no non-outliers
            return 1.4826 * np.median(abs(diff - np.median(diff))), np.nan, outlier_count, 1
        outlier_fraction = 1 - (len(diff_noout) / len(phot))
        return 1.4826 * np.median(abs(diff - np.median(diff))), 1.4826 * np.median(abs(diff_noout - np.median(diff_noout))), outlier_count, outlier_fraction
    else:
        return 1.4826 * np.median(abs(diff - np.median(diff)))


def rms_calc(phot, spec, outlier=False):
    """
    Calculate the RMS of the difference between photometric and spectroscopic redshifts
    """
    import numpy as np
    if len(phot) != len(spec):
        print('Lengths of the photometric and spectroscopic redshifts do not match')
        print(f'Phot: {len(phot)}, Spec: {len(spec)}')
        return
    phot = np.array(phot)
    spec = np.array(spec)
    diff = (phot - spec) / (1 + spec)
    if outlier:
        diff = diff[abs(diff) < 0.15]
        return np.sqrt(np.mean(diff**2))
    else:
        return np.sqrt(np.mean(diff**2))


def agn_template_loader(templates, agn_param, agn_dir, agn_temp_all, use_galaxy_templates=True):

    """
    Function to load AGN templates to the parameter file
    templates: list of templates to be added
    use_galaxy_templates: set to True to use galaxy templates as well
    """
    temp_param = 'templates/eazy_v1.3.spectra.param'  # basic parameter file, no agn templates
    last_id = 9  # id of the last template in the parameter file
    empty_param = 'templates/eazy_v1.3_empty.param'  # empty parameter file

    # opening the parameter files, and reading the contents
    with open(temp_param) as f:
        original_galaxy = f.read()

    with open(empty_param) as f:
        original_empty = f.read()

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
            current_id = 1 + i
            copy_original_galaxy = copy_original_galaxy + f'\n{current_id}   {agn_dir}{agn_temp_all[templates[i]]}   1.0 0 1.0    '
        open(agn_param, 'w').write(copy_original_galaxy)
        print(f'AGN templates added to the parameter file {agn_param}, no galaxy templates used')
        return


def check_template(agn_sed, agn_temp_all):
    """
    Check if all templates are wanted
    """
    import numpy as np
    if 'all' in agn_sed:
        return np.linspace(0, len(agn_temp_all) - 1, len(agn_temp_all), dtype=int).tolist()
    else:
        return agn_sed


def inv_check_template(agn_sed, template_key):
    """
    Check if the agn_sed is 'all' agn templates in directory through a length match. Try to avoid using 2 of the same template, as this will cause an unmentioned error
    """
    import glob

    if len(agn_sed) == 2: # if just square brackets
        return agn_sed

    agn_sed = [int(item) for item in agn_sed[1:-1].split(',')]  # convert string to list
    agn_temp_all = glob.glob(get_template_dict()[template_key] + '*')
    if len(list(agn_sed)) == len(agn_temp_all):
        return 'all'
    else:
        return agn_sed


def save_directory(output_location, field, test_title, id_key, template_key, agn_sed, use_galaxy_templates, z_step, t_combos, *args):
    """
    Function to produce a directory for saving the results
    """
    if not os.path.exists(output_location):
        os.makedirs(output_location)

    if not os.path.exists(f"{output_location}/{field}"):
        os.makedirs(f"{output_location}/{field}")

    if not os.path.exists(f"{output_location}/{field}/{test_title}"):
        os.makedirs(f"{output_location}/{field}/{test_title}")
    if args:
        return f"{output_location}/{field}/{test_title}/{field}_{test_title}_{id_key}_{template_key}_{agn_sed}_{use_galaxy_templates}_{z_step}_{t_combos}" + '_'.join([str(i) for i in args])
    else:
        return f"{output_location}/{field}/{test_title}/{field}_{test_title}_{id_key}_{template_key}_{agn_sed}_{use_galaxy_templates}_{z_step}_{t_combos}"


def save_key_data(output_location, field, test_title, key_input):
    """
    Function to produce a key for saving the results. Key data as a list of the same length and order given in headings
    """
    import pandas as pd
    key_data_file = f'{output_location}/{field}/{test_title}/{test_title}_data.csv'
    headings = ['id_key', 'zstep', 'template_key', 'agn_templates', 'galaxy templates', 't_combos', 'total_obj', 'mean_agn_frac',
                'spec_count', 'outlier_count', 'nmad_val', 'outlier_nmad', 'outlier_scatter', 'outlier_fraction',
                'bias', 'mean CRPS']

    key_data = pd.DataFrame(columns=headings)
    if not os.path.isfile(key_data_file):
        key_data.to_csv(key_data_file, index=False)


    key_data = pd.DataFrame([key_input], columns=headings)
    key_data.to_csv(key_data_file, mode='a', index=False, header=False)


def load_individual(output_location, field, test_title, id_key, template_key, agn_sed, use_galaxy_templates, z_step, t_combos, *args):
    """
    Function to produce a DataFrame from the directory, made to ensure that directories with different naming schemes are picked up
    """
    import pandas as pd
    import pickle
    output_directory = save_directory(output_location, field, test_title, id_key, template_key, agn_sed, use_galaxy_templates, z_step, t_combos, *args)
    try:
        with open(f'{output_directory}_individual_data.pkl', 'rb') as file:
            self = pickle.load(file)

        individual_data = pd.DataFrame(columns=['id', 'phot_redshift', 'chi2'])
        individual_data['id'] = self.idx
        individual_data['phot_redshift'] = self.zbest
        individual_data['chi2'] = self.chi2_best
        for i in range(self.fmodel.shape[1]):
            individual_data[f'band_{i}'] = self.fmodel[:,i]
        return individual_data
    except (FileNotFoundError, pickle.UnpicklingError, EOFError, AttributeError):
        try:
            return pd.read_csv(f'{output_directory}_individual_data.csv')
        except (FileNotFoundError, PermissionError):
            return pd.read_csv(
                f'{output_location}/{field}/{test_title}/induvidual_data_{field}_{id_key}_{z_step}_{agn_sed}_{use_galaxy_templates}.csv')


def load_key_data(output_location, field, test_title):
    """
    Function to load the key data
    """
    import pandas as pd
    return pd.read_csv(f'{output_location}/{field}/{test_title}/{test_title}_data.csv')


def load_self(main, output_location, field, test_title, id_key, template_key, agn_sed, use_galaxy_templates, z_step, t_combos, *args):
    """
    Function to load the self object
    """
    import pickle

    output_directory = save_directory(output_location, field, test_title, id_key, template_key, agn_sed, use_galaxy_templates, z_step, t_combos, *args)

    try:
        with open(f'{output_directory}_individual_data.pkl', 'rb') as file:
            return pickle.load(file)
    except (FileNotFoundError, pickle.UnpicklingError, EOFError, AttributeError):
        # need to load the individual data and then create the self object
        # Loading values that are present in catalog, done at the end of EAZY_test.py
        import pandas as pd
        import numpy as np
        columns_to_drop = ['id', 'phot_redshift', 'chi2']
        spec_data = pd.read_csv(f'{output_location}/{field}/spec_data.csv')

        reindex_list = list(main['id'])
        reindex_val = [x - 1 for x in reindex_list]  # cat is indexed from 1, not 0
        reindex_sel = spec_data['id'].isin(list(reindex_val))
        spec_data = spec_data[reindex_sel].reset_index(drop=True)

        filter_data = pd.read_csv(f'{output_location}/{field}/filter_data.csv')
        flux_data = pd.read_csv(f'{output_location}/{field}/flux_data.csv')

        filter_error = flux_data.filter(regex='^e_')
        filter_error = filter_error[reindex_sel].reset_index(drop=True)
        filter_flux = flux_data.filter(regex='^f_')
        filter_flux = filter_flux[reindex_sel].reset_index(drop=True)

        individual_df = load_individual(output_location, field, test_title, id_key, template_key, agn_sed, use_galaxy_templates, z_step, t_combos, *args)

        class self_load:
            def __init__(self, idx, ZSPEC, ZBEST, fnu, efnu, efmodel, pivot, zbest, chi2_best, NFILT):
                self.idx = np.array(main['id'].index)
                self.ZSPEC = np.array(spec_data['zspec'])
                self.fnu = np.array(filter_flux)
                self.efnu = np.array(filter_error)
                self.fmodel = np.array(individual_df.drop(columns=columns_to_drop))
                self.pivot = np.array(filter_data['pivot'])
                self.zbest = np.array(individual_df['phot_redshift'])
                self.chi2_best = np.array(individual_df['chi2'])
                self.NFILT = len(filter_data)

        return self_load(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


def normalizer(data):
    """
    Normalizes the data to a 0-1 scale
    """
    return (data - min(data)) / (max(data) - min(data))


def normalize_by_sum(data):
    """
    Normalizes the data by the sum
    """
    return data / sum(data)


def heaviside(data, threshold):
    """
    Heaviside function, 0 if less than, 1 if more than
    """
    import numpy as np
    data_copy = np.copy(data)
    for i in range(len(data)):
        if data_copy[i] >= threshold:
            data_copy[i] = 1
        else:
            data_copy[i] = 0
    return data_copy




def crps(lnp, zgrid_prob, zspec_prob, zmin=0.005):
    """
    Compute the Continuous Ranked Probability Score (CRPS), which is an assessment on the quality of the PDF
    """
    import numpy as np
    unlog_prob = np.exp(lnp)
    nobj = len(zspec_prob)
    crps_values = np.zeros(nobj)
    for i in range(nobj):
        if zspec_prob[i] <= zmin:
            crps_values[i] = np.nan
            continue
        cdf = normalizer(np.cumsum(unlog_prob[i]))
        crps_values[i] = np.trapz((cdf - heaviside(zgrid_prob, zspec_prob[i]))**2, zgrid_prob)
    return crps_values


def flux_to_mag(flux):
    """
    Converts flux to magnitude
    """
    import numpy as np
    return -2.5 * np.log10(flux) + 25


def redshifter(bands, redshift):
    """
    Redshifts the bands
    """
    return bands / (1 + redshift)


def parabola_fit(zgrid, lnp):
    """
    Analytic parabola fit from get_maxlnp_redshift in eazy.py, adjusted slightly
    """
    import numpy as np

    izmax = np.argmax(lnp) # index of the maximum value

    if izmax == 0 or izmax == len(lnp) - 1: # if the maximum value is at the edge
        return -99, -99

    # gets the 3 points around the maximum value
    x = np.array(zgrid[izmax-1 : izmax+2])
    y = np.array(lnp[izmax-1 : izmax+2])

    dx = np.diff(x).T
    dx2 = np.diff(x ** 2).T
    dy = np.diff(y).T

    c2 = (dy[1] / dx[1] - dy[0] / dx[0]) / (dx2[1] / dx[1] - dx2[0] / dx[0])
    c1 = (dy[0] - c2 * dx2[0]) / dx[0]
    c0 = y.T[0] - c1 * x.T[0] - c2 * x.T[0] ** 2

    zbest = -c1 / 2 / c2
    lnpmax = c2 * zbest ** 2 + c1 * zbest + c0

    return zbest, lnpmax


def prob_adder(lnps):
    """
    Adds the probabilities of the lnps
    """
    import numpy as np

    all_prob = np.zeros_like(lnps[0])
    for lnp in lnps:
        all_prob += np.exp(lnp)

    all_prob = normalize_by_sum(all_prob)

    return all_prob


def return_prob(lnp):
    """
    Returns the probability of the lnp
    """
    import numpy as np
    return np.exp(lnp)