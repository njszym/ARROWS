from pymatgen.core.composition import Composition
import numpy as np


def get_products(precursors, temperature, exp_data):

    # Precursors must be formatted as reduced formulae
    precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    precursors = sorted(precursors)

    # Formatting consistent with exp dictionary keys
    precursor_key = ', '.join(precursors)
    temp_key = '%s C' % int(temperature)

    # If precursors not sampled yet, return None
    if precursor_key not in exp_data.keys():
        return None, None

    # If temperature not sampled yet, return None
    if temp_key not in exp_data[precursor_key]['Temperatures'].keys():
        return None, None

    # Experimentally observed products
    products = exp_data[precursor_key]['Temperatures'][temp_key]['products']

    # Re-format into reduced formulae; ignore space group for now
    # Products are given as formula_spacegroup-number
    products = [cmpd_sg.split('_')[0] for cmpd_sg in products]
    products = [Composition(cmpd).reduced_formula for cmpd in products]

    # Weight fractions must sum to 100 (or 1, in which case re-normalize)
    weight_fracs = exp_data[precursor_key]['Temperatures'][temp_key]['product weight fractions']
    if np.isclose(sum(weight_fracs), 1, atol=0.01):
        weight_fracs = 100*np.array(weight_fracs)
    assert np.isclose(sum(weight_fracs), 100, atol=5.0), 'Weight fractions sum to %s, but they should sum to 100' % sum(weight_fracs)
    weight_fracs = [val/100 for val in weight_fracs]

    return products, weight_fracs

def get_xrd(precursors, temperature, exp_data):

    # Precursors must be formatted as reduced formulae
    precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    precursors = sorted(precursors)

    # Formatting consistent with exp dictionary keys
    precursor_key = ', '.join(precursors)
    temp_key = '%s C' % int(temperature)

    # Make sure experimental data is available
    assert exp_data[precursor_key]['Temperatures'][temp_key]['Experimentally Verified'] is True, 'No XRD available'

    # XRD
    x = exp_data[precursor_key]['Temperatures'][temp_key]['XRD']['x']
    y = exp_data[precursor_key]['Temperatures'][temp_key]['XRD']['y']

    return x, y
