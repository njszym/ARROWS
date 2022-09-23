from pymatgen.core.composition import Composition
import numpy as np


def get_products(precursors, temperature, exp_data):
    """
    Check if exp_data contains results from the specified
    precursors and synthesis temperature. If so, return the
    corresponding products and their weight fractions.

    Args:
        precursors (list): list of starting materials.
        temperature (int): synthesis temperature (deg. C).
        exp_data (dict): dictionary containing all available
            experimental synthesis outcomes.
    Returns:
        products (list): chemical formulae of the
            reaction products.
        weight_fracs (list): weight fractions of
            the reaction products.
    """

    # Format exp json keys (reduced formulae, ordered alphabetically)
    initial_keys = list(exp_data.keys()).copy()
    initial_length = len(initial_keys)
    for orig_key in initial_keys:
        if ', ' in orig_key:
            phasenames = orig_key.split(', ')
            formulae = [Composition(cmpd).reduced_formula for cmpd in phasenames]
            ordered_formulae = sorted(formulae)
            new_key = ', '.join(ordered_formulae)
            if new_key != orig_key:
                exp_data[new_key] = exp_data[orig_key]
                del exp_data[orig_key]
        elif ',' in orig_key:
            phasenames = orig_key.split(',')
            formulae = [Composition(cmpd).reduced_formula for cmpd in phasenames]
            ordered_formulae = sorted(formulae)
            new_key = ', '.join(ordered_formulae)
            if new_key != orig_key:
                exp_data[new_key] = exp_data[orig_key]
                del exp_data[orig_key]

    # Ensure no keys were lost
    assert len(exp_data.keys()) == initial_length, 'Something went wrong with Exp.json formatting'

    # Modifications must not be made to parent var
    solid_precursors = precursors.copy()

    # Gaseous phases need not be specified in json
    if 'O2' in solid_precursors:
        solid_precursors.remove('O2')
    if 'CO2' in solid_precursors:
        solid_precursors.remove('CO2')
    if 'H2O' in solid_precursors:
        solid_precursors.remove('NH3')
    if 'NH3' in solid_precursors:
        solid_precursors.remove('H2O')

    # Precursors must be formatted as reduced formulae
    solid_precursors = [Composition(cmpd).reduced_formula for cmpd in solid_precursors]
    solid_precursors = sorted(solid_precursors)

    # Formatting consistent with exp dictionary keys
    precursor_key = ', '.join(solid_precursors)
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
    assert len(products) == len(weight_fracs), 'Number of products (%s) is not equal to number of weight fractions (%s)' % (len(products), len(weight_fracs))
    weight_fracs = [val/100 for val in weight_fracs]

    return products, weight_fracs

def get_xrd(precursors, temperature, exp_data):
    """
    Retrieve X-ray diffraction data.

    Args:
        precursors (list): list of starting materials.
        temperature (int): synthesis temperature (deg. C).
    Returns:
        x (list): two-theta (degrees).
        y (list): intensity (raw, unscaled).
    """

    # Precursors must be formatted as reduced formulae
    precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    precursors = sorted(precursors)

    # Formatting consistent with exp dictionary keys
    precursor_key = ', '.join(precursors)
    temp_key = '%s C' % int(temperature)

    # Make sure experimental data is available
    assert exp_data[precursor_key]['Temperatures'][temp_key]['Experimentally Verified'] is True, 'No XRD available'

    # xy data for XRD pattern
    x = exp_data[precursor_key]['Temperatures'][temp_key]['XRD']['x']
    y = exp_data[precursor_key]['Temperatures'][temp_key]['XRD']['y']

    return x, y
