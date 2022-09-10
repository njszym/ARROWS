from pymatgen.core.composition import Composition
import pymatgen.analysis.phase_diagram as pd
from mp_api import MPRester
import numpy as np
import json


def get_hull_Ef(formula, cmpd_pd=None):
    """
    Retrieve formation energy of a given chemical formula
    according to its location on the convex hull.

    Note: at some point, I should include T-dependence
    """

    target_comp = Composition(formula)

    if cmpd_pd is None:
        all_entries = mpr.get_entries_in_chemsys([str(el) for el in target_comp.elements])
        cmpd_pd = pd.PhaseDiagram(all_entries)

    target_energy = cmpd_pd.get_hull_energy(target_comp)
    target_dict = target_comp.as_dict()
    competing_energies = 0.0
    sum_coeffs = 0.0

    for elem in target_dict.keys():
        elem_formula = '%s%s' % (elem, target_dict[elem])
        elem_comp = Composition(elem_formula)
        sum_coeffs += target_dict[elem]
        competing_energies += cmpd_pd.get_hull_energy(elem_comp)

    return (target_energy - competing_energies) / sum_coeffs

def get_entry_Ef(formula, temp, atmos='air', data_path='/Users/njszym/Research/MP/MP_stability.json'):
    """
    Retrieve entry formation energy from MP.
    Temperature-dependence is calculated via Bartel method.
    """

    # Partial pressures of gaseous species in air
    if atmos == 'air':
        p_O2 = 21200.
        p_CO2 = 4050.
        p_NH3 = 16.
        p_H2O = 2300.

    # Estimation based on 1e-6 partial pressure
    elif atmos == 'inert':
        p_O2 = 0.1
        p_CO2 = 0.1
        p_NH3 = 0.1
        p_H2O = 0.1

    else:
        raise Exception('Atmosphere must either be air or inert')

    target_comp = Composition(formula)
    ordered_formula = target_comp.alphabetical_formula
    condensed_formula = ordered_formula.replace(' ', '')

    with open(data_path) as fname:
        energy_data = json.load(fname)

    assert str(temp) in energy_data.keys(), """Invalid temperature.
    Only the following are allowed: %s""" % list(energy_data.keys())

    if condensed_formula in energy_data[str(temp)].keys():
        Ef = energy_data[str(temp)][condensed_formula]['Ef']
        elems = [str(el) for el in target_comp.elements]

        # Convert Ef to grand potential with open O2
        if 'O' in elems:
            n_O = target_comp.fractional_composition.as_dict()['O']
            dG = n_O*get_chempot_correction('O', temp, 21200)
            Ef -= dG

        # Rely on hull energy for gaseous species
        if Composition(formula).reduced_formula in ['CO2', 'O2', 'H2O', 'NH3', 'H3N']:
            return None

    else:
        Ef = None

    return Ef

def make_phase_diagram(entries, elems, temp, atmos='air', data_path='/Users/njszym/Research/MP/MP_stability.json'):

    # Partial pressures of gaseous species in air
    if atmos == 'air':
        p_O2 = 21200
        p_CO2 = 4050
        p_NH3 = 16.
        p_H2O = 2300.

    # Estimation based on 1e-6 partial pressure
    elif atmos == 'inert':
        p_O2 = 0.1
        p_CO2 = 0.1
        p_NH3 = 0.1
        p_H2O = 0.1

    else:
        raise Exception('Atmosphere must either be air or inert')

    zero_temp_pd = pd.PhaseDiagram(entries)

    if temp == 0.0:
        return zero_temp_pd

    with open(data_path) as fname:
        energy_data = json.load(fname)

    # Apply temperature corrections via Bartel method
    revised_entries = []
    for entry in entries:
        entry_dict = entry.as_dict()
        ordered_formula = entry.composition.alphabetical_formula
        condensed_formula = ordered_formula.replace(' ', '')
        if condensed_formula in energy_data[str(temp)].keys():
            E0 = energy_data['0'][condensed_formula]['Ef']
            Ef = energy_data[str(temp)][condensed_formula]['Ef']
            dE = entry.composition.num_atoms*(Ef - E0)
            entry_dict['energy'] += dE
        revised_entries.append(pd.PDEntry.from_dict(entry_dict))

    # Account for temperature-dependence of gaseous species
    if 'O' in elems:

        # O2
        O_comp = Composition('O')
        O_energ = zero_temp_pd.get_hull_energy(O_comp)
        O_energ += get_chempot_correction('O', temp, p_O2)
        O_entry = pd.PDEntry(O_comp, O_energ)
        revised_entries.append(O_entry)

        # CO2
        if 'C' in elems:
            CO2_comp = Composition('CO2')
            CO2_energ = zero_temp_pd.get_hull_energy(CO2_comp)
            CO2_energ += 3*get_chempot_correction('CO2', temp, p_CO2)
            CO2_entry = pd.PDEntry(CO2_comp, CO2_energ)
            revised_entries.append(CO2_entry)

        if 'H' in elems:
            H2O_comp = Composition('H2O')
            H2O_energ = zero_temp_pd.get_hull_energy(H2O_comp)
            H2O_energ += 3*get_chempot_correction('H2O', temp, p_H2O)
            H2O_entry = pd.PDEntry(H2O_comp, H2O_energ)
            revised_entries.append(H2O_entry)

    if ('N' in elems) and ('H' in elems):

        # NH3
        NH3_comp = Composition('NH3')
        NH3_energ = zero_temp_pd.get_hull_energy(NH3_comp)
        NH3_energ += 4*get_chempot_correction('NH3', temp, p_NH3)
        NH3_entry = pd.PDEntry(NH3_comp, NH3_energ)
        revised_entries.append(NH3_entry)

    return pd.PhaseDiagram(revised_entries)

def get_pd_dict(available_precursors, temperatures, atmos='air'):
    """
    Build a dictionary containing phase diagrams for the
    given chemical space throughout a range of temperatures.
    """

    mpr = MPRester('YdSP1Z8Tv8vfrLSaRCWdwp3EWvazS0vf')

    # Determine elements from precursors
    elems = []
    for cmpd in available_precursors:
        comp = Composition(cmpd)
        elems += [str(el) for el in comp.elements]
    elems = list(set(elems))

    # Get all entries in the chemical space
    entries = mpr.get_entries_in_chemsys(elems)

    # Build phase diagram dictionary
    pd_dict = dict()
    for temp in temperatures:
        pd_dict[temp] = make_phase_diagram(entries, elems, temp, atmos)

    return pd_dict


def make_custom_entry(struc, energy):

    return pd.PDEntry(struc.composition, energy)

def get_chempot_correction(element, temp, pres):
    """
    Get the normalized correction term Δμ for chemical potential of a gas
    phase consisting of element at given temperature and pressure,
    referenced to that in the standard state (T_std = 298.15 K,
    T_std = 1 bar). The gas phase is limited to be one of O2, N2, Cl2,
    F2, H2. Calculation formula can be found in the documentation of
    Materials Project website.

    Args:
        element (string): The string representing the element.
        temp (float): The temperature of the gas phase.
        pres (float): The pressure of the gas phase.

    Returns:
        The correction of chemical potential in eV/atom of the gas
        phase at given temperature and pressure.
    """

    EV_TO_KJ_PER_MOL = 96.4853

    if element not in ['O', 'N', 'Cl', 'F', 'H', 'CO2', 'NH3', 'H2O']:
        return 0
    std_temp = 298.15
    std_pres = 1E5
    ideal_gas_const = 8.3144598
    # Cp and S at standard state in J/(K.mol). Data from
    # https://janaf.nist.gov/tables/O-029.html
    # https://janaf.nist.gov/tables/N-023.html
    # https://janaf.nist.gov/tables/Cl-073.html
    # https://janaf.nist.gov/tables/F-054.html
    # https://janaf.nist.gov/tables/H-050.html
    Cp_dict = {'O': 29.376,
               'N': 29.124,
               'Cl': 33.949,
               'F': 31.302,
               'H': 28.836,
               'CO2': 37.129,
               'NH3': 35.640,
               'H2O': 33.22}

    S_dict = {'O': 205.147,
              'N': 191.609,
              'Cl': 223.079,
              'F': 202.789,
              'H': 130.680,
              'CO2': 213.79,
              'NH3': 192.80,
              'H2O': 194.10}
    Cp_std = Cp_dict[element]
    S_std = S_dict[element]
    PV_correction = ideal_gas_const * temp * np.log(pres / std_pres)
    TS_correction = - Cp_std * (temp * np.log(temp) - std_temp * np.log(std_temp)) \
        + Cp_std * (temp - std_temp) * (1 + np.log(std_temp)) \
        - S_std * (temp - std_temp)

    dG = PV_correction + TS_correction

    # Convert to eV/molecule unit.
    dG /= 1000 * EV_TO_KJ_PER_MOL

    # Normalize by number of atoms in the gas molecule
    if element == 'H2O':
        dG /= 3
    if element == 'CO2':
        dG /= 3
    if element == 'NH3':
        dG /= 4
    if element == 'O':
        dG /= 2

    return dG
