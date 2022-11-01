from pymatgen.core.composition import Composition
import pymatgen.analysis.phase_diagram as pd
from mp_api.client import MPRester
import numpy as np
import json


def get_hull_Ef(formula, cmpd_pd=None):
    """
    Retrieve formation energy of a given chemical formula
    according to its location on the convex hull.

    Args:
        formula (str): chemical formula.
        cmpd_pd (pymatgen phase diagram):
            phase diagram in the composition
            space containing formula.
    Returns:
        final_energ (float): formation energy
    """

    # Composition of interest
    target_comp = Composition(formula)

    # If no phase diagram is provided, create one based on the given composition
    if cmpd_pd is None:
        all_entries = mpr.get_entries_in_chemsys([str(el) for el in target_comp.elements])
        cmpd_pd = pd.PhaseDiagram(all_entries)

    # Get hull energy at specified composition
    target_energy = cmpd_pd.get_hull_energy(target_comp)

    # Get energies of elemental references
    sum_coeffs = 0.0
    competing_energies = 0.0
    target_dict = target_comp.as_dict()
    for elem in target_dict.keys():
        elem_formula = '%s%s' % (elem, target_dict[elem])
        elem_comp = Composition(elem_formula)
        sum_coeffs += target_dict[elem]
        competing_energies += cmpd_pd.get_hull_energy(elem_comp)

    # Calculate difference between hull energy and elemental references
    final_energ = (target_energy - competing_energies) / sum_coeffs

    return final_energ

def get_entry_Ef(formula, temp, atmos='air', data_path='arrows/energetics/MP_Energetics.json'):
    """
    Retrieve the formation energy of a given chemical formula.
    Zero-temperature energy calculated using DFT (from Materials Project).
    Temperature-dependence of the energy is estimated via the Bartel method.
    https://doi.org/10.1038/s41467-018-06682-4

    Args:
        formula (str): chemical formula.
        temp (int/float): temperature.
        atmos (str): air or inert.
        data_path (str): relative path to
            the file containing the 0 K
            energies from MP.
    Returns:
        Ef: formation energy.
    """

    # Round temperature to nearest 100 C
    T = round(temp, -2)

    # Partial pressures of gaseous species in air
    if atmos == 'air':
        p_O2 = 21200.
        p_N2 = 79033.
        p_CO2 = 4050.
        p_NH3 = 16.
        p_H2O = 2300.
        p_NO = 0.1

    # Estimation based on 1e-6 partial pressure
    elif atmos == 'inert':
        p_O2 = 0.1
        p_N2 = 0.1
        p_CO2 = 0.1
        p_NH3 = 0.1
        p_H2O = 0.1
        p_NO = 0.1

    else:
        raise Exception('Atmosphere must either be air or inert')

    # Composition of interest
    target_comp = Composition(formula)
    ordered_formula = target_comp.alphabetical_formula
    condensed_formula = ordered_formula.replace(' ', '')

    # Load energies from MP + Bartel
    with open(data_path) as fname:
        energy_data = json.load(fname)

    # Ensure energy data is available at specified temperature
    assert str(T) in energy_data.keys(), """Invalid temperature.
    Only the following are allowed: %s""" % list(energy_data.keys())

    if condensed_formula in energy_data[str(T)].keys():

        # Load formation energy
        Ef = energy_data[str(T)][condensed_formula]['Ef']
        elems = [str(el) for el in target_comp.elements]

        # Convert Ef to grand potential
        if 'O' in elems:
            n_O = target_comp.fractional_composition.as_dict()['O']
            dG = n_O*get_chempot_correction('O', T, p_O2)
            Ef -= dG

        # Rely on hull energy for gaseous species
        if Composition(formula).reduced_formula in ['CO2', 'O2', 'H2O', 'NH3', 'H3N']:
            return None

    else:
        Ef = None

    return Ef

def make_phase_diagram(entries, elems, temp, atmos='air', data_path='arrows/energetics/MP_Energetics.json'):
    """
    Build a pymatgen phase diagram object in the chemical space
    contained by elems, with energies calculated under
    the specified temperature and atmosphere.

    Args:
        entries (list): computed entries for all
            phases in the space.
        elems (list): elements bounding the space.
        temp (int/float): temperature.
        atmos (str): air or inert.
        data_path (str): relative path to
            the file containing the 0 K
            energies from MP.
    Returns:
        final_diag: pymatgen phase diagram
            containing all entries.
    """

    # Round temperature to nearest 100 C
    T = round(temp, -2)

    # Partial pressures of gaseous species in air
    if atmos == 'air':
        p_O2 = 21200
        p_N2 = 79033.
        p_CO2 = 4050
        p_NH3 = 16.
        p_H2O = 2300.
        p_NO = 0.1

    # Estimation based on 1e-6 partial pressure
    elif atmos == 'inert':
        p_O2 = 0.1
        p_N2 = 0.1
        p_CO2 = 0.1
        p_NH3 = 0.1
        p_H2O = 0.1
        p_NO = 0.1

    else:
        raise Exception('Atmosphere must either be air or inert')

    # Construct phase diagram at 0 K
    zero_temp_pd = pd.PhaseDiagram(entries)

    if T == 0:
        return zero_temp_pd

    # Load energies from MP + Bartel
    with open(data_path) as fname:
        energy_data = json.load(fname)

    # Apply temperature corrections via Bartel method
    revised_entries = []
    for entry in entries:
        entry_dict = entry.as_dict()
        ordered_formula = entry.composition.alphabetical_formula
        condensed_formula = ordered_formula.replace(' ', '')
        if condensed_formula in energy_data[str(T)].keys():
            E0 = energy_data['0'][condensed_formula]['Ef']
            Ef = energy_data[str(T)][condensed_formula]['Ef']
            dE = entry.composition.num_atoms*(Ef - E0)
            entry_dict['energy'] += dE
        revised_entries.append(pd.PDEntry.from_dict(entry_dict))

    # Account for temperature-dependence of common gaseous species
    if 'O' in elems:

        # O2
        O_comp = Composition('O')
        O_energ = zero_temp_pd.get_hull_energy(O_comp)
        O_energ += 2*get_chempot_correction('O', T, p_O2)
        O_entry = pd.PDEntry(O_comp, O_energ)
        revised_entries.append(O_entry)

        # CO2
        if 'C' in elems:
            CO2_comp = Composition('CO2')
            CO2_energ = zero_temp_pd.get_hull_energy(CO2_comp)
            CO2_energ += 3*get_chempot_correction('CO2', T, p_CO2)
            CO2_entry = pd.PDEntry(CO2_comp, CO2_energ)
            revised_entries.append(CO2_entry)

        # H2O
        if 'H' in elems:
            H2O_comp = Composition('H2O')
            H2O_energ = zero_temp_pd.get_hull_energy(H2O_comp)
            H2O_energ += 3*get_chempot_correction('H2O', T, p_H2O)
            H2O_entry = pd.PDEntry(H2O_comp, H2O_energ)
            revised_entries.append(H2O_entry)

        # NO
        if 'N' in elems:
            NO_comp = Composition('NO')
            NO_energ = zero_temp_pd.get_hull_energy(NO_comp)
            NO_energ += 2*get_chempot_correction('NO', T, p_NO)
            NO_entry = pd.PDEntry(NO_comp, NO_energ)
            revised_entries.append(NO_entry)

    if 'N' in elems:

        # N2
        N_comp = Composition('N')
        N_energ = zero_temp_pd.get_hull_energy(N_comp)
        N_energ += 2*get_chempot_correction('N', T, p_N2)
        N_entry = pd.PDEntry(N_comp, N_energ)
        revised_entries.append(N_entry)

        if 'H' in elems:

            # NH3
            NH3_comp = Composition('NH3')
            NH3_energ = zero_temp_pd.get_hull_energy(NH3_comp)
            NH3_energ += 4*get_chempot_correction('NH3', T, p_NH3)
            NH3_entry = pd.PDEntry(NH3_comp, NH3_energ)
            revised_entries.append(NH3_entry)

    final_diag = pd.PhaseDiagram(revised_entries)

    return final_diag

def get_pd_dict(available_precursors, temperatures, atmos='air'):
    """
    Build a dictionary of phase diagrams for the given chemical
    space throughout a range of temperatures.

    Args:
        available_precursors (list): a list of chemical formulae
            corresponding to available precursors.
        temp (list): temperatures at which to build the
            phase diagrams.
        atmos (str): air or inert.
    Returns:
        pd_dict: a dictionary of pymatgen phase
            diagram objects, one at each temperature.
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
    """
    Make your own PDEntry from a given structure
    and energy (0 K).

    Args:
        struc: pymatgen structure object.
        energy (float): calculated energy.
    Returns:
        PDEntry for that structure/energy.
    """

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

    # Conversion factor
    EV_TO_KJ_PER_MOL = 96.4853

    # Common gaseous species
    if element not in ['O', 'N', 'Cl', 'F', 'H', 'CO2', 'NH3', 'H2O', 'NO']:
        return 0

    # Constants
    std_temp = 298.15
    std_pres = 1E5
    ideal_gas_const = 8.3144598

    # Cp and S at standard state in J/(K.mol)
    # Data from https://janaf.nist.gov/tables
    Cp_dict = {'O': 29.376,
               'N': 29.124,
               'Cl': 33.949,
               'F': 31.302,
               'H': 28.836,
               'CO2': 37.129,
               'NH3': 35.640,
               'H2O': 33.22,
               'NO': 29.86}
    S_dict = {'O': 205.147,
              'N': 191.609,
              'Cl': 223.079,
              'F': 202.789,
              'H': 130.680,
              'CO2': 213.79,
              'NH3': 192.80,
              'H2O': 194.10,
              'NO': 210.76}

    # Some math
    Cp_std = Cp_dict[element]
    S_std = S_dict[element]
    PV_correction = ideal_gas_const * temp * np.log(pres / std_pres)
    TS_correction = - Cp_std * (temp * np.log(temp) - std_temp * np.log(std_temp)) \
        + Cp_std * (temp - std_temp) * (1 + np.log(std_temp)) \
        - S_std * (temp - std_temp)

    # Total free energy change
    dG = PV_correction + TS_correction

    # Convert to eV/molecule unit
    dG /= 1000 * EV_TO_KJ_PER_MOL

    # Normalize by number of atoms in the gas molecule
    if element == 'H2O':
        dG /= 3
    if element == 'CO2':
        dG /= 3
    if element == 'NH3':
        dG /= 4
    if element == 'NO':
        dG /= 2
    if element == 'O':
        dG /= 2
    if element == 'N':
        dG /= 2

    return dG
