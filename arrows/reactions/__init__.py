from pymatgen.core.composition import Composition
from arrows.reactions import balancer
from arrows import energetics
from itertools import combinations


def get_balanced_coeffs(reactants, products):
    """
    Get stoichiometric coefficients from the chemical
    reaction involving the specified reactants and products.

    Args:
        reactants (list): chemical formulae of the reactants.
        products (list): chemical formulae of the products.
    Returns:
        coeffs (list): stoichiometric coefficients of the
            reactants and products. If reaction cannot be
            balanced, a message (str) will be returned
            informing the user why it cannot be balanced.
    """

    coeffs = balancer.main(reactants, products)

    return coeffs


def get_rxn_energy(reactants, products, temp, cmpd_pd):
    """
    Calculate the reaction energy associated with the
    transformation from reactants to products.

    Args:
        reactants (list): chemical formulae of the reactants.
        products (list): chemical formulae of the products.
        temp (int/float): temperature.
        cmpd_pd (pymatgen phase diagram): phase diagram.
    Returns:
        Rxn energy in meV/atom.
    """

    # Convert products to list, if needed
    if isinstance(products, str):
        products = [products]

    # Get balanced coefficients for rxn
    full_coeffs = get_balanced_coeffs(reactants, products)

    # If rxn cannot be balanced, raise Exception
    if isinstance(full_coeffs, str):
        raise Exception(full_coeffs)

    # Average energy of reactants
    starting_sum = 0.0
    starting_energy = 0.0
    for (formula, coeff) in zip(reactants, full_coeffs[0]):

        # Get formation energy normalized per atom
        Ef = energetics.get_entry_Ef(formula, temp)
        if Ef is None: # If no MP entry exists
            Ef = energetics.get_hull_Ef(formula, cmpd_pd) # Use hull energy

        # Normalize formation energy per formula unit
        comp_dict = Composition(formula).as_dict()
        sum_elems = sum(comp_dict.values())
        starting_sum += coeff*sum_elems
        Ef = Ef*sum_elems

        starting_energy += coeff * Ef

    # Energy normalized per atom
    starting_energy /= starting_sum

    # Average energy of products
    end_sum = 0.0
    final_energy = 0.0
    for (formula, coeff) in zip(products, full_coeffs[1]):

        # Get formation energy normalized per atom
        Ef = energetics.get_entry_Ef(formula, temp)
        if Ef is None: # If no MP entry exists
            Ef = energetics.get_hull_Ef(formula, cmpd_pd) # Use hull energy

        # Normalize formation energy per formula unit
        comp_dict = Composition(formula).as_dict()
        sum_elems = sum(comp_dict.values())
        end_sum += coeff*sum_elems
        Ef = Ef*sum_elems

        final_energy += coeff * Ef

    # Energy normalized per atom
    final_energy /= end_sum

    # Return reaction energy in meV/atom
    return 1000*(final_energy - starting_energy)


def get_dG(initial_cmpds, initial_amounts, targets, allowed_byproducts, open_sys, pd_dict, temp):
    """
    Similar to get_rxn_energy, except this function allows
    linearly dependent precursors to be treated by simply averaging
    over their energies weighted by their amounts.

    Args:
        initial_cmpds (list): starting materials.
        initial_amounts (list): stoichiometric coefficients
            of the starting materials.
        target (list): target composition(s).
        allowed_byproducts (list): phases that may be
            included as secondary products.
        open_sys (bool): whether to allow oxidation.
        pd_dict (dict): a dictionary containing the phase
            diagrams at varied temperature.
        temp (int/float): current temperature.
    Returns:
        final_products (list): phases that are actually
            produced by the reaction.
        dG (float): assocaited rxn energy.
    """

    # Convert products to list, if needed
    if isinstance(targets, str):
        targets = [targets]

    # Phase diagram at specified temperature
    cmpd_pd = pd_dict[temp]

    # Get average energy of starting materials
    net_comp = {}
    total_atoms = 0.0
    initial_energy = 0.0
    for (formula, amount) in zip(initial_cmpds, initial_amounts):

        # Always used reduce formula
        formula = Composition(formula).reduced_formula

        # Get formation energy normalized per atom
        Ef = energetics.get_entry_Ef(formula, temp)
        if Ef is None: # If no MP entry exists
            Ef = energetics.get_hull_Ef(formula, cmpd_pd) # Use hull energy

        # Energy weighted by amount
        num_atoms = Composition(formula).num_atoms
        initial_energy += amount*num_atoms*Ef
        total_atoms += num_atoms

        # Build average composition dictionary
        comp_dict = Composition(formula).as_dict()
        for elem in comp_dict.keys():
            if elem in net_comp.keys():
                net_comp[elem] += amount*comp_dict[elem]
            else:
                net_comp[elem] = amount*comp_dict[elem]

    # Average energy of the precursors, normalized per atom
    initial_energy /= total_atoms

    # Chemical formula for the average composition
    avg_formula = Composition.from_dict(net_comp).reduced_formula

    # Check for balance with target(s) + byproduct(s)
    final_soln, final_products = None, None
    trial_soln = get_balanced_coeffs(targets, [avg_formula])
    if not isinstance(trial_soln, str):
        final_soln = trial_soln.copy()
        final_products = targets.copy()
    else:
        allowed_byproducts = [Composition(cmpd).reduced_formula for cmpd in allowed_byproducts]
        allowed_byproducts += ['O2', 'CO2'] # Allow gaseous evolution
        allowed_byproducts = list(set(allowed_byproducts))
        for num_byp in range(1, len(allowed_byproducts) + 1):
            possible_byproducts = combinations(allowed_byproducts, num_byp)
            for byp_set in possible_byproducts:
                all_products = targets + list(byp_set)
                trial_soln = get_balanced_coeffs(all_products, [avg_formula])
                if not isinstance(trial_soln, str): # If reaction can be balanced
                    final_soln = trial_soln.copy()
                    final_products = all_products.copy()

    # Enumerate through possible gaseous reactants
    # Use those that produce a balanced rxn
    # Here, we only consider O2 or CO2 uptake
    gaseous_reacs = []
    trial_soln = get_balanced_coeffs(targets, [avg_formula, 'O2'])
    if not isinstance(trial_soln, str):
        final_soln = trial_soln.copy()
        final_products = targets.copy()
        gaseous_reacs.append('O2')
    if final_soln == None:
        if open_sys:
            gaseous_byproducts = ['O2', 'CO2']
            for num_byp in range(1, len(gaseous_byproducts) + 1):
                possible_byproducts = combinations(gaseous_byproducts, num_byp)
                for byp_set in possible_byproducts:
                    all_products = targets + list(byp_set)
                    if 'O2' not in all_products:
                        trial_soln = get_balanced_coeffs(all_products, [avg_formula, 'O2'])
                        if not isinstance(trial_soln, str): # If reaction can be balanced
                            final_soln = trial_soln.copy()
                            final_products = all_products.copy()
                            gaseous_reacs.append('O2')

    trial_soln = get_balanced_coeffs(targets, [avg_formula, 'CO2'])
    if not isinstance(trial_soln, str):
        final_soln = trial_soln.copy()
        final_products = targets.copy()
        gaseous_reacs.append('CO2')
    if final_soln == None:
        if open_sys:
            gaseous_byproducts = ['O2', 'CO2']
            for num_byp in range(1, len(gaseous_byproducts) + 1):
                possible_byproducts = combinations(gaseous_byproducts, num_byp)
                for byp_set in possible_byproducts:
                    all_products = targets + list(byp_set)
                    if 'CO2' not in all_products:
                        trial_soln = get_balanced_coeffs(all_products, [avg_formula, 'CO2'])
                        if not isinstance(trial_soln, str): # If reaction can be balanced
                            final_soln = trial_soln.copy()
                            final_products = all_products.copy()
                            gaseous_reacs.append('CO2')

    trial_soln = get_balanced_coeffs(targets, [avg_formula, 'O2', 'CO2'])
    if not isinstance(trial_soln, str):
        final_soln = trial_soln.copy()
        final_products = targets.copy()
        gaseous_reacs.append('CO2')
        gaseous_reacs.append('O2')

    # If no combination of gaseous byproducts enable a balanced rxn, raise an Error
    assert final_soln != None, 'Precursor update went wrong. Reaction cannot be balanced with targets'

    net_comp = {}
    total_atoms = 0.0
    final_energy = 0.0
    for (formula, amount) in zip(final_products, final_soln[0]):

        # Always use reduced formula
        formula = Composition(formula).reduced_formula

        # Get formation energy normalized per atom
        Ef = energetics.get_entry_Ef(formula, temp)
        if Ef is None: # If no MP entry exists
            Ef = energetics.get_hull_Ef(formula, cmpd_pd) # Use hull energy

        # Energy weighted by amount
        num_atoms = Composition(formula).num_atoms
        final_energy += amount*num_atoms*Ef
        total_atoms += num_atoms

    if len(gaseous_reacs) > 0:
        # Subtract energy from gaseous reactant(s)
        for formula, amount in zip(gaseous_reacs, final_soln[1][1:]):
            Ef = energetics.get_entry_Ef(formula, temp)
            if Ef is None: # If no MP entry exists
                Ef = energetics.get_hull_Ef(formula, cmpd_pd) # Use hull energy
            num_atoms = Composition(formula).num_atoms
            final_energy -= amount*num_atoms*Ef

    # Average energy of the products, normalized per atom
    final_energy /= total_atoms

    # Change in energy (meV/atom) from precursors to products
    dG = 1000*(final_energy - initial_energy)

    return final_products, dG
