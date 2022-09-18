from pymatgen.core.composition import Composition
from pymatgen.analysis import interface_reactions
from arrows import reactions
from itertools import combinations
import numpy as np
import csv
import sys


def update_set(precursors, amounts, phase_diagram, temp, cutoff_energ, return_rxn=False):
    """
    Update a set of precursors based on the pairwise reaction with
    the highest TD driving force.
    """

    # Make dictionary for precursor amounts
    precursor_amounts = {}
    for (cmpd, coeff) in zip(precursors, amounts):
        precursor_amounts[cmpd] = coeff

    # All possible pairs of precursors
    all_pairs = combinations(precursors, 2)

    # Get most favorable reaction at each interface
    all_likely_rxns = []
    for pair in all_pairs:
        # Exclude reactions with gases
        if ('O2' not in pair) and ('CO2' not in pair):
            c1 = Composition(pair[0])
            c2 = Composition(pair[1])
            ir = interface_reactions.InterfacialReactivity(c1, c2, phase_diagram)
            kinks = list(ir.get_kinks())
            energies = [info[4] for info in kinks]
            all_likely_rxns.append(kinks[np.argmin(energies)])

    # Get reaction with largest TD driving force
    energies = [info[4] for info in all_likely_rxns]
    top_rxn = all_likely_rxns[np.argmin(energies)]
    driving_force = top_rxn[4]

    if driving_force < cutoff_energ:

        # Get pair of compounds that are expected to react
        reactants = [cmpd.reduced_formula for cmpd in top_rxn[3].reactants]
        products = [cmpd.reduced_formula for cmpd in top_rxn[3].products]

        if return_rxn:
            return [reactants, products, driving_force]

        # Get available amounts of each compound
        avail_amounts = []
        for cmpd in reactants:
            avail_amounts.append(precursor_amounts[cmpd])

        # Amount of each reactant consumed and product produced
        req_amounts = [abs(val) for val in top_rxn[3].coeffs if val < 0] # Sum to 1
        product_amounts = [abs(val) for val in top_rxn[3].coeffs if val > 0]

        # Calculate changes in reactant amounts and add in new products
        leftover_cmpds, leftover_amounts = calculate_amounts(reactants, avail_amounts, req_amounts, products, product_amounts)
        leftover_cmpds = [Composition(cmpd).reduced_formula for cmpd in leftover_cmpds]

        # Add rxn products to final set
        final_set, final_amounts = [], []
        for (cmpd, amount) in zip(leftover_cmpds, leftover_amounts):
            final_set.append(cmpd)
            if (cmpd in precursors) and (cmpd not in reactants):
                combined_amount = amount + precursor_amounts[cmpd]
                final_amounts.append(combined_amount)
            else:
                final_amounts.append(amount)

        # Add non-participating cmpds to final set
        for cmpd in precursors:
            if (cmpd not in reactants) and (cmpd not in leftover_cmpds):
                final_set.append(cmpd)
                final_amounts.append(precursor_amounts[cmpd])

        return final_set, final_amounts

    else:

        if return_rxn:
            return None
        else:
            return None, None


def calculate_amounts(reactants, avail_amounts, req_amounts, products, product_amounts):
    """
    Used to update compounds involved in a pairwise reaction.
    Should not be used to update an entire set!
    For this purpose, update_set() should be used instead.

    Args:
        reactants: cmpds involved in pairwise reaction
        avail_amounts: stoichiometric coefficients of each
        req_amounts: coefficients required for the pairwise reaction
        products: products formed in the pairwise reaction
        product_amounts: coefficients of the products
    Returns:
        final_set: products + remaining reactants
        final_amounts: coefficients for the final_set
    """

    # Make dictionary for reactant amounts
    reactant_amounts = {}
    for (cmpd, coeff) in zip(reactants, avail_amounts):
        reactant_amounts[cmpd] = coeff

    # Normalize so that coeffs sum to 1
    sum_req, sum_avail = sum(req_amounts), sum(avail_amounts)
    req_amounts = [v/sum_req for v in req_amounts]
    avail_amounts = [v/sum_avail for v in avail_amounts]

    # Check if both reactants are fully consumed
    consumed_index = 0
    if np.isclose(avail_amounts, req_amounts, atol=0.01).all():
        remove_cmpds = reactants
        scaling_factor = 1.0

    # Otherwise, check which reactant is consumed first
    else:
        indices = list(range(len(reactants)))
        index = 0
        for (cf_1, cf_2, cmpd) in zip(avail_amounts, req_amounts, reactants):
            if cf_1 < cf_2:
                remove_cmpds = [cmpd]
                scaling_factor = cf_1/cf_2
                consumed_index = index
            index += 1

        indices.remove(consumed_index)
        remaining_index = indices[0]
        remaining_cmpd = reactants[remaining_index]
        amount_consumed = scaling_factor * req_amounts[remaining_index]
        remaining_amount = avail_amounts[remaining_index] - amount_consumed
        actual_amount = remaining_amount * sum_avail # Back to non-normalized amount
        reactant_amounts[remaining_cmpd] = actual_amount # Update amount of remaining cmpd

    # Calculate how much products were produced based on how much reactants were consumed
    amounts_produced = []
    for (cmpd, coeff) in zip(products, product_amounts):
        avail = coeff*sum_avail*avail_amounts[consumed_index]
        req = sum_req*req_amounts[consumed_index]
        amounts_produced.append(avail/req)

    # Update reactant amounts and add new products
    for (cmpd, coeff) in zip(products, amounts_produced):
        if cmpd in reactants:
            reactant_amounts[cmpd] += coeff
        else:
            reactant_amounts[cmpd] = coeff

    # Remove consumed compound(s)
    remaining_set = list(set(reactants) - set(remove_cmpds))

    # Combine remaining reactants with new products
    final_set = remaining_set + products

    # Tabulate amounts of final products
    final_amounts = []
    for cmpd in final_set:
        final_amounts.append(round(reactant_amounts[cmpd], 3))

    return final_set, final_amounts


def retroanalyze(precursors, initial_amounts, products, final_amounts, pd_dict, temp, allowed_byproducts, open_sys=True, enforce_thermo=False, rxn_database=None, already_probed=[]):
    """
    Given a synthesis result, propose possible reaction pathways.

    Args:
        precursors: starting materials in the synthesis procedure
        initial_amounts: stoichiometric coefficients for precursor amounts
        products: observed reaction products
        final_amounts: stoichiometric coefficients for product amounts
        pd_dict: dictionary containing phase diagram objects in the relevant chemical space
        temp: temperature at which the synthesis was carried out
        open_sys: if True, system is open to gaseous byproducts (O2, CO2)
        enforce_thermo: if True, all proposed rxns must be thermodynamically favorable (dG < 0)
        rxn_database: dataset containing known rxn information
    """

    # Check if this synthesis route has already been probed
    current_precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    current_temp = [int(temp)]
    current_route = set(current_precursors + current_temp)
    if precursors in already_probed:
        return 'Reaction already probed.', None, None, None, []

    # Phase diagram at specified temperature
    phase_diagram = pd_dict[temp]

    # Ensure consistent formatting of chemical formulae
    precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    products = [Composition(cmpd).reduced_formula for cmpd in products]
    allowed_byproducts = [Composition(cmpd).reduced_formula for cmpd in allowed_byproducts]

    # Make dictionary for precursor amounts
    precursor_amounts = {}
    for (cmpd, coeff) in zip(precursors, initial_amounts):
        precursor_amounts[cmpd] = coeff

    # Make dictionary for product amounts
    product_amounts = {}
    for (cmpd, coeff) in zip(products, final_amounts):
        product_amounts[cmpd] = coeff

    intermediates = None
    if not rxn_database.is_empty:

        # Known *local* reactions sorted by temperature
        known_rxns = rxn_database.as_sorted_list(local=True)
        interm_set = None
        for first_rxn in known_rxns:

            # If reaction is known to occur at or below current temp
            if first_rxn[-1] <= temp:

                # Amounts consumed or produced
                pair = [Composition(cmpd).reduced_formula for cmpd in first_rxn[0]]
                prods = [Composition(cmpd).reduced_formula for cmpd in first_rxn[1]]

                # Ensure all reactants are present in the current precursor set
                if set(pair).issubset(set(precursors)):

                    bal_info = reactions.get_balanced_coeffs(pair, prods)

                    # Check for oxidation
                    ind = 0
                    solid_pair = pair.copy()
                    possible_oxidants = [['O2'], ['CO2'], ['O2', 'CO2']]
                    while isinstance(bal_info, str):
                        assert ind < len(possible_oxidants), 'Pairwise rxn (%s) cannot be balanced' % ', '.join(solid_pair)
                        full_pair = solid_pair + possible_oxidants[ind]
                        bal_info = reactions.get_balanced_coeffs(full_pair, prods)
                        pair = full_pair.copy()
                        ind += 1
                    req_amounts, amounts_formed = bal_info[0], bal_info[1]

                    # Available amounts
                    avail_amounts = []
                    for cmpd in pair:
                        if cmpd not in ['O2', 'CO2']:
                            avail_amounts.append(precursor_amounts[cmpd])
                        else:
                            avail_amounts.append(1000.0) # Unlimited

                    # Calculate changes in reactant amounts and add in new products
                    leftover_cmpds, leftover_amounts = calculate_amounts(pair, avail_amounts, req_amounts, prods, amounts_formed)
                    leftover_cmpds = [Composition(cmpd).reduced_formula for cmpd in leftover_cmpds]

                    # Add back in the compounds that weren't involved in the pairwise reaction
                    interm_set, interm_amounts = leftover_cmpds, leftover_amounts
                    for cmpd in precursors:
                        if (cmpd not in interm_set) and (cmpd not in pair):
                            interm_set.append(cmpd)
                            interm_amounts.append(precursor_amounts[cmpd])

                    intermediates = []
                    initial_amounts = []
                    for cmpd, amt in zip(interm_set, interm_amounts):
                        if cmpd not in ['O2', 'CO2', 'H3N', 'H2O']:
                            intermediates.append(cmpd)
                            initial_amounts.append(amt)

                    # Consider intermediates as new precursor set
                    precursors = intermediates.copy()

                    # Make dictionary for updated precursor amounts
                    for (cmpd, coeff) in zip(precursors, initial_amounts):
                        precursor_amounts[cmpd] = coeff

    # Check again if this synthesis route has already been probed
    current_precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    current_temp = [int(temp)]
    current_route = set(current_precursors + current_temp)
    if current_route in already_probed:
        return 'Reaction already probed.', None, None, None, []

    # For for inert pairs of compounds
    inert_pairs = []
    possible_pairs = combinations(precursors, 2)
    for pair in possible_pairs:
        if set(pair).issubset(set(products)):
            if ('O2' not in pair) and ('CO2' not in pair):
                inert_pairs.append(frozenset(pair))

    # Check for existing compounds with increased amounts
    amount_changes = {}
    tol = 0.01 # Allow some tolerance
    for (start_cmpd, start_coeff) in zip(precursors, initial_amounts):
        for (end_cmpd, end_coeff) in zip(products, final_amounts):
            if start_cmpd == end_cmpd:
                if end_coeff > (start_coeff + tol):
                    amount_changes[start_cmpd] = end_coeff - start_coeff

    # Check for new compounds produced
    novel_cmpds = list(set(products) - set(precursors))
    for cmpd in novel_cmpds:
        amount_changes[cmpd] = product_amounts[cmpd]

    # Net products minus gaseous species
    observ_products = list(amount_changes.keys())
    observ_products = [cmpd for cmpd in observ_products if cmpd not in ['O2', 'CO2']]

    # If precursors == products, no further analysis necessary
    if len(observ_products) == 0:
        if intermediates == None:
            return 'No reactions occured', None, None, None, inert_pairs
        else:
            return 'Only known intermediate reactions occured', None, None, intermediates, inert_pairs

    # Otherwise, explore suspected reaction pathways
    else:

        # Possible pairwise reactants, including precursors and products
        interm_sets = []
        all_cmpds = list(precursors) + list(observ_products)
        for cmpd in all_cmpds:
            # Check for decomposition (single reactant)
            interm_sets.append([cmpd])
        interm_sets += list(combinations(all_cmpds, 2))
        interm_sets = [list(pair) for pair in interm_sets]
        # Gaseous species may participate in pairwise rxns
        if open_sys:
            for solid_set in interm_sets.copy():
                w_O2 = solid_set + ['O2']
                w_CO2 = solid_set + ['CO2']
                w_both = solid_set + ['O2', 'CO2']
                interm_sets.append(w_O2)
                interm_sets.append(w_CO2)
                interm_sets.append(w_both)

        # Possible sets of products formed via the reactions above
        product_sets = []
        for n in range(1, len(observ_products) + 1):
            product_sets += list(combinations(observ_products, n))
        product_sets = [list(pair) for pair in product_sets]
        # Gaesous byproducts may evolve from rxns
        for solid_set in product_sets.copy():
            w_O2 = solid_set + ['O2']
            w_CO2 = solid_set + ['CO2']
            w_both = solid_set + ['O2', 'CO2']
            product_sets.append(w_O2)
            product_sets.append(w_CO2)
            product_sets.append(w_both)
        # Include allowed byproducts
        solid_byproducts = list(set(allowed_byproducts) - {'O2', 'CO2'}) # Already included
        byproduct_sets = []
        for n in range(1, len(solid_byproducts) + 1):
            byproduct_sets += list(combinations(solid_byproducts, n))
        for existing_set in product_sets.copy():
            for byp_set in byproduct_sets:
                w_byp = list(existing_set) + list(byp_set)
                product_sets.append(w_byp)

        # Check for compostional balance between reactant pairs and product sets
        suspected_rxns = []
        for prod_set in product_sets:
            for pair in interm_sets:
                check_bal = reactions.get_balanced_coeffs(list(pair), list(prod_set))
                if not isinstance(check_bal, str): # If balanced
                    if enforce_thermo:
                        energ = reactions.get_rxn_energy(list(pair), list(prod_set), temp, phase_diagram)
                        if set(pair) != set(prod_set):
                            if energ < 0: # Only include TD favorable rxns
                                if frozenset(pair) not in rxn_database.inert_pairs(temp): # Exclude pairs that do not react at current T
                                    suspected_rxns.append([list(pair), list(prod_set)])
                    else:
                        if frozenset(pair) not in rxn_database.inert_pairs(temp): # Exclude pairs that do not react at current T
                            if set(pair) != set(prod_set):
                                suspected_rxns.append([list(pair), list(prod_set)])

        # Tabulate which compounds are consumed or produced
        consumed_precursors, pairwise_products = [], []
        for rxn_info in suspected_rxns:
            consumed_precursors += rxn_info[0]
            pairwise_products += rxn_info[1]

        # Products not yet accounted for by any pairwise reactions
        mystery_products = list(set(observ_products) - set(pairwise_products))

        # Check which products may be formed in multiple ways
        redundant_products = []
        for cmpd in observ_products:
            if pairwise_products.count(cmpd) > 1:
                redundant_products.append(cmpd)

        # Products with known origin, excluding gaseous phases
        known_products = list(set(pairwise_products) - set(redundant_products) - {'O2', 'CO2'})

        # Finalize and return messages to user
        if (len(mystery_products) == 0) and (len(redundant_products) == 0):
            mssg = 'Reaction pathway fully determined.'
            return mssg, suspected_rxns, known_products, intermediates, inert_pairs
        elif len(known_products) > 0:
            mssg = 'Reaction pathway partially determined.'
            return mssg, suspected_rxns, known_products, intermediates, inert_pairs
        elif len(known_products) > 0:
            mssg = 'Inert pairs discovered.'
            return mssg, suspected_rxns, None, intermediates, inert_pairs
        else:
            mssg = 'No reactions discovered'
            return mssg, suspected_rxns, None, intermediates, inert_pairs


def sort_by_two(data, key_1=-1, key_2=-2):
    """
    Sort zipped list (data) by key_1, then key_2.
    Both key_1 and key_2 should be ints, respresenting
    the indicies of the items to be sorted by.
    """

    first_keys = list(set([vec[key_1] for vec in data]))

    sorted_data = []
    for val in sorted(first_keys):
        subgroup = [vec for vec in data if vec[key_1] == val]
        sorted_data += sorted(subgroup, key=lambda x: x[key_2])

    return sorted_data


def inform_user(temp, precursors, products, amounts, mssg, sus_rxn_info, known_products, intermediates):
    print('\nTemperature: %s' % temp)
    precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    products = [Composition(cmpd).reduced_formula for cmpd in products]
    print('Precursors: %s' % ', '.join(precursors))
    if intermediates != None:
        intermediates = [Composition(cmpd).reduced_formula for cmpd in intermediates]
        intermediates = list(set(intermediates) - {'O2', 'CO2', 'H3N', 'H2O'})
        print('Intermediates: %s' % ', '.join(intermediates))
    print('Products: %s' % ', '.join(products))
    amounts = [str(round(amt, 2)) for amt in amounts]
    print('Amounts: %s' % ', '.join(amounts))
    print(mssg)
    if sus_rxn_info != None:
        print('Suspected reactions:')
        for sus_rxn in sus_rxn_info:
            reacs = [Composition(cmpd).reduced_formula for cmpd in sus_rxn[0]]
            prods = [Composition(cmpd).reduced_formula for cmpd in sus_rxn[1]]
            print('%s == %s' % (' + '.join(reacs), ' + '.join(prods)))
    if known_products != None:
        known_products = [Composition(cmpd).reduced_formula for cmpd in known_products]
        print('Products with known origin: %s' % ', '.join(known_products))


class rxn_database:

    """
    Note:
    Inert (non-reaction) temperatures are strict (rxn occurs > T)
    Whereas reaction temperatures are non-strict (rxn occurs <= T)
    In other words, lower T < rxn T <= upper T
    """

    def __init__(self):
        self.known_rxns = {}

    def load(self, filepath='PairwiseRxns.csv'):

        with open(filepath) as f:

            # Iterate through each line
            for i, line in enumerate(f.readlines()):

                # Skip header
                if i == 0:
                    continue

                # Load pairwise reactants
                reacs = line.split(',')[0].split(' + ')
                reacs = frozenset([Composition(cmpd).reduced_formula for cmpd in reacs])

                # Load associated products (if any)
                prods = line.split(', ')[1].split(' + ')
                if 'None' not in prods:
                    prods = frozenset([Composition(cmpd).reduced_formula for cmpd in prods])
                else:
                    prods = None

                # Extreme bounds
                lower_T, upper_T = 0, 2000

                # Load temperature data
                mssg = line.split(',')[2]
                if 'Reacts between' in mssg:
                    T_range = mssg.split()[-2]
                    lower_T = int(T_range.split('-')[0])
                    upper_T = int(T_range.split('-')[1])
                elif 'Reacts below' in mssg:
                    upper_T = int(mssg.split()[-2])
                elif 'Does not react' in mssg:
                    lower_T = int(mssg.split()[-2])

                # Add rxn data to dictionary
                self.known_rxns[reacs] = [[prods, [lower_T, upper_T], 'Global']]

    def update(self, mssg, sus_rxn_info, known_products, inert_pairs, temp):

        # Check for updates
        is_updated = False

        # Set lower bounds on rxn temperatures
        for reacs in inert_pairs:

            # Check inert pairs may take part in a reaction
            # that was not yet compelte (reactants leftover)
            all_sus_reacs = []
            for sus_rxn in sus_rxn_info:
                sus_reacs = frozenset([Composition(cmpd).reduced_formula for cmpd in sus_rxn[0] if \
                    Composition(cmpd).reduced_formula not in ['O2', 'CO2']])
                all_sus_reacs.append(sus_reacs)

            reacs = frozenset([Composition(cmpd).reduced_formula for cmpd in list(reacs) if \
                Composition(cmpd).reduced_formula not in ['O2', 'CO2']])

            if reacs not in all_sus_reacs:

                # Check if inert pair is already known
                if reacs in self.known_rxns.keys():
                    for i, report in enumerate(self.known_rxns[reacs]):
                        # Only update if temperature falls within existing bounds
                        if (temp > self.known_rxns[reacs][i][1][0]) and (temp < self.known_rxns[reacs][i][1][1]):
                            self.known_rxns[reacs][i][1][0] = temp
                            is_updated = True
                else:
                    # Use 2,000 as a hard upper limit on all rxns
                    self.known_rxns[reacs] = [[None, [temp, 2000], 'Local']]
                    is_updated = True

        if mssg == 'Reaction pathway fully determined.':

            # Set upper bounds on rxn temperatures
            for sus_rxn in sus_rxn_info:

                # Use frozenset; order does not matter; hashable
                reacs = frozenset([Composition(cmpd).reduced_formula for cmpd in sus_rxn[0] if \
                    Composition(cmpd).reduced_formula not in ['O2', 'CO2']])
                prods = frozenset([Composition(cmpd).reduced_formula for cmpd in sus_rxn[1]])

                # Check if any info is available for these reactants
                new_products = True
                if reacs in self.known_rxns.keys():

                    for i, report in enumerate(self.known_rxns[reacs]):

                        # If products are new, update the entry
                        if report[0] is None:
                            new_products = False
                            self.known_rxns[reacs][i][2] = 'Local'
                            self.known_rxns[reacs][i][0] = prods
                            self.known_rxns[reacs][i][1][1] = temp
                            is_updated = True

                        else:

                            if report[0] == prods:

                                new_products = False

                                # Only update temperature if new T < old T
                                if temp < self.known_rxns[reacs][i][1][1]:
                                    self.known_rxns[reacs][i][2] = 'Local'
                                    self.known_rxns[reacs][i][1][1] = temp
                                    is_updated = True

                # If these reactants are new, add them to the database
                else:
                    new_products = False
                    self.known_rxns[reacs] = [[prods, [0, temp], 'Local']]
                    is_updated = True

                # If reactants are known but products are new, create a new report
                if new_products:
                    inert_bound = 0
                    for report in self.known_rxns[reacs]:
                        # Use previous reports regarding inert pairs to find lower bound
                        if (report[0] is None) and (report[1][0] > inert_bound):
                            inert_bound =  report[1][0]
                        # Upper bound of previous rxn may be used as lower bound of current one
                        elif (report[0] is not None) and (report[1][1] > inert_bound):
                            inert_bound =  report[1][1]
                    self.known_rxns[reacs].append([prods, [inert_bound, temp], 'Local'])
                    is_updated = True

        if mssg == 'Reaction pathway partially determined.':

            # Set upper bounds on rxn temperatures
            for sus_rxn in sus_rxn_info:

                # Use frozenset; order does not matter; hashable
                reacs = frozenset([Composition(cmpd).reduced_formula for cmpd in sus_rxn[0] if \
                    Composition(cmpd).reduced_formula not in ['O2', 'CO2']])
                prods = frozenset([Composition(cmpd).reduced_formula for cmpd in sus_rxn[1]])

                # Check if suspected reaction produces known phase
                for known_phase in known_products:

                    known_phase = Composition(known_phase).reduced_formula

                    # If so, this reaction is reliable. Add it to rxn database
                    if known_phase in prods:

                        # Check if any info is available for these reactants
                        new_products = True
                        if reacs in self.known_rxns.keys():

                            for i, report in enumerate(self.known_rxns[reacs]):

                                # If products are new, update the entry
                                if report[0] is None:
                                    new_products = False
                                    self.known_rxns[reacs][i][2] = 'Local'
                                    self.known_rxns[reacs][i][0] = prods
                                    self.known_rxns[reacs][i][1][1] = temp
                                    is_updated = True

                                else:

                                    if report[0] == prods:

                                        new_products = False

                                        # Only update temperature if new T < old T
                                        if temp < self.known_rxns[reacs][i][1][1]:
                                            self.known_rxns[reacs][i][2] = 'Local'
                                            self.known_rxns[reacs][i][1][1] = temp
                                            is_updated = True

                        # If these reactants are new, add them to the database
                        else:
                            new_products = False
                            self.known_rxns[reacs] = [[prods, [0, temp], 'Local']]
                            is_updated = True

                        # If reactants are known but products are new, create a new report
                        if new_products:
                            inert_bound = 0
                            for report in self.known_rxns[reacs]:
                                # Use previous reports regarding inert pairs to find lower bound
                                if (report[0] is None) and (report[1][0] > inert_bound):
                                    inert_bound =  report[1][0]
                                # Upper bound of previous rxn may be used as lower bound of current one
                                elif (report[0] is not None) and (report[1][1] > inert_bound):
                                    inert_bound =  report[1][1]
                            self.known_rxns[reacs].append([prods, [inert_bound, temp], 'Local'])
                            is_updated = True

        return is_updated

    def inert_pairs(self, temp):
        non_reacs = []
        for reacs in self.known_rxns.keys():
            all_inert = True
            for report in self.known_rxns[reacs]:
                if report[1][0] < temp:
                    all_inert = False
            if all_inert:
                non_reacs.append(reacs)
        return non_reacs

    def make_global(self):
        for reacs in self.known_rxns.keys():
            for i, report in enumerate(self.known_rxns[reacs]):
                self.known_rxns[reacs][i][2] = 'Global'

    def as_dict(self):
        return self.known_rxns

    def as_sorted_list(self, local=False):
        rxn_list = []
        for reacs in self.known_rxns.keys():
            for i, report in enumerate(self.known_rxns[reacs]):
                if local:
                    if self.known_rxns[reacs][i][2] == 'Local':
                        rxn_list.append([reacs, self.known_rxns[reacs][i][0], self.known_rxns[reacs][i][1][0], self.known_rxns[reacs][i][1][1]])
                else:
                    rxn_list.append([reacs, self.known_rxns[reacs][i][0], self.known_rxns[reacs][i][1][0], self.known_rxns[reacs][i][1][1]])
        # Sort by temperature
        sorted_rxns = sorted(rxn_list, key=lambda x: x[-1])
        return sorted_rxns

    def print_info(self):
        print('\nKnown reactions:')
        for reacs in self.known_rxns.keys():
            for i, report in enumerate(self.known_rxns[reacs]):
                reactants = ' + '.join(reacs)
                if self.known_rxns[reacs][i][0] != None:
                    prods = ' + '.join(self.known_rxns[reacs][i][0])
                else:
                    prods = 'Unknown'
                lower_temp = self.known_rxns[reacs][i][1][0]
                upper_temp = self.known_rxns[reacs][i][1][1]
                rxn_str = '%s == %s @ %s-%s C' % (reactants, prods, lower_temp, upper_temp)
                print(rxn_str)

    def save(self, to='PairwiseRxns.csv'):
        with open(to, 'w+') as datafile:
            csv_writer = csv.writer(datafile)
            csv_writer.writerow(['Pairwise reactants', 'Pairwise Products', 'Temperature Range'])
            for reacs in self.known_rxns.keys():
                for i, report in enumerate(self.known_rxns[reacs]):
                    reactants = [Composition(ph).reduced_formula for ph in reacs]
                    reactants = ' + '.join(reactants)
                    if self.known_rxns[reacs][i][0] != None:
                        products = self.known_rxns[reacs][i][0]
                        products = [Composition(ph).reduced_formula for ph in products]
                        products = ' + '.join(products)
                        products = ' %s' % products
                    else:
                        products = ' None'
                    lower_temp = self.known_rxns[reacs][i][1][0]
                    upper_temp = self.known_rxns[reacs][i][1][1]
                    if (lower_temp > 0.0) and (upper_temp < 2000.0):
                        temp_bounds = ' Reacts between %s-%s C' % (lower_temp, upper_temp)
                    elif (lower_temp > 0.0) and (upper_temp == 2000.0):
                        temp_bounds = ' Does not react at or below %s C' % lower_temp
                    elif (lower_temp == 0.0) and (upper_temp < 2000.0):
                        temp_bounds = ' Reacts below %s C' % upper_temp
                    csv_writer.writerow([reactants, products, temp_bounds])

    @property
    def is_empty(self):
        if len(self.known_rxns) == 0:
            return True
        else:
            return False


def pred_evolution(precursors, initial_amounts, rxn_database, greedy, min_T, allow_oxidation):

    # Placeholder for now
    temp = 1000.0

    # Ensure consistent formatting of chemical formulae
    precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]

    # Make dictionary for precursor amounts
    precursor_amounts = {}
    for (cmpd, coeff) in zip(precursors, initial_amounts):
        precursor_amounts[cmpd] = coeff

    final_cmpds, final_amounts = None, None

    if rxn_database.is_empty:
        sys.stdout.flush()
        final_cmpds, final_amounts = precursors, initial_amounts

    else:

        # Evolve set until we have insufficient rxn information
        while precursors != None:

            # Save for the end
            final_cmpds, final_amounts = precursors.copy(), initial_amounts.copy()

            # Known reactions sorted by temperature
            known_rxns = rxn_database.as_sorted_list()
            interm_set = None

            """
            Check whether:
            1) All pairs are accounted for by global reactions
            2) No reaction temperatures are degenerate

            Unless --greedy is specified, in which case
            low-T rxns are always assumed to occur.
            """

            # Decomposition reactions
            possible_pairs = list(combinations(precursors, 1))

            # Oxidation reactions
            possible_pairs += [(ph, 'O2') for ph in precursors]
            possible_pairs += [(ph, 'CO2') for ph in precursors]

            # Pairwise reactions
            possible_pairs += list(combinations(precursors, 2))

            # Known reactions
            known_pairs = [info[0] for info in known_rxns]
            known_temps = [info[-1] for info in known_rxns]

            found_greedy = False
            all_known, first_rxn, degen = True, None, False
            min_rxn_temp = max(known_temps) + 100
            for pair in possible_pairs:
                if frozenset(pair) not in known_pairs:
                    if not found_greedy and not greedy:
                        # Don't penalize the absence of known oxidation/decomposition rxns
                        # These are already implicit in the known reactions (assuming same atmosphere)
                        # For example: if MnO reacts with Li2O @ 600 C, it takes precedent over oxidation
                        if ('O2' not in pair) and ('CO2' not in pair) and (len(pair) > 1):
                            all_known = False
                else:
                    rxn = known_rxns[known_pairs.index(frozenset(pair))]
                    rxn_temp = rxn[-1]
                    # If rxn temp is equal to the lower bound and greedy is True
                    if (rxn_temp == min_T) and (greedy is True):
                        # If rxn temp is lower than others *and* products are known
                        if (rxn_temp < min_rxn_temp) and (rxn[1] != None):
                            all_known = True
                            found_greedy = True
                            min_rxn_temp = rxn_temp
                            first_rxn = rxn
                            degen = False
                        elif (rxn[-1] == min_rxn_temp) and (rxn[1] != None):
                            degen = True
                    # Otherwise, all pairwise combinations must be known
                    else:
                        rxn = known_rxns[known_pairs.index(frozenset(pair))]
                        rxn_temp = rxn[-1]
                        # If rxn temp is lower than others *and* products are known
                        if (rxn_temp < min_rxn_temp) and (rxn[1] != None):
                            min_rxn_temp = rxn_temp
                            first_rxn = rxn
                            degen = False
                        elif (rxn[-1] == min_rxn_temp) and (rxn[1] != None):
                            degen = True

            # If both criteria are satisfied, update the set accordingly
            if (all_known == True) and (degen == False) and (first_rxn != None):

                # If reaction is known to occur at or below current temp
                if first_rxn[-1] <= temp:

                    # Amounts consumed or produced
                    pair = [Composition(cmpd).reduced_formula for cmpd in first_rxn[0]]
                    prods = [Composition(cmpd).reduced_formula for cmpd in first_rxn[1]]
                    bal_info = reactions.get_balanced_coeffs(pair, prods)

                    # Check for oxidation
                    ind = 0
                    solid_pair = pair.copy()
                    possible_oxidants = [['O2'], ['CO2'], ['O2', 'CO2']]
                    while isinstance(bal_info, str):
                        assert ind < len(possible_oxidants), 'Pairwise rxn (%s) cannot be balanced' % ', '.join(solid_pair)
                        full_pair = solid_pair + possible_oxidants[ind]
                        bal_info = reactions.get_balanced_coeffs(full_pair, prods)
                        pair = full_pair.copy()
                        ind += 1
                    req_amounts, amounts_formed = bal_info[0], bal_info[1]

                    # Available amounts
                    avail_amounts = []
                    for cmpd in pair:
                        if cmpd not in ['O2', 'CO2']:
                            avail_amounts.append(precursor_amounts[cmpd])
                        else:
                            avail_amounts.append(1000.0) # Unlimited

                    # Calculate changes in reactant amounts and add in new products
                    leftover_cmpds, leftover_amounts = calculate_amounts(pair, avail_amounts, req_amounts, prods, amounts_formed)
                    leftover_cmpds = [Composition(cmpd).reduced_formula for cmpd in leftover_cmpds]

                    # Add back in the compounds that weren't involved in the pairwise reaction
                    interm_set, interm_amounts = leftover_cmpds, leftover_amounts
                    for cmpd in precursors:
                        if (cmpd not in interm_set) and (cmpd not in pair):
                            interm_set.append(cmpd)
                            interm_amounts.append(precursor_amounts[cmpd])

                    # Consider intermediates as new precursor set
                    precursors, initial_amounts = [], []
                    for (cmpd, coeff) in zip(interm_set, interm_amounts):
                        cmpd_formula = Composition(cmpd).reduced_formula
                        # Exclude gaseous byproducts
                        if cmpd_formula not in ['O2', 'CO2', 'H3N', 'H2O']:
                            precursors.append(cmpd_formula)
                            initial_amounts.append(coeff)

                    # Make dictionary for updated precursor amounts
                    precursor_amounts = {}
                    for (cmpd, coeff) in zip(precursors, initial_amounts):
                        precursor_amounts[cmpd] = coeff

                # Otherwise, exit loop
                else:
                    precursors = None

            # Otherwise, exit loop
            else:
                precursors, initial_amounts = None, None

    final_cmpds = [Composition(cmpd).reduced_formula for cmpd in final_cmpds]

    return final_cmpds, final_amounts
