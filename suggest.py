from arrows import energetics, reactions, pairwise, exparser
from pymatgen.core.composition import Composition
from itertools import combinations
import numpy as np
import json
import time
import csv
import sys
import os


def load_rxn_data(fname='Rxn_TD.csv', explore=False):
    """
    Loads pre-calculated list of precursors and their
    associated reaction energetics.

    Args:
        fname (str): filename containing rxn info
        explore (bool): whether to prioritize reactions with
            maximal driving force (dG) to form the target phase,
            or to prioritze precursor sets with the most number
            of new interfaces.
    Returns:
        sorted_rxn_info (list): a ranked list of precursor sets
            along with their associated info including stoichiometry,
            predicted intermediates, expected products, interfaces,
            and calculated driving force (dG) at the last known
            rxn step (from the predicted intermediates).
    """

    # Load reaction data
    sorted_rxn_info = []
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file)
        i = 0
        for row in csv_reader:
            if i != 0:
                reactants = row[0].split(' + ')
                reactants = [Composition(cmpd).reduced_formula for cmpd in reactants]
                interfaces = [set(pair) for pair in combinations(reactants, 2)]
                num_interfaces = len(interfaces)
                amounts = [float(v) for v in row[1].split(' + ')]
                products = row[2].split(' + ')
                energ = float(row[3])
                expec_yield = 0.0
                # Reactants saved twice to preserve original info; second may be updated
                sorted_rxn_info.append([reactants, amounts, reactants, amounts, products, expec_yield, interfaces, num_interfaces, energ])
            i += 1

    # Exploration: prioritize no. of new interfaces
    if explore:
        sorted_rxn_info = sorted(sorted_rxn_info, key=lambda x: (-x[-2], x[-1]))

    # Exploitation: prioritize maximal dG
    else:
        sorted_rxn_info = sorted(sorted_rxn_info, key=lambda x: (x[-1], -x[-2]))

    return sorted_rxn_info

def update_ranking(rxn_database, sorted_rxn_info, explore=False):
    """
    Update the ranking of precursor sets based on newly learned
    information in the pairwise reaction database.

    Args:
        rxn_database (dict): a dictionary where each key is a
            pair of reactants, and rxn_database[key] contains
            the expected reaction temperature and products
            associated with those reactants.
        sorted_rxn_info (list): a ranked list of precursor sets
            along with their associated info including stoichiometry,
            predicted intermediates, expected products, interfaces,
            and calculated driving force (dG) at the last known
            rxn step (from the predicted intermediates).
        explore (bool): whether to prioritize reactions with
            maximal driving force (dG) to form the target phase,
            or to prioritze precursor sets with the most number
            of new interfaces.
    Returns:
        sorted_rxn_info (list): a ranked list of precursor sets
            along with their associated info including stoichiometry,
            predicted intermediates, expected products, interfaces,
            and calculated driving force (dG) at the last known
            rxn step (from the predicted intermediates).
    """

    # Evolve all precursor sets using the latest information
    evolved_rxn_info = []
    for rxn_ind, starting_rxn in enumerate(sorted_rxn_info):
        original_set = starting_rxn[0]
        original_amounts = starting_rxn[1]
        starting_materials = starting_rxn[2]
        starting_amounts = starting_rxn[3]
        new_products = starting_rxn[4]
        new_materials, new_amounts = pairwise.pred_evolution(starting_materials, starting_amounts, rxn_database, greedy, min(temps), allow_oxidation)
        if set(new_materials) != set(starting_materials):
            if verbose:
                print('\nPredicted evolution: %s --> %s' % (' + '.join(sorted(starting_materials)), ' + '.join(sorted(new_materials))))
            if target_product in new_materials:
                ind = new_materials.index(target_product)
                expec_yield = new_amounts[ind]
                if not reward_partial_yield:
                    if expec_yield != 1.0:
                        expec_yield = 0.0
                new_products = [target_product]
                energ = 0.0
            else:
                expec_yield = 0.0
                """
                Sometimes, reactions occur that cause you to deviate from the target composition.
                For example, you lose some gaseous species that was meant to participate in the synthesis reaction.
                In such cases, the desired reaction cannot be balanced. These are therefore excluded from any further consideration.
                """
                try:
                    new_products, energ = reactions.get_dG(new_materials, new_amounts, target_product, allowed_byproducts, open_sys, pd_dict, min(temps))
                except:
                    continue
            all_interfaces = [frozenset(pair) for pair in combinations(new_materials, 2)]
            known_interfaces = list(rxn_database.as_dict().keys())
            new_interfaces = set(all_interfaces) - set(known_interfaces)
            new_interfaces = [set(interf) for interf in new_interfaces]
            num_interfaces = len(new_interfaces)
            evolved_rxn_info.append([original_set, original_amounts, new_materials, new_amounts, new_products, expec_yield, new_interfaces, num_interfaces, energ])
        else:
            evolved_rxn_info.append(starting_rxn)

    # Exploration priority: yield, no. of interfaces, dG
    if explore:
        sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (-x[-4], -x[-2], x[-1]))

    # Exploration priority: yield, dG, no. of interfaces
    else:
        sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (-x[-4], x[-1], -x[-2]))

    return sorted_rxn_info


def load_pairwise_rxns(rxn_database, sorted_rxn_info, fname='PairwiseRxns.csv', explore=False):
    """
    Load pairwise reaction from existing csv file and
    use that information to update the precursor ranking.

    Args:
        rxn_database (dict): a dictionary where each key is a
            pair of reactants, and rxn_database[key] contains
            the expected reaction temperature and products
            associated with those reactants.
        sorted_rxn_info (list): a ranked list of precursor sets
            along with their associated info including stoichiometry,
            predicted intermediates, expected products, interfaces,
            and calculated driving force (dG) at the last known
            rxn step (from the predicted intermediates).
        explore (bool): whether to prioritize reactions with
            maximal driving force (dG) to form the target phase,
            or to prioritze precursor sets with the most number
            of new interfaces.
    Returns:
        rxn_database (dict): an *updated* dictionary where each key is
            a pair of reactants, and rxn_database[key] contains
            the expected reaction temperature and products
            associated with those reactants.
        sorted_rxn_info (list): a ranked list of precursor sets
            along with their associated info including stoichiometry,
            predicted intermediates, expected products, interfaces,
            and calculated driving force (dG) at the last known
            rxn step (from the predicted intermediates).
    """


    # Load existing data
    rxn_database.load(filepath=fname)

    # Update rxn databse using existing data
    sorted_rxn_info = update_ranking(rxn_database, sorted_rxn_info, explore=False)

    return rxn_database, sorted_rxn_info

def check_redundancy(interm_phases, interm_wts, known_interm):
    """
    Check whether the current set of intermediate phases and
    their associated weight fractions have been encountered
    previously. If so: reaction pathway is redundant.

    Args:
        interm_phases (list): chemical formulae of the
            observed intermediate phases.
        interm_wts (list): weight fractions associated
            with interm_phases.
        known_interm (dict): a dictionary where each key is
            a set of previously encountered intermediates,
            and known_interm[key] contains their weight fractions.
    Returns:
        redundant (bool): True if the intermediate phases and weights
            have been encountered previously.
    """

    # Check if predicted intermediates have already been sampled
    redundant = False
    if interm_phases in known_interm.keys():
        if known_interm[interm_phases]['Success'] == False:
            for past_wts in known_interm[interm_phases]['Amounts']:
                similar = np.isclose(interm_wts, past_wts, atol=0.1) # 10% wf tolerance
                if False not in similar:
                    redundant = True
        else:
            for past_wts in known_interm[interm_phases]['Amounts']:
                similar = np.isclose(interm_wts, past_wts, atol=0.1) # 10% wf tolerance
                if False not in similar:
                    redundant = True
                    highT_products, highT_amounts = exparser.get_products(precursors, increasing_temps[-1], exp_data)
                    if len(highT_products) == 1:
                        if highT_products[0] == target_product:
                            print('\nRedundant success')
                            print('Products: %s' % target_product)
                            print('Amounts: 1.0')

    return redundant

def print_final_mssg(precursors, T, increasing_temps, sorted_rxn_info, exp_data, batch_size=1, verbose=False):
    """
    Prints a message informing the user any suggested exps
    and closing the script.

    Args:
        precursors (list): a list of suggested precursors
            for the next experimental iteration.
        T (int/float): suggested temperature.
        increasing_temps (list): a list of all temperatures
            that may be tested in this campaign.
        sorted_rxn_info (list): a ranked list of precursor sets
            along with their associated info including stoichiometry,
            predicted intermediates, expected products, interfaces,
            and calculated driving force (dG) at the last known
            rxn step (from the predicted intermediates).
        exp_data (dict): a dictionary containing all
            available experimental results.
        batch_size (int): number of experiments to suggest.
        verbose (bool): whether to print the full ranking
            of precursor sets.
    Returns:
        None
    """

    if verbose:
         print('\nCurrent Ranking:')
         for rxn in sorted_rxn_info:
             print(', '.join(rxn[0]))

    if batch_size == 1:
        print('\n-- Suggested experiment --')
        print('Precursors: %s' % precursors)
        print('Temperature: %s C' % T)

    else:
        print('\n-- Suggested experiments --')
        num_suggest = 0
        print('%s:' % int(num_suggest+1))
        print('Precursors: %s' % precursors)
        print('Temperature: %s C' % T)
        num_suggest += 1
        r_ind = 1
        while (num_suggest < batch_size) and (r_ind < len(sorted_rxn_info)):
            rxn = sorted_rxn_info[r_ind]
            precursors = rxn[0]
            products, final_wts = exparser.get_products(precursors, min(increasing_temps), exp_data)
            if products is None:
                print(int(num_suggest+1))
                print('Precursors: %s' % precursors)
                print('Temperature: %s C' % T)
                num_suggest += 1
            r_ind += 1

    sys.exit()

def update_interm(interm_phases, interm_wts, known_interm, redundant, increasing_temps, exp_data):
    """
    Update the dictionary of known intermediate phases
    based on the most recently observed set of intermediates.

    Args:
        interm_phases (list): chemical formulae of the
            observed intermediate phases.
        interm_wts (list): weight fractions associated
            with interm_phases.
        known_interm (dict): a dictionary where each key is
            a set of previously encountered intermediates,
            and known_interm[key] contains their weight fractions.
        redundant (bool): True if the intermediate phases and weights
            have been encountered previously.
        increasing_temps (list): a list of all temperatures
            that may be tested in this campaign.
        exp_data (dict): a dictionary containing all
            available experimental results.
    Returns:
        known_interm (dict): an *updated* dictionary where each key
            is a set of previously encountered intermediates,
            and known_interm[key] contains their weight fractions.
    """

    if not redundant:
        if interm != None:
            if interm_phases in known_interm.keys():
                known_interm[interm_phases]['Amounts'].append(interm_wts)
                highT_products, highT_amounts = exparser.get_products(precursors, increasing_temps[-1], exp_data)
                if len(highT_products) == 1:
                    if highT_products[0] == target_product:
                        known_interm[interm_phases]['Success'] = True
            else:
                known_interm[interm_phases] = {}
                known_interm[interm_phases]['Amounts'] = [interm_wts]
                known_interm[interm_phases]['Success'] = False
                highT_products, highT_amounts = exparser.get_products(precursors, increasing_temps[-1], exp_data)
                if len(highT_products) == 1:
                    if highT_products[0] == target_product:
                        known_interm[interm_phases]['Success'] = True

    return known_interm

def update_probed_rxns(precursors, T, interm, products, probed_rxns):
    """
    Update a list of previously observed reactions based
    on the most recently tested one.

    Args:
        precursors (list): a list of precursors.
        T (int/float): temperature.
        interm (list): observed intermediate phases.
        products (ilst): observed rxn products.
        probed_rxns (list): previously tested rxns.
    Returns:
        probed_rxns (list): an *updated* list of
            previously tested rxns.
    """

    # Save precursors to probed routes
    current_precursors = [Composition(cmpd).reduced_formula for cmpd in precursors]
    current_temp = [int(T)]
    current_route = set(current_precursors + current_temp)
    probed_rxns.append(current_route)

    # Also save intermediates (if there are any) to probed routes
    if interm != None:
        current_interm = [Composition(cmpd).reduced_formula for cmpd in interm]
        current_temp = [int(T)]
        current_route = set(current_interm + current_temp)
        probed_rxns.append(current_route)

    # Include products as probed routes as well, though this neglects kinetics (finite rxn times)
    current_products = [Composition(cmpd).reduced_formula for cmpd in products]
    current_temp = [int(T)]
    current_route = set(current_products + current_temp)
    probed_rxns.append(current_route)

    return probed_rxns


if __name__ == '__main__':

    verbose = False # Whether to print all info
    explore = False # Default to exploitation
    all = False # Whether to stop when phase pure target obtained
    enforce_thermo = False # Consider pairwise rxns with dG > 0
    greedy = False # Whether to assume low-T rxns occur first
    reward_partial_yield = False # Whether to reward non-pure results
    batch_size = 1 # Number of suggested experiments per batch
    for arg in sys.argv:
        if '--verbose' in arg:
            verbose = True
        if '--explore' in arg:
            explore = True
        if '--all' in arg:
            all = True
        if '--enforce_thermo' in arg:
            enforce_thermo = True
        if '--greedy' in arg:
            greedy = True
        if '--partial_yield' in arg:
            reward_partial_yield = True
        if '--batch' in arg:
            batch_size = int(arg.split('=')[1])

    # Load settings
    with open('Settings.json') as f:
        settings = json.load(f)
    available_precursors = settings['Precursors']
    allow_oxidation = False
    if settings['Allow Oxidation'] == 'True':
        allow_oxidation = True
        available_precursors.append('O2')
        available_precursors.append('CO2')
    target_product = Composition(settings['Target']).reduced_formula
    allowed_byproducts = settings['Allowed Byproducts']
    temps = settings['Temperatures']
    open_sys = settings['Open System']

    # Load experimental rxn data (if any exist)
    if 'Exp.json' in os.listdir('.'):
        with open('Exp.json') as f:
            exp_data = json.load(f)
        exp_data = exp_data['Universal File']
    else:
        print('No experimental data found. Starting from scratch.')
        exp_data = None

    # Build phase diagrams
    pd_dict = energetics.get_pd_dict(available_precursors, temps)

    # Ensure reaction data exists
    assert 'Rxn_TD.csv' in os.listdir('.'), 'No reaction data found. Please run gather_rxns first.'

    # Load reaction data
    sorted_rxn_info = load_rxn_data('Rxn_TD.csv', explore)

    # Pairwise rxn database
    rxn_database = pairwise.rxn_database()

    # Load existing pairwise rxn data
    if 'PairwiseRxns.csv' in os.listdir('.'):
        rxn_database, sorted_rxn_info = load_pairwise_rxns(rxn_database, sorted_rxn_info, 'PairwiseRxns.csv', explore)

    # Temperature ordering (low to high)
    increasing_temps = sorted(temps)

    # Iterate through each rxn
    probed_rxns = []
    known_interm = {}
    num_rxns = len(sorted_rxn_info)
    for i in range(num_rxns):

        # Print info
        sys.stdout.flush()

        # Try rxn with highest dG
        try:
            rxn = sorted_rxn_info[0]

        # Unless all unique reactions have already been sampled
        except IndexError:
            print('\nAll unique reactions sampled.')
            sys.exit()

        # Get predicted intermediates
        interm = sorted(list(zip(rxn[2], rxn[3])))
        interm_phases = tuple([ph[0] for ph in interm])
        interm_coeffs = tuple([ph[1] for ph in interm])

        # Convert stoichiometry to weight fraction
        net_weight = sum([cf*Composition(ph).weight for cf, ph in zip(interm_coeffs, interm_phases)])
        interm_wts = [cf*Composition(ph).weight/net_weight for cf, ph in zip(interm_coeffs, interm_phases)]

        # Check if predicted intermediates have already been sampled
        redundant = check_redundancy(interm_phases, interm_wts, known_interm)

        # Track changes to rxn database
        updated = False

        # Iterate through each temperature
        for T in increasing_temps:

            # Skip if redundant
            if redundant:
                if verbose:
                    print('\nRedundant: %s @ %s C' %  (', '.join(rxn[0]), T))
                continue

            # Starting materials
            precursors, initial_amounts = rxn[0], rxn[1]

            # Parse experimental reaction data
            products, final_wts = exparser.get_products(precursors, T, exp_data)

            # If precursors or temperature not sampled yet, suggest new experiment(s)
            if products is None:
                print_final_mssg(precursors, T, increasing_temps, sorted_rxn_info, exp_data, batch_size, verbose)

            # Formulate intermediates
            interm = sorted(list(zip(products, final_wts)))
            interm_phases = tuple([ph[0] for ph in interm])
            interm_wts = tuple([ph[1] for ph in interm])

            # Check for redundant intermediates tha form at low T
            if T == min(temps):
                redundant = check_redundancy(interm_phases, interm_wts, known_interm)
                known_interm = update_interm(interm_phases, interm_wts, known_interm, redundant, increasing_temps, exp_data)

            # Perform reaction pathway analysis
            mssg, sus_rxn_info, known_products, interm, inert_pairs = pairwise.retroanalyze(precursors, initial_amounts, products, final_wts,
                pd_dict, T, allowed_byproducts, open_sys, enforce_thermo, rxn_database, probed_rxns)

            # If this reaction has already been probed, skip it
            if mssg == 'Reaction already probed.':
                if verbose:
                    print(mssg)
                continue

            # Update probed rxns with new data
            probed_rxns = update_probed_rxns(precursors, T, interm, products, probed_rxns)

            # No need to go further if intermediate reactions are already known
            if mssg == 'Only known intermediate reactions occured.':
                pairwise.inform_user(T, precursors, products, final_wts, mssg, sus_rxn_info, known_products, interm)
                continue

            # Add known reactions to the database
            is_updated = rxn_database.update(mssg, sus_rxn_info, known_products, inert_pairs, T)

            # Print messages if verbose
            if verbose:
                pairwise.inform_user(T, precursors, products, final_wts, mssg, sus_rxn_info, known_products, interm)
                rxn_database.print_info()

            # Check whether new reactions were found
            if is_updated:
                rxn_database.save()
                updated = True

        # Inform reaction database that precursors are changing
        rxn_database.make_global()

        # Remove the rxn that's just been tested
        del sorted_rxn_info[0]

        # If phase pure target obtained, halt campaign (unless --all specified)
        if (len(products) == 1) and (products[0] == target_product):
            if not all:
                print('Phase pure target obtained. Halting campaign.')
                print('Successful synthesis route: %s @ %s C' % (', '.join(precursors), int(T)))
                sys.exit()

        # If the pairwise rxn database was modified, update the precursor ranking accordingly
        if updated:
            sorted_rxn_info = update_ranking(rxn_database, sorted_rxn_info, explore)

    print('\nAll possible reactions sampled.')


