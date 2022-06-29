from arrows import energetics, reactions, pairwise, exparser
from pymatgen.core.composition import Composition
from itertools import combinations
import numpy as np
import json
import time
import csv
import sys
import os


if __name__ == '__main__':

    verbose = False # Whether to print all info
    explore = True # Default to exploration
    all = False # Stop phase pure target obtained
    enforce_thermo = False # Allow rxns with dG > 0
    greedy = False # Assume low-T rxns always occur first
    reward_yield = True # Target partial yield
    for arg in sys.argv:
        if '--verbose' in arg:
            verbose = True
        if '--exploit' in arg:
            explore = False
        if '--all' in arg:
            all = True
        if '--all' in arg:
            all = True
        if '--enforce_thermo' in arg:
            enforce_thermo = True
        if '--greedy' in arg:
            greedy = True
        if '--pure' in arg:
            reward_yield = False

    # Load settings
    with open('Settings.json') as f:
        settings = json.load(f)
    available_precursors = settings['Precursors']
    if settings['Allow Oxidation'] == 'True':
        available_precursors.append('O2')
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

    assert 'Rxn_TD.csv' in os.listdir('.'), 'No reaction data found. Please run gather_rxns first.'

    # Load reaction data
    sorted_rxn_info = []
    with open('Rxn_TD.csv') as csv_file:
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

    # Exploration: prioritize no. of interfaces
    if explore:
        sorted_rxn_info = sorted(sorted_rxn_info, key=lambda x: (-x[-2], x[-1]))

    # Exploitation: prioritize maximal dG
    else:
        sorted_rxn_info = sorted(sorted_rxn_info, key=lambda x: (x[-1], -x[-2]))

    # Pairwise rxn database
    rxn_database = pairwise.rxn_database()

    # Possible temperature ordering
    increasing_temps = sorted(temps)
    decreasing_temps = increasing_temps.copy()
    decreasing_temps.reverse()

    # Keep track of highest TD driving force
    best_dG = 0.0

    # Iterate through each rxn
    probed_rxns = []
    num_rxns = len(sorted_rxn_info)
    for i in range(num_rxns):

        # Print info
        sys.stdout.flush()

        # Test rxn with highest dG
        rxn = sorted_rxn_info[0]

        # Check for updates
        updated = False

        # Iterate through each temperature, high to low
        for T in decreasing_temps:

            # Starting materials
            precursors, initial_amounts = rxn[0], rxn[1]

            # If no experiments have been performed yet, suggest first set
            if exp_data is None:
                print('-- Suggested experiment --')
                print('Precursors: %s' % precursors)
                print('Temperature: %s C' % T)
                sys.exit()

            # Parse experimental reaction data
            products, final_amounts = exparser.get_products(precursors, T, exp_data)

            # If precursors or temperature not sampled yet, suggest current experiment
            if products is None:
                if verbose:
                    print('Current Ranking:')
                    for rxn in sorted_rxn_info:
                        print(rxn)
                print('-- Suggested experiment --')
                print('Precursors: %s' % precursors)
                print('Temperature: %s C' % T)
                sys.exit()

            # If target product is made phase pure, finish run
            if len(products) == 1:
                sole_product = Composition(products[0]).reduced_formula
                if sole_product == Composition(target_product).reduced_formula:
                    if not all:
                        print('-- Optimal  synthesis route identified --')
                        print('Precursors: %s' % precursors)
                        print('Temperature: %s' % T)
                        sys.exit()

            # If no reaction occured, no need to analyze further
            if set(products) == set(precursors):
                continue

            # Perform reaction pathway analysis
            mssg, sus_rxn_info, known_products, interm, inert_pairs = pairwise.retroanalyze(precursors, initial_amounts, products, final_amounts,
                pd_dict, T, allowed_byproducts, open_sys, enforce_thermo, rxn_database, probed_rxns)

            # Add known reactions to the database
            is_updated = rxn_database.update(mssg, sus_rxn_info, known_products, inert_pairs, T)

            # Print messages if verbose
            if verbose:
                pairwise.inform_user(T, precursors, products, final_amounts, mssg, sus_rxn_info, known_products, interm)
                rxn_database.print_info()

            # Check whether new reactions were found
            if is_updated:
                updated = True

        # Re-do analysis from low to high T, now including known rxns
        for T in increasing_temps:

            # Starting materials
            precursors, initial_amounts = rxn[0], rxn[1]

            # Parse experimental reaction data
            products, final_amounts = exparser.get_products(precursors, T, exp_data)

            # Perform reaction pathway analysis
            mssg, sus_rxn_info, known_products, interm, inert_pairs = pairwise.retroanalyze(precursors, initial_amounts, products, final_amounts,
                pd_dict, T, allowed_byproducts, open_sys, enforce_thermo, rxn_database, probed_rxns)

            if mssg == 'Reaction already probed.':
                if verbose:
                    print(mssg)
                continue

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

            # Add known reactions to the database
            is_updated = rxn_database.update(mssg, sus_rxn_info, known_products, inert_pairs, T)

            # Print messages if verbose
            if verbose:
                pairwise.inform_user(T, precursors, products, final_amounts, mssg, sus_rxn_info, known_products, interm)
                rxn_database.print_info()

            # Check whether new reactions were found
            if is_updated:
                rxn_database.save()
                updated = True

        # Inform reaction database that precursors are changing
        rxn_database.make_global()

        # Remove rxn that's already been tested
        del sorted_rxn_info[0]

        # Check whether the reaction database was updated
        if updated:

            # If so, evolve all precursor sets using the latest information
            evolved_rxn_info = []
            for rxn_ind, starting_rxn in enumerate(sorted_rxn_info):
                original_set = starting_rxn[0]
                original_amounts = starting_rxn[1]
                starting_materials = starting_rxn[2]
                starting_amounts = starting_rxn[3]
                new_products = starting_rxn[4]
                new_materials, new_amounts = pairwise.pred_evolution(starting_materials, starting_amounts, rxn_database, greedy, min(temps))
                if set(new_materials) != set(starting_materials):
                    if verbose:
                        print('\nPredicted evolution: %s --> %s' % (starting_materials, new_materials))
                    if target_product in new_materials:
                        ind = new_materials.index(target_product)
                        expec_yield = new_amounts[ind]
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

            if reward_yield:

                # Exploration: prioritize no. of interfaces in evolved sets
                if explore:
                    sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (-x[-4], -x[-2], x[-1]))

                # Exploitation: prioritize expected target yield and maximal dG from evolved sets
                else:
                    sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (-x[-4], x[-1], -x[-2]))

            else:

                # Exploration: prioritize no. of interfaces in evolved sets
                if explore:
                    sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (-x[-2], x[-1]))

                # Exploitation: maximal dG from evolved sets
                else:
                    sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (x[-1], -x[-2]))

    print('All possible reactions sampled.')


