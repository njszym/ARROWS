from arrows import energetics, reactions, pairwise, exparser
from pymatgen.core.composition import Composition
from itertools import combinations
import numpy as np
import json
import csv
import sys
import time


if __name__ == '__main__':

    # Load experimental rxn data
    with open('Exp_Data.json') as f:
        exp_data = json.load(f)
    exp_data = exp_data['Universal File']

    # User-defined exploration vs. exploitation
    explore = True # Default to exploration
    all = False # Stop phase pure target obtained
    for arg in sys.argv:
        if '--exploit' in arg:
            explore = False
        if '--all' in arg:
            all = True

    # Specify available precursors and desired product phase
    available_precursors = ['Y2O3', 'Y2C3O9', 'BaO', 'BaO2', 'BaCO3', 'Cu2O', 'CuO', 'CuCO3', 'BaCuO2', 'Ba2Cu3O6', 'Y2Cu2O5']
    target_product = 'Ba2YCu3O7' # or Y Ba2 Cu3 O6.5
    target_product = Composition(target_product).reduced_formula
    allowed_byproducts = [] # No byproducts allowed
    temps = [600, 700, 800, 900] # Temperatures to sample
    open_sys = True # Open system (air)
    enforce_thermo = False # Allows rxns with dG > 0

    # Build phase diagrams
    pd_dict = energetics.get_pd_dict(available_precursors, temps)

    # Load reaction data
    sorted_rxn_info = []
    with open('Simulated_Rxns.csv') as csv_file:
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

            # Parse experimental reaction data
            products, final_amounts = exparser.get_products(precursors, T, exp_data)

            # If target product is made phase pure, finish run
            if len(products) == 1:
                sole_product = Composition(products[0]).reduced_formula
                if sole_product == Composition(target_product).reduced_formula:
                    if not all:
                        print('Optimum synthesis route identified:')
                        print('Precursors: %s' % precursors)
                        print('Temperature: %s' % T)
                        sys.exit()

            # If no reaction occured, no need to analyze further
            if set(products) == set(precursors):
                continue

            # Perform reaction pathway analysis
            mssg, sus_rxn_info, known_products, interm, inert_pairs = pairwise.retroanalyze(precursors, initial_amounts, products, final_amounts,
                pd_dict, T, open_sys, enforce_thermo, rxn_database, probed_rxns)

            """
            # Print user messages
            pairwise.inform_user(T, precursors, products, final_amounts, mssg, sus_rxn_info, known_products, interm)
            """

            # Add known reactions to the database
            is_updated = rxn_database.update(mssg, sus_rxn_info, known_products, inert_pairs, T)
            """
            rxn_database.print_info()
            """

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
                pd_dict, T, open_sys, enforce_thermo, rxn_database, probed_rxns)

            if mssg == 'Reaction already probed.':
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

            # Print user messages
            pairwise.inform_user(T, precursors, products, final_amounts, mssg, sus_rxn_info, known_products, interm)

            # Add known reactions to the database
            is_updated = rxn_database.update(mssg, sus_rxn_info, known_products, inert_pairs, T)
            rxn_database.print_info()

            # Check whether new reactions were found
            if is_updated:
                updated = True

        # Inform reaction database that precursors are changing
        rxn_database.make_global()

        # Remove rxn that's already been tested
        del sorted_rxn_info[0]

        # Check whether the reaction database was updated
        if updated:

            # If so, evolve all precursor sets using the latest information
            evolved_rxn_info = []
            for starting_rxn in sorted_rxn_info:
                original_set = starting_rxn[0]
                original_amounts = starting_rxn[1]
                starting_materials = starting_rxn[2]
                starting_amounts = starting_rxn[3]
                new_products = starting_rxn[4]
                new_materials, new_amounts = pairwise.pred_evolution(starting_materials, starting_amounts, rxn_database)
                if set(new_materials) != set(starting_materials):
                    print('\nPredicted evolution: %s --> %s' % (starting_materials, new_materials))
                    if target_product in new_materials:
                        ind = new_materials.index(target_product)
                        expec_yield = new_amounts[ind]
                        new_products = [target_product]
                        energ = 0.0
                    else:
                        expec_yield = 0.0
                        new_products, energ = reactions.get_dG(new_materials, new_amounts, target_product, open_sys, pd_dict, T)
                    all_interfaces = [frozenset(pair) for pair in combinations(new_materials, 2)]
                    known_interfaces = list(rxn_database.as_dict().keys())
                    new_interfaces = set(all_interfaces) - set(known_interfaces)
                    new_interfaces = [set(interf) for interf in new_interfaces]
                    num_interfaces = len(new_interfaces)
                    evolved_rxn_info.append([original_set, original_amounts, new_materials, new_amounts, new_products, expec_yield, new_interfaces, num_interfaces, energ])
                else:
                    evolved_rxn_info.append(starting_rxn)

            # Exploration: prioritize no. of interfaces in evolved sets
            if explore:
                sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (-x[-4], -x[-2], x[-1]))

            # Exploitation: prioritize expected target yield and maximal dG from evolved sets
            else:
                sorted_rxn_info = sorted(evolved_rxn_info, key=lambda x: (-x[-4], x[-1], -x[-2]))


