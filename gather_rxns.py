from arrows import energetics, reactions, searcher
import csv


if __name__ == '__main__':

    # Specify available precursors, desired products, and allowed byproducts
    available_precursors = ['Y2O3', 'Y2C3O9', 'BaO', 'BaO2', 'BaCO3', 'Cu2O', 'CuO', 'CuCO3', 'BaCuO2', 'Ba2Cu3O6', 'Y2Cu2O5']
    target = 'Y Ba2 Cu3 O6.5'
    allowed_byproducts = ['O2', 'CO2']
    temps = [600, 700, 800, 900, 1000]

    # Build phase diagrams
    pd_dict = energetics.get_pd_dict(available_precursors, temps)

    # Tabulate precursor sets that balance to produce target
    balanced_sets = searcher.get_precursor_sets(available_precursors, target, allowed_byproducts)

    # Calculate reaction energies (at max T)
    rxn_info = []
    for (reactants, products) in balanced_sets:
        precursor_amounts = reactions.get_balanced_coeffs(reactants, products)[0]
        precursor_amounts = [round(val, 3) for val in precursor_amounts]
        rxn_energ = reactions.get_rxn_energy(reactants, products, max(temps), pd_dict[max(temps)])
        rxn_info.append([reactants, precursor_amounts, products, rxn_energ])

    # Sort rxns from most to least TD favorable
    sorted_info = sorted(rxn_info, key=lambda x: x[-1])

    # Save initial reaction data to csv file
    with open('Rxn_Data.csv', 'w+') as datafile:
        csv_writer = csv.writer(datafile)
        csv_writer.writerow(['Precursors', 'Amounts', 'Products', 'Reaction energy (meV/atom)'])
        for rxn in sorted_info:
            precursors = ' + '.join(rxn[0])
            amounts = ' + '.join([str(v) for v in rxn[1]])
            products = ' + '.join(rxn[2])
            energ = rxn[3]
            csv_writer.writerow([precursors, amounts, products, energ])
