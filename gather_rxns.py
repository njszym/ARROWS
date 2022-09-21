from arrows import energetics, reactions, searcher
import json
import csv


if __name__ == '__main__':

    # Load settings
    with open('Settings.json') as f:
        settings = json.load(f)
    available_precursors = settings['Precursors']
    if settings['Allow Oxidation'] == 'True':
        available_precursors.append('O2')
    target = settings['Target']
    allowed_byproducts = settings['Allowed Byproducts']
    temps = settings['Temperatures']
    if 'Max Precursors' in settings.keys():
        max_pc = settings['Max Precursors']
    else:
        max_pc = None

    # Build phase diagrams
    pd_dict = energetics.get_pd_dict(available_precursors, temps)

    # Tabulate precursor sets that balance to produce target
    balanced_sets = searcher.get_precursor_sets(available_precursors, target, allowed_byproducts, max_pc)

    # Calculate reaction energies (at min T)
    rxn_info = []
    for (reactants, products) in balanced_sets:
        precursor_amounts = reactions.get_balanced_coeffs(reactants, products)[0]
        precursor_amounts = [round(val, 3) for val in precursor_amounts]
        rxn_energ = reactions.get_rxn_energy(reactants, products, min(temps), pd_dict[min(temps)])
        rxn_info.append([reactants, precursor_amounts, products, rxn_energ])

    # Sort rxns from most to least TD favorable
    sorted_info = sorted(rxn_info, key=lambda x: x[-1])

    # Save initial reaction data to csv file
    with open('Rxn_TD.csv', 'w+') as datafile:
        csv_writer = csv.writer(datafile)
        csv_writer.writerow(['Precursors', 'Amounts', 'Products', 'Reaction energy (meV/atom)'])
        for rxn in sorted_info:
            precursors = ' + '.join(rxn[0])
            amounts = ' + '.join([str(v) for v in rxn[1]])
            products = ' + '.join(rxn[2])
            energ = rxn[3]
            csv_writer.writerow([precursors, amounts, products, energ])
