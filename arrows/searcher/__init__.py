from arrows import reactions
from itertools import combinations


def get_precursor_sets(available_precursors, target_products, allowed_byproducts=[], max_pc=None):

    # Ensure proper formatting
    if isinstance(target_products, str):
        target_products = [target_products]
    elif isinstance(target_products, list):
        target_products = target_products
    else:
        raise Exception("""Target products must be formatted as string or list""")

    if isinstance(allowed_byproducts, str):
        allowed_byproducts = [allowed_byproducts]
    elif isinstance(allowed_byproducts, list):
        allowed_byproducts = allowed_byproducts
    else:
        raise Exception("""Byproducts must be formatted as string or list""")

    # Get elements in chemical space
    elem_list = []
    for cmpd in available_precursors:
        elems = reactions.balancer.parseElems(cmpd)
        elem_list += elems
    elems = list(set(elem_list))

    if not max_pc:
        # Limit set by Gibbs phase rule
        max_pc = len(elems)

    balanced_sets = []
    for num_pc in range(2, max_pc + 1):
        possible_sets = combinations(available_precursors, num_pc)
        for pc_set in possible_sets:
            trial_soln = reactions.get_balanced_coeffs(pc_set, target_products)
            if not isinstance(trial_soln, str): # If reaction can be balanced
                balanced_sets.append([list(pc_set), target_products])
            else:
                for num_byp in range(1, len(allowed_byproducts) + 1):
                    possible_byproducts = combinations(allowed_byproducts, num_byp)
                    for byp_set in possible_byproducts:
                        all_products = target_products + list(byp_set)
                        trial_soln = reactions.get_balanced_coeffs(pc_set, all_products)
                        if not isinstance(trial_soln, str): # If reaction can be balanced
                            balanced_sets.append([list(pc_set), all_products])

    return balanced_sets

