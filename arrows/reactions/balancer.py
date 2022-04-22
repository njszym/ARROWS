from pymatgen.core.composition import Composition
import numpy as np


class RxnBalance(object):

    def __init__(self, reactants, products):

        self.reacs = reactants
        self.prods = products

    def balance(self):

        reac_elems = []
        for cmpd in self.reacs:
            reac_elems += parseElems(cmpd)
        reac_elems = sorted(list(set(reac_elems)))

        prod_elems = []
        for cmpd in self.prods:
            prod_elems += parseElems(cmpd)
        prod_elems = sorted(list(set(prod_elems)))

        # Normalize first product to 1
        reactants = list(self.reacs)
        norm_product = self.prods[0]
        if len(self.prods) > 1:
            for cmpd in self.prods[1:]:
                reactants.append(cmpd)

        # Elements in reactants and products must match
        if set(reac_elems) != set(prod_elems):
            return [np.array([0]*len(reactants)), 1000]

        elem_vec = []
        for cmpd in reactants:
            elem_vec.extend(parseElems(cmpd))
        elem_vec = list(set(elem_vec))

        reac_mat = []
        i = 0
        for cmpd in reactants:
            reac_mat.append([])
            reac_mat[i] = [0]*len(elem_vec)
            cmpd_comp = Composition(cmpd).as_dict()
            for elem in cmpd_comp.keys():
                j = 0
                for check_elem in elem_vec:
                    if elem == check_elem:
                        reac_mat[i][j] = cmpd_comp[elem]
                    j += 1
            i += 1

        rank = np.linalg.matrix_rank(reac_mat)
        if rank < len(reactants): # linearly dependent
            return [np.array([0]*len(reactants)), 1000]

        prod_vec = [0]*len(elem_vec)
        cmpd_comp = Composition(norm_product).as_dict()
        for elem in cmpd_comp.keys():
            j = 0
            for check_elem in elem_vec:
                if elem == check_elem:
                    prod_vec[j] = cmpd_comp[elem]
                j += 1

        A = np.array(reac_mat).transpose()
        b = np.array(prod_vec)

        return np.linalg.lstsq(A, b, rcond=None)[0], np.linalg.lstsq(A, b, rcond=None)[1]

def parseElems(formula):
    if '(' in formula:
        cmpd_name = ''
        rform = Composition(formula)
        for (e, n) in rform.items():
            cmpd_name += str(e) + str(int(n))
        formula = cmpd_name
    letters_only = ''.join([letter for letter in formula if letter.isalpha()])
    index = -1
    elems = []
    for letter in letters_only:
        if letter.isupper():
            elems.append(letter)
            index += 1
        else:
            elems[index] += letter
    return list(set(elems))


def main(reactants, products):

    num_reacs = len(reactants)

    if isinstance(products, str):
        products = [products]
    elif isinstance(products, list):
        products = products
    else:
        raise Exception("""Products must be formatted as string or list""")

    rxn_balancer = RxnBalance(reactants, products)
    soln = rxn_balancer.balance()

    reactant_coeffs = soln[0].flatten()[:num_reacs]
    product_coeffs = np.concatenate(([1.0], soln[0].flatten()[num_reacs:]))

    # Ensure all precursors participate
    if (reactant_coeffs > 1e-6).all():

        # Ensure all byproducts participate
        # Should be produced, not consumed
        if (product_coeffs[1:] < -1e-6).all():

            # Ensure equation is balanced
            tol = 1e-6 # Allow some tolerance
            if (len(soln[1]) == 0) or (soln[1] < tol):
                return [reactant_coeffs, abs(product_coeffs)]

            else:
                return 'Reaction cannot be balanced'

        else:
            return 'Not all byproducts are not formed'

    else:
        return 'Not all precursors participate in the reaction'

