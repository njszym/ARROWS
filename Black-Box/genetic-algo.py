from geneticalgorithm import geneticalgorithm as ga
from pymatgen.core.composition import Composition
from pymatgen.core import Structure
from itertools import combinations
from arrows import exparser
import numpy as np
import json
import csv
import sys
import os


open('log', 'w+')


"""
Copy experimental json into current dir.
"""
with open('Exp.json') as f:
    exp_data = json.load(f)
exp_data = exp_data['Universal File']

already_tried = []

def get_yield(X):
    """
    Function to be minimized.
    If the precursor set is balanced, YBCO yield is returned.
    Otherwise, 0.0 is returned.
    """

    # Temperature intervals of 100
    X[-1] = int(X[-1]*100)

    Y = [str(int(v)) for v in X]
    Y_str = ''.join(Y)
    if Y_str in already_tried:
        return 1.0
    else:
        already_tried.append(Y_str)

    sys.stdout.flush()

    # Get each item from input vector
    Ba2Cu3O6, BaCO3, BaCuO2, BaO, BaO2, Cu2O, CuCO3, CuO, Y2C3O9, Y2Cu2O5, Y2O3, temp = X

    # Ensure precursor set isn't empty
    if sum(X[:-1]) == 0:
        return 1.0

    # Convert one-hot input vector to a list of precursors
    precursor_set = set()

    Y2O3 = round(Y2O3, 0)
    if Y2O3 == 1.0:
        precursor_set.add('Y2O3')

    Y2C3O9 = round(Y2C3O9, 0)
    if Y2C3O9 == 1.0:
        precursor_set.add('Y2C3O9')

    BaO = round(BaO, 0)
    if BaO == 1.0:
        precursor_set.add('BaO')

    BaO2 = round(BaO2, 0)
    if BaO2 == 1.0:
        precursor_set.add('BaO2')

    BaCO3 = round(BaCO3, 0)
    if BaCO3 == 1.0:
        precursor_set.add('BaCO3')

    CuO = round(CuO, 0)
    if CuO == 1.0:
        precursor_set.add('CuO')

    Cu2O = round(Cu2O, 0)
    if Cu2O == 1.0:
        precursor_set.add('Cu2O')

    CuCO3 = round(CuCO3, 0)
    if CuCO3 == 1.0:
        precursor_set.add('CuCO3')

    BaCuO2 = round(BaCuO2, 0)
    if BaCuO2 == 1.0:
        precursor_set.add('BaCuO2')

    Ba2Cu3O6 = round(Ba2Cu3O6, 0)
    if Ba2Cu3O6 == 1.0:
        precursor_set.add('Ba2Cu3O6')

    Y2Cu2O5 = round(Y2Cu2O5, 0)
    if Y2Cu2O5 == 1.0:
        precursor_set.add('Y2Cu2O5')

    precursor_set = [Composition(cmpd).reduced_formula for cmpd in precursor_set]

    products, final_amounts = exparser.get_products(precursor_set, int(temp), exp_data)
    if products is not None:
        target_yield = 0.0
        for cmpd, amnt in zip(products, final_amounts):
            if cmpd == 'Ba2YCu3O7': # Target (YBCO)
                target_yield += amnt
        with open('log', 'a') as f:
            f.write('%s\n' % target_yield)
        return -target_yield
    else:
        return 1.0

if __name__ == '__main__':

    # Specify available precursors, desired products, and allowed byproducts
    available_precursors = ['Y2O3', 'Y2C3O9', 'BaO', 'BaO2', 'BaCO3', 'Cu2O', 'CuO', 'CuCO3', 'BaCuO2', 'Ba2Cu3O6', 'Y2Cu2O5']
    target_products = 'Y Ba2 Cu3 O6.5'
    allowed_byproducts = ['O2', 'CO2']

    num_boolean = len(available_precursors)
    num_vars = num_boolean + 1 # Precursors + temperature
    varbounds = np.array([[0, 1]]*num_boolean + [[6, 9]])
    vartypes = np.array([['int']]*num_boolean + [['int']])

    """
    May vary these hyperparameters.
    For the benchmark reported in our paper,
    we varied the population size between 5 and 50
    and averaged the results.
    """
    algorithm_params = {'max_num_iteration': 20000,\
        'population_size': 10,\
        'mutation_probability': 0.25,\
        'elit_ratio': 0,\
        'crossover_probability': 0.75,\
        'parents_portion': 0.3,\
        'crossover_type':'uniform',\
        'max_iteration_without_improv': None,\
        'function_timeout': 60.0}

    model = ga(function=get_yield, dimension=num_vars, variable_type_mixed=vartypes, variable_boundaries=varbounds, algorithm_parameters=algorithm_params)

    model.run()
