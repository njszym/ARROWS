from sklearn.ensemble import RandomForestRegressor
from pymatgen.core import Composition
from sklearn import preprocessing
from arrows import exparser
import numpy as np
import random
import json
import csv


"""
Copy experimental json into current dir.
"""
# Load experimental rxn data
with open('Exp.json') as f:
    exp_data = json.load(f)
exp_data = exp_data['Universal File']

# Load precursor sets
all_cmpds = []
precursor_sets = []
with open('Rxn_TD.csv') as csv_file:
    csv_reader = csv.reader(csv_file)
    for i, row in enumerate(csv_reader):
        if i != 0:
            reactants = row[0].split(' + ')
            reactants = [Composition(cmpd).reduced_formula for cmpd in reactants]
            precursor_sets.append(reactants)
            for cmpd in reactants:
                all_cmpds.append(cmpd)

# List of unique phases
all_cmpds = list(set(all_cmpds))
zeroes = [0.0]*len(all_cmpds)

# Assign each phase a unique onehot vector
onehot_dict = {}
for i, cmpd in enumerate(all_cmpds):
    onehot_vec = zeroes.copy()
    onehot_vec[i] = 1.0
    onehot_dict[cmpd] = np.array(onehot_vec)

all_inputs = []
for pset in precursor_sets:
    onehot_vec = np.array(zeroes.copy())
    for cmpd in pset:
        onehot_vec += onehot_dict[cmpd]
    all_inputs.append(list(onehot_vec))

# Convert input to output
def get_output(input):
    """
    Input: precursors & temperature
    Output: YBCO yield
    """

    precursor_set = []
    for i, val in enumerate(input):
        if val == 1.0:
            precursor_set.append(all_cmpds[i])

    temp = 900.0

    products, final_amounts = exparser.get_products(precursor_set, int(temp), exp_data)
    target_yield = 0.0
    for cmpd, amnt in zip(products, final_amounts):
        if cmpd == 'Ba2YCu3O7': # Target (YBCO)
            target_yield += amnt

    return target_yield

all_outputs = []
for input in all_inputs:
    all_outputs.append(get_output(input))

"""
Run many campaigns, each with
a different starting batch.
"""
all_logs = []
for iteration in range(100):

    # Define and shuffle input/outputs
    full_input = all_inputs.copy()
    full_output = all_outputs.copy()

    # Shuffle data
    merged_data = list(zip(full_input, full_output))
    random.shuffle(merged_data)
    full_input = [data[0] for data in merged_data]
    full_output = [data[1] for data in merged_data]

    # Start with 5 datapoints
    starting_index = 1
    known_input = full_input[:starting_index]
    known_output = full_output[:starting_index]
    leftover_input = full_input[starting_index:]
    leftover_output = full_output[starting_index:]

    # No. iterations to identify optima
    optima_indices = []

    # Iteratively add more datapoints as necessary
    for index in range(starting_index, len(full_input)):

        print(index)

        # Train
        model = RandomForestRegressor()
        model.fit(known_input, known_output)

        # Predict
        pred_output = model.predict(leftover_input)

        # Pure greedy selection
        acq_fn = pred_output
        max_index = np.argmax(acq_fn)
        select_input = leftover_input[max_index]
        select_output = leftover_output[max_index]

        print(select_input, select_output)

        # If actual maximum is found, halt optimization
        if select_output == 1.0:
            print('Optimum found')
            optima_indices.append(index)

        # Add to known input
        known_input = np.concatenate([known_input, [select_input]])

        # Add to known output
        known_output = np.append(known_output, [select_output])

        # Remove from leftover input
        new_input = []
        for i, vec in enumerate(leftover_input):
            if i != max_index:
                new_input.append(vec)
        leftover_input = new_input.copy()

        # Remove from leftover output
        new_output = []
        for i, val in enumerate(leftover_output):
            if i != max_index:
                new_output.append(val)
        leftover_output = new_output.copy()

        if len(leftover_input) == 0:
            break

    num = 0
    num_optim_found = []
    for i in range(188):
        if i in optima_indices:
            num += 1
        num_optim_found.append(num)

    all_logs.append(num_optim_found)

np.save('BO', np.array(all_logs))

