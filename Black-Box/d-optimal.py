from sklearn.preprocessing import MultiLabelBinarizer
import pandas as pd
import cvxpy as cp
import numpy as np
from itertools import product

# How many experiments to propose in D-optimal design
n_exps = 40

# Define precursor sets and temperatures
precursor_sets = [['Y2O3', 'BaO2', 'CuCO3'], ['Y2O3', 'BaO', 'CuCO3'], ['BaO2', 'CuCO3', 'Y2Cu2O5'], ['Y2(CO3)3', 'BaO2', 'CuCO3'], ['Y2O3', 'BaCO3', 'CuCO3'],
    ['Y2(CO3)3', 'BaO', 'CuCO3'], ['BaO', 'CuCO3', 'Y2Cu2O5'], ['Y2(CO3)3', 'BaCO3', 'CuCO3'], ['BaCO3', 'CuCO3', 'Y2Cu2O5'], ['Y2O3', 'BaO2', 'Cu2O'],
    ['Y2O3', 'CuCO3', 'BaCuO2'], ['BaO2', 'Cu2O', 'Y2Cu2O5'], ['Y2O3', 'BaO', 'BaO2', 'Cu2O'], ['Y2O3', 'BaO2', 'BaCO3', 'Cu2O'], ['Y2O3', 'BaO2', 'CuO'],
    ['BaO2', 'CuO', 'Y2Cu2O5'], ['BaO2', 'BaCO3', 'Cu2O', 'Y2Cu2O5'], ['BaO', 'BaO2', 'Cu2O', 'Y2Cu2O5'], ['Y2(CO3)3', 'CuCO3', 'BaCuO2'], ['Y2(CO3)3', 'BaO2', 'Cu2O'],
    ['Y2(CO3)3', 'BaO2', 'BaCO3', 'Cu2O'], ['Y2(CO3)3', 'BaO', 'BaO2', 'Cu2O'], ['Y2(CO3)3', 'BaO2', 'CuO'], ['Y2O3', 'BaCO3', 'CuO'], ['BaCO3', 'CuO', 'Y2Cu2O5'],
    ['Y2O3', 'BaO2', 'Cu2O', 'BaCuO2'], ['Y2(CO3)3', 'BaCO3', 'CuO'], ['Y2O3', 'BaO', 'CuO'], ['BaO', 'CuO', 'Y2Cu2O5'], ['Y2(CO3)3', 'BaO2', 'Cu2O', 'BaCuO2'],
    ['BaCO3', 'Cu2O', 'Ba2(CuO2)3', 'Y2Cu2O5'], ['Y2(CO3)3', 'BaO', 'CuO'], ['Y2O3', 'BaCO3', 'Cu2O', 'Ba2(CuO2)3'], ['BaO', 'Cu2O', 'Ba2(CuO2)3', 'Y2Cu2O5'],
    ['BaO2', 'Ba2(CuO2)3', 'Y2Cu2O5'], ['Y2(CO3)3', 'BaCO3', 'Cu2O', 'Ba2(CuO2)3'], ['BaCO3', 'Ba2(CuO2)3', 'Y2Cu2O5'], ['Y2O3', 'BaO', 'Cu2O', 'Ba2(CuO2)3'],
    ['Y2(CO3)3', 'BaO', 'Cu2O', 'Ba2(CuO2)3'], ['BaO', 'Ba2(CuO2)3', 'Y2Cu2O5'], ['Y2O3', 'CuO', 'BaCuO2'], ['BaCuO2', 'Y2Cu2O5'], ['Y2(CO3)3', 'CuO', 'BaCuO2'],
    ['Y2(CO3)3', 'Cu2O', 'BaCuO2', 'Ba2(CuO2)3'], ['Y2O3', 'Cu2O', 'BaCuO2', 'Ba2(CuO2)3'],
    ['Y2(CO3)3', 'Ba2(CuO2)3'], ['Y2O3', 'Ba2(CuO2)3']]
temperatures = np.array([600, 700, 800, 900]).reshape(-1, 1)

# Apply MultiLabelBinarizer
encoder = MultiLabelBinarizer()
precursors_encoded = encoder.fit_transform(precursor_sets)

# Create a list of tuples where each tuple is a combination of a precursor set and a temperature
design_space = list(product(range(len(precursor_sets)), temperatures.flatten()))

# Construct design matrix X
X = np.hstack((np.ones((len(design_space), 1)), precursors_encoded[np.array(design_space)[:, 0]], np.array([temp for _, temp in design_space]).reshape(-1, 1)))

# Scale the temperature axis from 0 to 1
X[:, -1] = (X[:, -1] - X[:, -1].min()) / (X[:, -1].max() - X[:, -1].min())

# Define the optimization problem
w = cp.Variable(len(X), nonneg=True)
objective = cp.Maximize(cp.trace(X.T @ cp.diag(w) @ X))
constraints = [cp.sum(w) == n_exps, 0 <= w, w <= 1]

# Solve the problem
problem = cp.Problem(objective, constraints)
problem.solve(solver=cp.CVXOPT)

if problem.status == 'optimal':
    # Round the solution and get the selected experiments
    selected_indices = np.where(np.round(w.value) > 0.5)[0]
    selected_experiments_encoded = pd.DataFrame(X[selected_indices, :], columns=['Intercept'] + list(encoder.classes_) + ['Temperature'])
    selected_experiments_real = pd.DataFrame([design_space[i] for i in selected_indices], columns=['Precursor Set Index', 'Temperature'])
    selected_experiments_real['Precursor Set'] = [precursor_sets[i] for i in selected_experiments_real['Precursor Set Index']]
    print(selected_experiments_real)
else:
    print('The problem does not have an optimal solution.')
