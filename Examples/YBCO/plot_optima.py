import matplotlib.pyplot as plt
import numpy as np



# Track YBCO yield
num_rxns = 0.0
num_optim = []
with open('log') as f:
    lines = list(f.readlines())
    for i, line in enumerate(lines):
        if 'Products:' in line:
            if 'Ba2YCu3O7' in line:

                products = line.split()[1:]
                products = [entry.replace(',', '') for entry in products]
                ybco_ind = products.index('Ba2YCu3O7')
                amounts = lines[i+1].split()[1:]
                amounts = [entry.replace(',', '') for entry in amounts]
                curr_yield = 100*float(amounts[ybco_ind])

                if curr_yield == 100.0:
                    if len(num_optim) == 0:
                        num_optim.append(1.0)
                    else:
                        num_optim.append(num_optim[-1] + 1.0)
                else:
                    if len(num_optim) == 0:
                        num_optim.append(0.0)
                    else:
                        num_optim.append(num_optim[-1])

            else:
                if len(num_optim) == 0:
                    num_optim.append(0.0)
                else:
                    num_optim.append(num_optim[-1])

plt.figure()

plt.plot([0, 188], [10, 10], 'b-', linewidth=0.85, linestyle='dashed')

plt.plot([0, 188], [0, 10], 'r-', linewidth=0.85, linestyle='dashed', label='Random')

steps = list(range(len(num_optim)))
plt.plot(steps, num_optim, 'g-', label='ARROWS')

plt.xlabel('No. of Experiments', fontsize=21, labelpad=8)
plt.ylabel('Optima Found', fontsize=21, labelpad=16)

plt.xticks([0, 40, 80, 120, 160], fontsize=18)
plt.yticks(fontsize=18)

plt.xlim(0, 188)
plt.ylim(0, 11.35)

plt.legend(prop={'size': 19}, loc='lower right')

plt.text(10, 10.25, 'Total Optima', fontsize=19, color='blue')

plt.tight_layout()
plt.savefig('Log.png', dpi=400)
plt.close()

