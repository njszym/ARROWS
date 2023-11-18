# ARROWS

![Alt text](Logo.png?raw=true)

A package designed to guide solid-state synthesis experiments toward maximal yield of a target phase. As described in the [corresponding publication](https://doi.org/10.1038/s41467-023-42329-9), this package operates under two key assumptions:
* Solid-state reactions tend to occur in pairs. In other words, a multi-phase mixture will generally react two phases at a time.
* The most  effective reaction pathway is the one that leads to maximal thermodynamic driving force (ΔG) at the target-forming step.

Based on these two assumptions, ARROWS prioritizes sets of precursors with large ΔG to form a user-specified target when starting a new experimental campaign. After the user has performed these experiments and fed their outcomes back to ARROWS, the package determines whether the attempted synthesis was successful, and if not, it learns which pairwise reactions formed detrimental intermediate phases that consumed the available free energy and therefore the target's formation. In subsequent experimental iterations, ARROWS proposes new sets of precursors that it expects to avoid such intermediates and therefore maintain a larger thermodynamic driving force to form the desired target.

The pairwise reactions learned by ARROWS are stored in a local file that can be transferred between various experimental campaigns, enabling the algorithm to become more efficient as this reaction database grows.

## Installation

First, make sure [git LFS](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage) is installed. This is necessary to download the json file containing temperature-dependent energetics from the Materials Project (MP).

Once LFS is installed, clone the repository:

```
git clone https://github.com/njszym/ARROWS.git
```

Then, to install all required modules, navigate to the cloned directory and execute:

```
pip install . --user
```

## Setting up a new campaign

To set up a new experimental campaign, a ```Settings.json``` file is required. An example is provided below:

```
{
    "Precursors" : ["Y2O3", "Y2C3O9", "BaO", "BaO2", "BaCO3", "Cu2O", "CuO", "CuCO3", "BaCuO2", "Ba2Cu3O6", "Y2Cu2O5"],
    "Target" : "Ba2YCu3O7",
    "Allowed Byproducts" : ["O2", "CO2"],
    "Temperatures" : [600, 700, 800, 900],
    "Open System" : "True",
    "Allow Oxidation" : "True"
}
```

Each flag should be customized to fit the desired experiments.

```Precursors```: A list of all available precursors that may be used throughout the synthesis experiments.

```Target```: The desired synthesis product.

```Allowed Byproducts```: Any phases that are allowed to form, in addition to the target phase. These often include gaseous phases (O2, CO2, NH3, H2O), but may also be some water-soluble products that are easy to remove post hoc (e.g., NaCl).

```Temperatures```: A list of temperatures that will be sampled during the synthesis experiments.

```Open System```: Set to True to account for the possibility of gaseous byproducts (O2 or CO2).

```Allow Oxidation```: Set to True if O2 may be included as a reactant. Otherwise, set to False if oxidation states may only be fixed or reduced (more common for high-temperature experiments under a reducing atmosphere). 

## Preparing possible precursor sets

Once the ```Settings.json``` file has been created, a list of possible precursor sets that may form the target phase can be automatically generated as follows:

```
python gather_rxns.py
```

This will generate a file called ```Rxn_TD.csv``` containing the information associated with each synthesis reaction. For example:

```
Precursors,Amounts,Products,Reaction energy (meV/atom)
Y2O3 + BaO2 + CuCO3,0.5 + 2.0 + 3.0,Y Ba2 Cu3 O6.5 + O2 + CO2,-684.5725652364363
Y2O3 + BaO + CuCO3,0.5 + 2.0 + 3.0,Y Ba2 Cu3 O6.5 + CO2,-667.0656866322632
BaO2 + CuCO3 + Y2Cu2O5,2.0 + 2.0 + 0.5,Y Ba2 Cu3 O6.5 + O2 + CO2,-594.9394718943933
Y2C3O9 + BaO2 + CuCO3,0.5 + 2.0 + 3.0,Y Ba2 Cu3 O6.5 + O2 + CO2,-591.8966565917668
...
```

Here, the first column reperesents the possible precursor sets. The second column represents the expected products (target phase + byproducts). The third column represents the DFT-calculated change in the free energy associated with that reaction (at the upper bound of the specified temperature range).


## Running experiments

Now that all possible precursor sets have been written to ```Rxn_TD.csv```, the first round of experiments can be suggested. This can be accomplished as follows:


```
python suggest.py
```

This will generate an output informing the user which precursors and temperature should be tested. For example:

```
-- Suggested experiment --
Precursors: ['Y2O3', 'BaO', 'BaO2', 'Cu2O']
Temperature: 900 C
```

After running the specified experiments, the results should be added to the ```Exp.json``` file. An example illustrating the format is given below:

```
{"Universal File":
    {
    "BaO2, CuCO3, Y2(CO3)3":
        {
        "Precursor stoichiometry": [1.0, 3.0, 1.0], 
        "Temperatures":
            {
            "600 C":
                {
                "products": ["Y2O3_199", "BaCO3_62", "CuO_15"],
                "product weight fractions": [30, 35.25, 34.75]
                },
            "700 C":
                {
                "products": ["BaCO3_62", "CuO_15", "Y2O3_199"],
                "product weight fractions": [49, 42.5, 8.5]
                }
            }
        }
    }
}
```

By running ```suggest.py``` again, ARROWS will learn from these results an output a new suggestion. This process can be repeated until (1) the target phase is formed without any impurities, or (2) all possible experiments have been exhausted.

## Pairwise reaction data

As ARROWS learns from the experimental results, it will build a database of pairwise reactions. This database will be saved in a file called ```PairwiseRxns.csv```, which will look something like:


```
Pairwise reactants,Pairwise Products,Temperature Range
Y2O3 + BaCuO2, BaY2CuO5, Reacts between 800-900 C
Y2Cu2O5 + Ba2(CuO2)3, O2 + CuO + Ba2YCu3O7, Reacts between 800-900 C
Ba2(CuO2)3 + BaCO3, O2 + CO2 + BaCuO2, Reacts between 600-800 C
Y2O3 + Ba2(CuO2)3, O2 + Ba2YCu3O7, Reacts between 600-700 C
Y2O3 + BaCO3, None, Does not react at or below 800 C
O2 + Cu2O, CuO, Reacts below 600 C
...
```

To utilize previously obtained pairwise reaction data for a new experimental campaign, place the corresponding ```PairwiseRxns.csv``` into the working directory before running ```suggest.py```. The reaction database for this campaign will include all of the known pairwise reactions from the csv file.

Caution: previous pairwise reaction data should only be used when working under the same atmosphere. Products may vary depending on the partial pressures of various gaseous species.

## Custom options

When suggesting new experiments, several options can be specified at run time: ```suggest.py --options```

Where ```--options``` includes the following:

```--verbose```: Print out detailed information associated with the pairwise reaction analysis.

```--explore```: This will prioritize precursor sets with the most new interfaces, as opposed to those with maximal driving force for target formation.

```--enforce_thermo```: If this option is specified, only thermodynamically favorable pairwise reactions will be considered.

```--greedy```: If a known pairwise reaction occurs below the minimum temperature, assume it *always* occurs first when the two reactants are present in other precursor sets.

```--partial_yield```: If specified, prioritize precursor sets that are predicted to form any yield of the target phase, even if not pure.

```--all```: Explore all possible synthesis routes, even after an optimal one has been identified.

```--batch_size=N```: Suggest a batch of (```N```) experiments at each iteration. By default, ARROWS runs sequentially by suggesting only one experiment at a time.
