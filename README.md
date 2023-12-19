# Micrococcal Nuclease Simulation

This repository contains a Python-based simulation that models the action of Micrococcal Nuclease (MNase) on DNA 11nm fibres. The core components of the simulation include:

- **Fibre Class:** Represents the 11nm DNA fibre structure, including nucleosomes and linker DNA segments.

- **MNase Class:** Simulates the stochastic cleavage action of MNase on DNA fibres, randomly selecting nucleosomes and linkers for enzymatic digestion.

- **Main Simulation Logic:** Orchestrates the simulation over a series of time steps, continuously updating the state of the DNA fibres based on MNase action.

## Key Features:

- **Stochastic Modelling:** Utilizes Monte Carlo methods to realistically simulate the random nature of enzymatic cleavage on DNA.

- **Object-Oriented Approach:** Employs classes to effectively represent complex biological entities and their interactions.

- **Statistical Analysis:** The output, consisting of the lengths of DNA fibres post-MNase treatment, is suitable for further statistical analysis to understand the distribution and average effects of the nuclease.

## Objectives

* [x] Initialize git repository
* [x] Create boilerplate files
  * [x] README.md
  * [x] .gitignore
  * [x] main.py
  * [x] model/\_\_init\_\_.py
* [x] Define classes
  * [x] Fibre (data structure)
    * [x] Define attributes: n_nuclesomes, left_linker, right_linker
    * [x] Implement constructor (\_\_init\_\_ method)
    * [x] Implement method to determine number of nucleotides in Fibre
  * [x] MNase (namespace)
    * [x] Implement method to select a random nucleosome index
    * [x] Implement method to select either the left or right linker to cleave
    * [x] Implement method to cleave the fibre
* [ ] Implement Visualization
  * [ ] Install matplotlib
  * [ ] Define a function to visualize distribution of fibre lengths as a histogram
  * [ ] Define a function to visualize fibre length over time
