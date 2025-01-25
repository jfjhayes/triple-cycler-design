# triple-cycler-design
Earth-Mars-Venus triple cycler trajectory optimisation using MATLAB's Genetic Algorithm and fmincon heuristic and constrained optimisers in support of ENG5041P: Individual Project 5 at the James Watt School of Engineering, University of Glasgow.

Run by downloading the ZIP file, extracting it, and placing it in your MATLAB workspace directory. To run, you will require MATLABs parallel computing and global optimisation toolboxes and a Lambert problem solver not included in this repository. This project consists of a genetic algorithm (GA) and a constrained optimisation phase, with several helper functions and plotting tools included.

It is a little clunky for swapping between outbound (Earth-Mars) and inbound (Mars-Earth) cyclers, but high-level changes can be made in defineParams.m and objectiveFunction.m. Currently, it is only configured for cyclers with a period equal to two EVM synodic periods (12.8 yrs) with four flybys. 
