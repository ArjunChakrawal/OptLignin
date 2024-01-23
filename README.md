# OptLignin
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10553503.svg)](https://doi.org/10.5281/zenodo.10553503)
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/cloudposse.svg?style=social&label=Follow%20@ArjunChakrawal)](https://twitter.com/ArjunChakrawal)

This repository contains the data and scripts used in the publication "Modelling optimal ligninolytic activity during plant litter decomposition" by  Chakrawal et al., 2024.

In our study, we investigated optimal control problems related to ligninolytic activity during plant litter decomposition. We utilized a MATLAB-based tool developed by Viktor Leek at the Division of Vehicular Systems, Linköping University. For detailed installation guidelines and examples, please visit [Yoptimization](https://www.yoptimization.com/). To ensure the reproducibility of our work, we have archived the specific Yoptimization version used in our repository. It's important to note that another necessary component for Yoptimization is CasADi; specifically, we utilized the [casadi-windows-matlabR2016a-v3.5.5 version](https://github.com/casadi/casadi/releases/download/3.5.5/casadi-windows-matlabR2016a-v3.5.5.zip).

## Folder structure

- **Yop benchmark**: The Yop webpage has many worked out [examples](https://www.yoptimization.com/examples). We have included additional examples relevant for microbial ecology, in the Yop benchmark folder, demonstrating the capabilities of the toolbox.
	- Examples from Lenhart and Workman, Optimal control applied to biological models
	- Manzoni et al. 2022
	- Xuezhong Wang_Solving optimal control problems with MATLAB

	**OptLignin1.0**: The scripts used to create figures in the main text and supplementary materials are provided in this folder. The file names are self-explanatory. To perform model-data fitting for each study in the database, use the `run_all_ocp.m` file.
	- **data**: The digitized data is available as Excel files in the data folder.
	- **est_params**: The `est_params` folder contains Excel and text files documenting the estimated parameters for individual studies.
	- **fig**: The `fig` folder includes model-data fitting figures for all the litter bags from various studies.
	- **LSQ_fitting**: The `LSQ_fitting` folder contains scripts for model-data fitting, particularly for time-invariant rate parameters.
	- **results**: The `results` folder contains figures reported in the main text and supplementary materials.
