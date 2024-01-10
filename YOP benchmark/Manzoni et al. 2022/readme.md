# Microbial growth optimization for finite food source

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

<!-- ![Image](sche.svg) -->

<!-- <img src="sche.svg" width="400" height="150">  -->


This repository contains scripts to solve microbial growth as an optimal control problem and is a part of our publication [Manzoni et al. 2022](https://doi.org/10.3389/fevo.2023.1094269). 

We formulated an optimal control problem as the maximization of microbial growth rate with the dynamics of resource consumption as constraints (see figure). The microbial growth denoted by G is assumed to be a function of the uptake rate denoted by u of resources (e.g., carbon). Further, it is assumed that resources are lost due to other biotic/abiotic processes lumped into one term following first order kinetics and rate constant $\gamma$. The optimal control problem can be written as follows,


$\textrm{\textbf{Objective}}:  maximize \int_0^T G(u) dt$

$\textrm{subject to}  \ \ \dot{C}= -u - \gamma \ C + \mu \ G$, $C(0) =1$ and $C(T)=0$.

Growth rate is given by $G = e_{max} \ \beta \left(\frac{u-r}{u+\beta}\right)$ where $e_{max}$, $\beta$, and $r$ are maximum growth efficiency, half saturation constant, and maintenance rate. $\mu$ is the fraction of growth recycled as necromass in the resource pool. This formulates a free terminal time and fixed state endpoint optimal control problem with $C$ as the state and $u$ as the control variables. 

The above problem is solved using adapted forward-backward sweep from [Lenhart and Workman (2007)](https://doi.org/10.1201/9781420011418) and [*bvp4c*](https://se.mathworks.com/help/matlab/ref/bvp4c.html) function and [*YOP toolbox*](https://www.yoptimization.com/) from MATLAB. 


Contact

[![Twitter](https://img.shields.io/badge/Arjun_Chakrawal-%231DA1F2.svg?style=for-the-badge&logo=Twitter&logoColor=white)](https://twitter.com/ArjunChakrawal)
<!-- [![Twitter](https://img.shields.io/badge/ManzoniLab-%231DA1F2.svg?style=for-the-badge&logo=Twitter&logoColor=white)](https://twitter.com/ManzoniLab) -->

We welcome collaborations.
