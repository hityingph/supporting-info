Supporting information for: “[Variable thermal transport in black, blue, and violet phosphorene from extensive atomistic simulations with a neuroevolution potential](https://www.sciencedirect.com/science/article/abs/pii/S0017931022011498)”, Penghua Ying et al., *Int. J. Heat Mass Tran.*, **2023**, 202, 123681. DOI: [10.1016/j.ijheatmasstransfer.2022.123681](https://doi.org/10.1016/j.ijheatmasstransfer.2022.123681)

The input scripts for HNEMD and NEMD simulations of monolayer phosphorene structures using [GPUMD](https://github.com/brucefan1983/GPUMD) package. Specifically, the [GPUMD 3.3.1 version](https://github.com/brucefan1983/GPUMD/releases/tag/v3.3.1) was used. Complete input and output files for the NEP training and testing are freely available at a Zenodo repository (http://dx.doi.org/10.5281/zenodo.6575727). 

**Contents of repository**

- HNEMD_diffusive: HNEMD simulations to calculate the diffusive thermal conductivity.
- NEMD_ballistic: NEMD simulations to calculate the ballistic thermal conductivity.
- Potential: The machine learning NEP model.

In each directory, the corresponding input files and output results for molecular dynamics simulations were provided. Additionaly, two python codes in root directory were also given to reproduce the Fig.5 and Fig.6. To run the python codes, the [thermo](https://github.com/hityingph/thermo) package should be loaded firstly.

