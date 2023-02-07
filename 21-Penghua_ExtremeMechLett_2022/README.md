Supporting information for: “[Atomistic insights into the mechanical anisotropy and fragility of monolayer fullerene networks using quantum mechanical calculations and machine-learning molecular dynamics simulations](https://www.sciencedirect.com/science/article/abs/pii/S235243162200205X)”, Penghua Ying et al., *Extreme Mech. Lett.*, **2023**, 58, 101929. DOI: [10.1016/j.eml.2022.101929](https://doi.org/10.1016/j.eml.2022.101929)

The input scripts for tensile simulations of monolayer qHPF structures using molecular dynamics in [GPUMD](https://github.com/brucefan1983/GPUMD) package. Specifically, the [GPUMD 3.3.1 version](https://github.com/brucefan1983/GPUMD/releases/tag/v3.3.1) was used. Complete input and output files for the NEP training and testing are freely available at a Zenodo repository (http://dx.doi.org/10.5281/zenodo.7018573). 

**Contents of repository**

- Size: Size effect on tensile behaviour of monolayer qHPF.
- Strain_rate: Strain rate effect on tensile behaviour of monolayer qHPF.
- Temperature: Temperature effect on tensile behaviour of monolayer qHPF.
- Potential: The machine learning NEP model.

In each directory, the corresponding input files and output results for molecular dynamics simulations were provided. Additionaly, the python codes were also given to reproduce the Fig.7, Fig.8, and Fig.9. To run the python codes, the [thermo](https://github.com/hityingph/thermo) package should be loaded firstly.

