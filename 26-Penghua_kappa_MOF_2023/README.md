Supporting Information for: "Sub-micrometer Phonon Mean Free Paths in Metal-Organic Frameworks Revealed by Machine-Learning Molecular Dynamics Simulations", Penghua Ying et al., In progress

This document contains the details for the input files used for lattice thermal conductivity calculations of MOF-5, ZIF-8, and HKUST-1 using molecular dynamics (MD) simulations methods implemented in the [GPUMD, version 3.7](https://github.com/brucefan1983/GPUMD) package. The complete input and output files for the NEP training and testing are freely accessible at this GitLab repository: [https://gitlab.com/brucefan1983/nep-data](https://gitlab.com/brucefan1983/nep-data).

**Repository Contents:**

- **HNEMD:** Files related to [Homogeneous Non-Equilibrium Molecular Dynamics](https://gpumd.org/theory/heat_transport.html#hnemd-method)
- **EMD:** Files related to [Equilibrium Molecular Dynamics](https://gpumd.org/theory/heat_transport.html#emd-method)
- **NEMD:** Files related to [Non-Equilibrium Molecular Dynamics](https://gpumd.org/theory/heat_transport.html#nemd-method)

To install the GPUMD package and run a simulation, follow the instructions in the [GPUMD documentation](https://gpumd.org/). After executing the MD simulations with the provided input files, use the following tutorials to post-process the output raw data:

- [HNEMD tutorial](https://github.com/brucefan1983/GPUMD/blob/master/examples/empirical_potentials/thermal_transport_nemd_and_hnemd/diffusive/Thermal%20Transport%20-%20Diffusive.ipynb)
- [NEMD tutorial](https://github.com/brucefan1983/GPUMD/blob/master/examples/empirical_potentials/thermal_transport_nemd_and_hnemd/ballistic/Thermal%20Transport%20-%20Ballistic.ipynb) 
- [EMD tutorial](https://github.com/brucefan1983/GPUMD/blob/master/examples/empirical_potentials/thermal_transport_emd/Thermal%20Transport%20-%20EMD.ipynb)