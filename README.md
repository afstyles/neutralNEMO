# neutralNEMO (v0.9.4, early release)

Calculate neutral surfaces for NEMO ([Madec et al., 2019](https://zenodo.org/doi/10.5281/zenodo.1464816)) model runs using neutralNEMO!

This package makes it easier to use the highly-optimized algorithms of the [neutralocean](https://github.com/geoffstanley/neutralocean/tree/main) package ([Stanley et al., 2021](https://doi.org/10.1029/2020MS002436)) with model outputs from NEMO.

The package contains routines to correctly load temperature and salinity data alongside the necessary horizontal and vertical grid information. 

A steamlined routine is also included to calculate an approximately neutral density surface which is fixed to a specific location in time and space. Additional options to calculate the Veronis density and potential density surfaces are included.

# Installation

neutralNEMO can be installed using conda:

```
conda install -c conda-forge neutralnemo
```

or alternatively, neutralNEMO can be installed through pip

```
pip install neutralNEMO
```

# The future of neutralNEMO

This is currently an early release of neutralNEMO and I am keen to here how others get on with the package when they use it on their model runs.  If it works well for you, please spread the word and let me know as it will make my day! If you do encounter any problems, please contact me or raise an issue, I will try to help in any way I can.

# Documentation
Documentation and examples can be found at [https://neutralnemo.readthedocs.io](https://neutralnemo.readthedocs.io)
