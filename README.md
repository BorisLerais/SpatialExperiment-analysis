The repository contains the scripts used for data analysis of tissue imaging via multiplex immuno-fluorescence, adapted from a published pipeline (Windhager et al. 2023). 

The pipeline should be easily adapted to most tissue imaging technologies (IMC, MIBI, CoDEX, etc.)

This repository contains :
- A example jupyter notebook (Example_notebook.ipynb) which showcases some steps of the pipeline using a publicly available dataset.
- The list of all functions used, as a readable notebook (functions.ipynb) and as a sourceable r file (functions.r)
- For those that want to follow along the scripts, the packages needed to install are available in a r file version (install.r) or using conda to set up the environment (environment.ylm)
- Finally the groovy script used in qupath for gating a composite classifier is available in Sequential_gating.groovy