# 03_analyse_free_energies

## Instructions

This paper utilized [cinnabar](https://github.com/OpenFreeEnergy/cinnabar), [alchemlyb](https://github.com/alchemistry/alchemlyb) and standard Python data analysis and plotting libraries in order to generate the accuracy and the precision statisics presented for the EQ and NEQ protocols. The `free_energy_analysis.ipynb` Jupyter notebook contains the code used in calculation of these metrics as well as plotting code for the post-processed EQ and NEQ free energy results which are located in the remaining subdirectories. Note that the Docker image provided does not contain these dependencies in order to keep the image size viable.