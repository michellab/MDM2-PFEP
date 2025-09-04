# MDM2-PFEP-Paper

Instructions for reproducing results from MDM2 alchemical side-chain mutation study.

## Installation instructions

### General installation

This paper utilizes the [BioSimSpace](https://biosimspace.openbiosim.org/) residue-of-interest (ROI) functionallity for creating perturbable protein molecules. This functionality is available in all BioSimSpace versions including and after `2024.2`. To install a general BioSimSpace installation:

```bash
conda create -n openbiosim -c conda-forge -c openbiosim biosimspace
conda activate openbiosim
```

### Exact paper version

The paper itself uses `2024.2` version of BioSimSpace for the perturbed molecule generation. If you wish to use this specific version, you have several options:

1. Running provided docker image (recommended)
2. Installing older version via conda

#### Running provided docker image

1. Install [Docker](https://www.docker.com/) into your machine
2. Pull the provided docker image: `docker pull akalpokas/biosimspace_2024.2:latest` (repository: https://hub.docker.com/r/akalpokas/biosimspace_2024.2)
3. Interactively run the pulled image: `docker run -it akalpokas/biosimspace_2024.2:latest`

You will be placed in the root directory and base conda environment. 

4. You can access BioSimSpace by activating the `openbiosim` conda environment: `conda activate openbiosim`
5. And finally run `ipython` and import BioSimSpace: `import BioSimSpace as BSS`

#### Installing older version via conda

> [!WARNING]
> There is no gurantee of this method working on every machine, **especially in the future** when some BioSimSpace dependencies might become unavailable due to version deprecation.

1. Install conda/miniconda/miniforge on your machine. For example miniforge is available here: https://github.com/conda-forge/miniforge
2. Install the specific BioSimSpace version: `conda create -n openbiosim -c conda-forge -c openbiosim "biosimspace==2024.2.0"`
3. Access BioSimSpace by activating the conda enviroment: `conda activate openbiosim`
4. And run `ipython` and import BioSimSpace: `import BioSimSpace as BSS`
___

## Input files

### Repository Structure

The input files in this repository is split into 3 different parent folders:
1. [raw_inputs](raw_inputs/) - documents how the raw protein and ligand input files were transformed into input simulation files.
2. [prepared_datasets](prepared_datasets/) - documents how the prepared input files were taken, assembled, and ran via GROMACS.
3. [analysis](analysis/) - documents how the post-processed simulation data was used to generate error statistics and plots.