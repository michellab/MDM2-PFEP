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

You will be placed in the `/home/ubuntu/` directory and in the `openbiosim` conda environment. 

5. And finally run `ipython` and import BioSimSpace: `import BioSimSpace as BSS`

#### Running provided Jupyter notebooks

If you wish to interact with the Jupyter notebooks obtained from the this GitHub repository, mount the local GitHub repository as a volume in the docker container:
```bash
docker run -p 10000:8888 -it -v /absolute/path/to/git/directory/:/home/ubuntu/mdm2 akalpokas/biosimspace_2024.2:latest
```

> [!NOTE]
> You can choose to change the port numbers on your machine if you have other services running on port `10000` for example.

Now once inside the container, change the directory to the mounted volume `mdm2` and launch the jupyter lab server:
```bash
jupyter-lab --no-browser --port=8888 --ip=0.0.0.0 --NotebookApp.token='' --NotebookApp.password=''
```

Now you should be able to access the Jupyter lab server outside the docker container on your host machine:
http://127.0.0.1:10000/lab
or
http://localhost:10000/lab


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
1. [01_prepare_inputs](01_prepare_inputs/) - documents how the raw protein and ligand input files were transformed into input simulation files.
2. [02_run_prepared_inputs](02_run_prepared_inputs/) - documents how the prepared input files were taken, assembled, and ran via GROMACS.
3. [03_analyse_free_energies](03_analyse_free_energies/) - documents how the post-processed simulation data was used to generate error statistics and plots.