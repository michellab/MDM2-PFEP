# MDM2-PFEP-Paper
Input files for MDM2 alchemical side-chain mutation study.

___

## Repository Structure
```
└── datasets
    ├── FF14SB_TIP3P
    │   └── input_structures
    │       ├── am_open
    │       │   ├── e23g
    │       │   ├── i19g
    │       │   ├── q18g
    │       │   ├── t16g
    │       │   ├── t16g_i19g
    │       │   └── v14g
    │       ├── apo_closed
    │       │   ├── i19g
    │       │   ├── q18g
    │       │   ├── t16g
    │       │   ├── t16g_i19g
    │       │   └── v14g
    │       └── nutlin_closed
    │           ├── e23g
    │           ├── i19g
    │           ├── q18g
    │           ├── t16g
    │           ├── t16g_i19g
    │           └── v14g
    ├── FF99SB-ILDN_OPC
    │   ├── input_structures
    │   │   ├── am_open
    │   │   │   ├── i19g
    │   │   │   ├── q18g
    │   │   │   ├── t16g
    │   │   │   ├── t16g_i19g
    │   │   │   └── v14g
    │   │   ├── apo_closed
    │   │   │   ├── i19g
    │   │   │   ├── q18g
    │   │   │   ├── t16g
    │   │   │   ├── t16g_i19g
    │   │   │   └── v14g
    │   │   └── nutlin_closed
    │   │       ├── i19g
    │   │       ├── q18g
    │   │       ├── t16g
    │   │       ├── t16g_i19g
    │   │       └── v14g
    │   └── solvent_parameters
    └── parameter_files
        ├── EQ
        │   ├── double_mutations
        │   └── single_mutations
        └── NEQ
            ├── double_mutations
            └── single_mutations
```

Each individual mutation folder ([e.g V14G mutation with FF14SB-TIP3P forcefield](datasets/FF14SB_TIP3P/input_structures/am_open/v14g)) contains input coordinate (`mdm2_ions.gro`), input topology (`mdm2_ions.top`) and any input equilibration positional restraints (`posre.itp` or `posre_lig.itp`) used to carry out that specific perturbation. These are fully solvated systems with ions added and are ready to be used in simulation without further preparation.

The [parameters folder](datasets/parameter_files/) contains the input simulation parameter files as well as scripts used to set up and carry out perturbation.

___

## Reproducing A Particular Simulation
In order to reproduce a particular perturbation with a given protocol and forcefield, the following steps can be taken:
1. Create a template folder for a given perturbation which contains both the system input and simulation paramter files.
2. Use the template folder to create replicate simulations for a particular perturbation.
3. Execute the scripts in the replicate folder to run the simulation.

### EQ - Single Mutation
For example, to set up a V14G EQ perturbation with FF14SB-TIP3P protocol:
```shell

cd datasets/FF14SB_TIP3P/input_structures/am_open/

# Create a template folder that will be used to store all of the input files required
cp -r v14g template_v14g

# Copy over the scripts and parameter files of interest
# NOTE: `nvt.mdp` and `npt.mdp` file "define" section might need to be uncommented/adjusted depending on whether you wish to run
# an apo or holo type of simulation. 
# If you wish to run an apo simulation, "define" section only needs "-DPOSRES" part.
# If you wish to run a holo simulation, "define" section needs "-DPOSRES -DPOSRES_LIG" parts to apply positional restraints to
# both ligand and the protein during the equilibration.
cp ../../../parameter_files/EQ/single_mutations/* template_v14g/

# The template directory is now fully prepared, now we can create replicates of the simulation without overwriting any files
mkdir v14g_repl_1
cd v14g_repl_1/

# -l 11 flag is used to tell the script to set up 11 λ-window calculation, this will create 11 folders numbered from 0 to 10
# and containing a script 'run_fep.sh' which can be used to reproduce the sequence in which the protocol was ran.
../template_v14g/mutagen.sh -l 11

# For example
cd lambda_0
./run_fep.sh
```

### EQ - Double Mutation
For example, to set up a T16G-I19G EQ perturbation with FF99SB*-ILDN-OPC protocol:
```shell

cd datasets/FF99SB-ILDN_OPC/input_structures/apo_closed/

# As before, we create a template folder
cp -r t16g_i19g template_t16g_i19g

# This time we copy over the files needed for the double mutations
cp ../../../parameter_files/EQ/double_mutations/* template_t16g_i19g/

mkdir t16g_i19g_repl_1
cd t16g_i19g_repl_1/

# This time we will pass -l 21 flag to the script in order to create a perturbation with 21 λ-windows. Same script as above can then
# be used in order to reproduce this mutation
../template_t16g_i19g/mutagen.sh -l 21
```

### NEQ - Single Mutation
For example to set up a I19G NEQ perturbation with FF14SB-TIP3P protocol:
```shell

cd datasets/FF14SB_TIP3P/input_structures/am_open/

# As before with EQ protocol, for the NEQ protocol we will also create a template folder
cp -r i19g template_i19g

# Copy over the scripts and parameter files of interest
cp ../../../parameter_files/NEQ/single_mutations/* template_i19g

# As before, now we can create our replicate folders
mkdir i19g_repl_1
cd i19g_repl_1/

# This time all we need to do is call the script from the template directory without passing it any parameters.
# This will create two folders 'state_A' and 'state_B' as well as the associated 'run_state_a.sh' and 'run_state_b.sh' scripts
# which can be used in order to reproduce this perturbation with NEQ protocol.
../template_i19g/mutagen.sh

# For example
./run_state_a.sh
```

### NEQ - Double Mutation
For example to set up a T16G-I19G NEQ perturbation with FF14SB-TIP3P protocol:
```shell

cd datasets/FF14SB_TIP3P/input_structures/nutlin_closed/

cp -r t16g_i19g template_t16g_i19g

# Copy over the scripts and parameter files of interest, these contain different TI run parameters compared to the single mutation input files
cp ../../../parameter_files/NEQ/double_mutations/* template_t16g_i19g

# As before, now we can create our replicate folders
mkdir t16g_i19g_repl_1
cd t16g_i19g_repl_1/

# Once we call the script below, we can run the simulations as before
../template_t16g_i19g/mutagen.sh

# For example
./run_state_a.sh
```