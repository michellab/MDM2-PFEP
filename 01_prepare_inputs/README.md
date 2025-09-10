# 01_prepare_inputs

This directory contains files and instructions for reproducing the parametrisation of the input files that were used in the simulations. Specifically it contains:
- [build_perturbable_gromacs_topologies.ipynb](build_perturbable_gromacs_topologies.ipynb) notebook showing the code used to generate perturbable protein systems with FF14SB and FF99SN-ILDN forcefields. **This notebook is intended to used with the Jupyter lab environment launched from the docker image as specified in the main repository directory.**
- [raw_proteins](raw_proteins/) wild-type and mutant MDM2 input structures used in generation of the GROMACS inputs
- [raw_ligands](raw_ligands/) input structures and parameters of the Nutlin-3a and AM-7209 ligands (parametrised using [acpype](https://github.com/alanwilter/acpype)) used in generation of the GROMACS inputs.

### FF99SB-ILDN_OPC example
Run the shown cells in the [build_perturbable_gromacs_topologies.ipynb](build_perturbable_gromacs_topologies.ipynb) notebook to generate the GROMACS inputs. Then enter the output directory for this system and rename the created files (`mdm2_am_open_v14g.gro` and `mdm2_am_open_v14g.top`) to `mdm2.gro` and `mdm2.top` in order to be compatible with scripts. Copy the files from `output_systems/FF99SB-ILDN_OPC/solvate_opc/` to this directory and run `solvate_opc.sh` to solvate the system. Now you should have `mdm2_ions.gro` and `mdm2_ions.top` files in your directory. Copy the `am_ligand.gro` file from `raw_ligands` and now you can generate the restraints used for the equilibration.

```bash
printf '1 \n q \n' | gmx make_ndx -f mdm2_ions.gro -o protein.ndx
printf '0 \n q \n' | gmx make_ndx -f am_ligand.gro -o ligand.ndx

echo "19" | gmx genrestr -f mdm2_ions.gro -n protein.ndx -o posre.itp
echo "0" | gmx genrestr -f mdm2_ions.gro -n ligand.ndx -o posre_lig.itp
```

Now add the following #ifdef statements to the `mdm2_ions.top` file to allow the use of restraints:

```
#ifdef POSRES
#include "posre.itp"
#endif
```

```
#ifdef POSRES_LIG
#include "posre_lig.itp"
#endif
```

### FF14SB_TIP3P example
As before, run the shown cells in the [build_perturbable_gromacs_topologies.ipynb](build_perturbable_gromacs_topologies.ipynb) notebook to generate the GROMACS inputs. Then rename the created files (`mdm2_apo_closed_t16g_i19g.gro` and `mdm2_apo_closed_t16g_i19g.top`) to `mdm2_ions.gro` and `mdm2_ions.top` (FF14SB simulations do not need to be post processed). Restraints for this system can be generated with:

```bash
printf '1 \n q \n' | gmx make_ndx -f mdm2_ions.gro -o protein.ndx
echo "21" | gmx genrestr -f mdm2_ions.gro -n protein.ndx -o posre.itp
```

And again the add the #ifdef statement to the `mdm2_ions.top` file:

```
#ifdef POSRES
#include "posre.itp"
#endif
```