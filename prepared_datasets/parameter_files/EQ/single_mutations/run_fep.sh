#!/bin/bash

# energy minimization
cp *.itp em/
cd em || exit
gmx grompp -f em.mdp -c mdm2_ions.gro -p mdm2_ions.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em || exit 1
cd ..

# nvt restrained
cp em/em.gro em/mdm2_ions.top em/*.itp nvt/
cd nvt || exit
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p mdm2_ions.top -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt || exit 1
cd ..

# npt restrained
cp nvt/nvt.gro nvt/nvt.cpt nvt/mdm2_ions.top nvt/*.itp npt/
cd npt || exit
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p mdm2_ions.top -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt || exit 1
cd ..

# production run
cp npt/npt.gro npt/npt.cpt npt/mdm2_ions.top npt/*.itp md_fep/
cd md_fep || exit
gmx grompp -f md_fep.mdp -c npt.gro -t npt.cpt -p mdm2_ions.top -o md.tpr -maxwarn 1
gmx mdrun -deffnm md || exit 1
