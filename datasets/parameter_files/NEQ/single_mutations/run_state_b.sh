#!/bin/bash
gmx_max_warn=1

root_dir=$(pwd)

for state in B; do
  cd state_"$state" || exit
  
  # energy minimization
  cd em || exit
  gmx grompp -f em.mdp -c mdm2_ions.gro -p mdm2_ions.top -o em.tpr -maxwarn "$gmx_max_warn"
  gmx mdrun -deffnm em || exit 1
  cd ..

  # nvt restrained
  cp em/em.gro em/mdm2_ions.top em/*.itp nvt/
  cd nvt || exit
  gmx grompp -f nvt.mdp -c em.gro -r em.gro -p mdm2_ions.top -o nvt.tpr -maxwarn "$gmx_max_warn"
  gmx mdrun -v -deffnm nvt || exit 1
  cd ..
  
  # npt restrained
  cp nvt/nvt.gro nvt/nvt.cpt nvt/mdm2_ions.top nvt/*.itp npt/
  cd npt || exit
  gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p mdm2_ions.top -o npt.tpr -maxwarn "$gmx_max_warn"
  gmx mdrun -v -deffnm npt || exit 1
  cd ..

  # production run
  cp npt/npt.gro npt/npt.cpt npt/mdm2_ions.top npt/*.itp md_prod/
  cd md_prod || exit
  gmx grompp -f md_prod.mdp -c npt.gro -t npt.cpt -p mdm2_ions.top -o md.tpr -maxwarn "$gmx_max_warn"
  gmx mdrun -v -deffnm md || exit 1

  # non-eq fep
  cd snapshots || exit

  # extract frames between 1 ns and final ns every 100 ps
  echo "System" | gmx trjconv -f ../md.trr -s ../md.tpr -o frame_.gro -sep -ur compact -pbc mol

  # build tpr files for each
  for file in ./*.gro; do
    gmx grompp -f ../non_eq_fep/ti.mdp -c "${file%.*}".gro -p ../mdm2_ions.top -o ../non_eq_fep/"${file%.*}".tpr -maxwarn  $(( "$gmx_max_warn" + 1 )) 
  done
  cd ..

  rm -r snapshots

  # run TI
  cd non_eq_fep || exit 1
  for file in ./*.tpr; do
    gmx mdrun -v -s "$file" -deffnm "${file%.*}" 
  done

  # return back to root_dir
  cd "$root_dir" || exit
done
