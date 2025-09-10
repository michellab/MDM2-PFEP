#!/bin/bash
set -e

# get the script directory and change to current directory
template_dir="$(dirname "$(readlink -f "$0")")"

cp "$template_dir"/mdm2_ions.gro "$template_dir"/mdm2_ions.top "$template_dir"/*.itp .

for lambda in {0..1}; do
  if [[ $lambda -eq 0 ]]; then
    state="A"
    lambda_delta=8e-6
  elif [[ $lambda -eq 1 ]]; then
    state="B"
    lambda_delta=-8e-6
  fi
  
  mkdir state_"$state"
  cd state_"$state"
  mkdir em nvt npt md_prod
  mkdir md_prod/snapshots md_prod/non_eq_fep

  # em
  cp "$template_dir"/em.mdp em/
  cp ../mdm2_ions.gro ../mdm2_ions.top ../*.itp em/
  sed -i "s/= LAMBDA_STATE/= $lambda/" em/em.mdp

  # nvt
  cp "$template_dir"/nvt.mdp nvt/
  sed -i "s/= LAMBDA_STATE/= $lambda/" nvt/nvt.mdp

  # npt
  cp "$template_dir"/npt.mdp npt/
  sed -i "s/= LAMBDA_STATE/= $lambda/" npt/npt.mdp

  # md_prod
  cp "$template_dir"/md_prod.mdp md_prod/
  sed -i "s/= LAMBDA_STATE/= $lambda/" md_prod/md_prod.mdp

  # non-eq fep
  cp "$template_dir"/ti.mdp md_prod/non_eq_fep 
  sed -i "s/= LAMBDA_STATE/= $lambda/" md_prod/non_eq_fep/ti.mdp
  sed -i "s/= LAMBDA_DELTA/= $lambda_delta/" md_prod/non_eq_fep/ti.mdp
  cd ..
done

cp "$template_dir"/run_state_a.sh "$template_dir"/run_state_b.sh .
