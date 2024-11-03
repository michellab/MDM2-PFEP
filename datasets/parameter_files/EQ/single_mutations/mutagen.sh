#!/bin/bash

# get the script directory and change to current directory
template_dir="$(dirname "$(readlink -f "$0")")"
work_dir="$(pwd)"

while getopts "h:l:" opt; do
  case $opt in
 l)
    echo "Will use $OPTARG λ-windows."
    # need to decrement the user input value by one because in scripts λ-indexing starts at 0
    # so if we want to run 11 λ-windows, our loop needs to run from values of 0 to 10 inclusive.
    lambda_windows=$((OPTARG - 1))
    ;;
  h)
    echo "DESCRIPTION

          Small convenience script to set up MDM2 protein FEP simulation sequence."
    echo "OPTIONS

      -l                         (./mutagen.sh -l)
          Number of λ-windows to use for the simulation. This option is required and is used to modify input files."
    exit 0
    ;;
  \?)
    echo "Invalid option: -$OPTARG" >&2
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument." >&2
    exit 1
    ;;
  esac
done

cd "$work_dir" || exit


for ((i = 0; i <= "$lambda_windows"; i++)); do
  mkdir "lambda_$i"

  ## em
  mkdir -p "lambda_$i"/em
  cp "$template_dir"/*.itp "lambda_$i"/
  cp "$template_dir"/mdm2_ions.top "$template_dir"/mdm2_ions.gro "$template_dir"/em.mdp "lambda_$i"/em || exit
  cd "lambda_$i"/em || exit
  sed -i "s/= LAMBDA_STATE/= $i/" em.mdp
  cd ..

  ## nvt
  mkdir nvt
  cp "$template_dir"/nvt.mdp nvt || exit
  cd nvt || exit
  sed -i "s/= LAMBDA_STATE/= $i/" nvt.mdp
  cd ..

  ## npt
  mkdir npt
  cp "$template_dir"/npt.mdp npt || exit
  cd npt || exit
  sed -i "s/= LAMBDA_STATE/= $i/" npt.mdp
  cd ..

  ## md prod for fep run
  mkdir md_fep
  cd md_fep || exit
  cp "$template_dir"/md_fep.mdp . || exit
  sed -i "s/= LAMBDA_STATE/= $i/" md_fep.mdp

  cd ..
  cp "$template_dir"/run_fep.sh .
  echo "Setup for λ-window $i completed."
  cd "$work_dir" || exit
done
