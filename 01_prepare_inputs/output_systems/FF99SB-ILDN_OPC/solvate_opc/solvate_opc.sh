echo "; opc water - externally added use only with opc.itp
EP_opc     0     0.0000     -1.358     A    1.78179743628e-01    0.000000e+00
HW_opc     1     1.0080      0.679     A    0.00000000000e+00    0.000000e+00
OW_opc     8    16.0000      0.000     A    3.16655208196e-01    8.903586e-01
Cl          17      35.45    0.0000  A   4.40104e-01  4.18400e-01
Na          11      22.99    0.0000  A   3.32840e-01  1.15897e-02

[ moleculetype ]
; molname       nrexcl
NA              1

[ atoms ]
; id    at type         res nr  residu name     at name  cg nr  charge
1       Na              1       NA              NA       1      1.00000


[ moleculetype ]
; molname       nrexcl
CL              1

[ atoms ]
; id    at type         res nr  residu name     at name  cg nr  charge
1       Cl              1       CL              CL       1      -1.00000
" > tmp_opc.atoms

# Now use awk to locate [ moleculetype ] in the file and insert the content of tmp_opc.atoms before it.
awk '/\[ moleculetype \]/ && !p {system("cat tmp_opc.atoms"); p=1} {print}' mdm2.top > mdm2_mod.top

# Remove the temporary file
rm tmp_opc.atoms

# Run the gromacs commands to solvate the system with OPC water
gmx editconf -f mdm2.gro -o box.gro -c -d 1.0 -bt dodecahedron
gmx solvate -cp box.gro -o mdm2_solv.gro -p mdm2_mod.top -cs opc.gro
mv mdm2_mod.top mdm2_solv.top
sed -i '/\[ system \]/i #include "opc.itp"' mdm2_solv.top
gmx grompp -f ions.mdp -c mdm2_solv.gro -p mdm2_solv.top -o ions.tpr -maxwarn 1
echo "OPC" | gmx genion -s ions.tpr -o mdm2_ions.gro -p mdm2_solv.top -pname NA -nname CL -neutral -conc 0.15
mv mdm2_solv.top mdm2_ions.top