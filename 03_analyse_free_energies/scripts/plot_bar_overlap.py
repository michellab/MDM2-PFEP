import glob
import matplotlib.pyplot as plt
import alchemlyb
from alchemlyb.visualisation import plot_mbar_overlap_matrix
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.estimators import BAR,MBAR
from copy import deepcopy

"""
The purpose of this script is to use alchemlyb MBAR estimator to plot an overlap matrix
for data that was generated using GROMACS BAR calculation. This is done by loading the reduced
potential data into a pandas dataframe (as per https://alchemlyb.readthedocs.io/en/latest/tutorial.html#parsing-the-free-energy-data).
The NaN values between non-adjacent λ-windows (e.g. 0.0 and 0.2 where the schedule is [0.0, 0.1, 0.2, ...]) are replaced with extremely
high values to mimick effectively 0 overlap between non-adjacent λ-windows. This then allows the dataframe to be fit using MBAR
estimator and allows the overlap matrix to be plotted. For robustness, the dG values calculated using BAR (with unmodified dataframe)
and MBAR (using modified dataframe) printed.

Instructions:
    run this script (python3 plot_bar_overlap.py) in a mutation directory containing all of the lambda subdirectories (lambda_0, lambda_1
    lambda_2 ...)
"""

# Extract the reduced potential data
lambda_dirs = glob.glob(f"lambda*/md_fep/md.xvg")
sorted_dirs = sorted(lambda_dirs, key=lambda x: int(x.split("/")[-3].split("_")[-1]))
data_list = [extract_u_nk(xvg, T=298) for xvg in sorted_dirs]
u_nk = alchemlyb.concat(data_list)

u_nk_bar = deepcopy(u_nk)

# Replace NaN values with extremely high values in u_nk dataframe
u_nk_mbar = u_nk.fillna(1000000)

# Fit the MBAR estimator to the data
mbar = MBAR()
mbar.fit(u_nk_mbar)

# Plot the overlap matrix
fig, ax = plt.subplots(figsize=(8, 8))
ax = plot_mbar_overlap_matrix(mbar.overlap_matrix, ax=ax)

# Save the overlap matrix
fig.savefig('mbar_overlap_matrix.png', dpi=300, bbox_inches='tight')

# Compare BAR and MBAR dG values
bar = BAR()
bar.fit(u_nk_bar)

print(f"BAR dG: {bar.delta_f_.loc[0.0, 1.0]} kT/mol")
print(f"MBAR dG: {mbar.delta_f_.loc[0.0, 1.0]} kT/mol")
