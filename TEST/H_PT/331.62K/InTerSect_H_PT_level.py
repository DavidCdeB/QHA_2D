#

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys
from sympy import *
import sympy as sym
import os
from itertools import chain
import pickle as pl


# Intial candidates for fit, per FU: - thus, the E vs V input data has to be per FU
E0_init = -941.510817926696  
V0_init = 63.54960592453 
B0_init = 76.3746233515232
B0_prime_init = 4.05340727164527


def BM(x, a, b, c, d):
         return  a + b*x + c*x**2 + d*x**3

def P(x, b, c, d):
     return 4.3597482E+3 * (-b - 2*c*x - 3 *d*x**2)

def H(x, a, b, c, d):
     return  a + b*x + c*x**2 + d*x**3


filefolder_Calcite_I_SG_167 = '../../Files_Outputs/Calcite_I/G_PT'
filefolder_Calcite_II_SG_14 = '../../Files_Outputs/Calcite_II/G_PT'
filefolder_Calcite_I_SG_167_EL_plus_E0 = '../../Files_Outputs/Calcite_I'
filefolder_Calcite_II_SG_14_EL_plus_E0 = '../../Files_Outputs/Calcite_II'
filefolder_Calcite_I_SG_167_ET = '../../Files_Outputs/Calcite_I'
filefolder_Calcite_II_SG_14_ET = '../../Files_Outputs/Calcite_II'


filefolder_energetics = 'F_vs_V_331.62K'

# Calcite I (Red triangles): 
V_C_I, E_C_I = np.loadtxt(os.path.join(filefolder_Calcite_I_SG_167, filefolder_energetics, 'F_vs_V_331.62K.dat'), skiprows = 1).T
V_C_I, E_plus_E0_C_I = np.loadtxt(os.path.join(filefolder_Calcite_I_SG_167_EL_plus_E0, 'EL_plus_E0_vs_V', 'EL_plus_E0_vs_V.dat'), skiprows = 1).T
V_C_I, ET_C_I = np.loadtxt(os.path.join(filefolder_Calcite_I_SG_167_ET, 'ET_vs_V', 'ET_vs_V_331.62K', 'ET_vs_V_331.62K.dat'), skiprows = 1).T

# 14 (Empty grey triangles):
V_14, E_14 = np.loadtxt(os.path.join(filefolder_Calcite_II_SG_14, filefolder_energetics, 'F_vs_V_331.62K.dat'), skiprows = 1).T
V_14, E_plus_E0_14 = np.loadtxt(os.path.join(filefolder_Calcite_II_SG_14_EL_plus_E0, 'EL_plus_E0_vs_V', 'EL_plus_E0_vs_V.dat'), skiprows = 1).T
V_14, ET_14 = np.loadtxt(os.path.join(filefolder_Calcite_II_SG_14_ET, 'ET_vs_V', 'ET_vs_V_331.62K', 'ET_vs_V_331.62K.dat'), skiprows = 1).T

init_vals = [E0_init, V0_init, B0_init, B0_prime_init]

popt_C_I, pcov_C_I = curve_fit(BM, V_C_I, E_C_I, p0=init_vals)
popt_14, pcov_14 = curve_fit(BM, V_14, E_14, p0=init_vals)

# Linspace for plotting the fitting curves:
V_C_I_lin = np.linspace(V_C_I[0], V_C_I[-1], 10000)
V_14_lin = np.linspace(V_14[0], V_14[-1], 10000)

fig_handle = plt.figure()

# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='black', label='Cubic fit Calcite I' )
p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='Cubic fit Calcite II')

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
fontP.set_size('15')

plt.legend((p1, p2, p5, p6), ("Calcite I", "Cubic fit Calcite I", "Calcite II", 'Cubic fit Calcite II'), prop=fontP)

global V0, B0, B0_prime
E0   =     popt_C_I[0] 
V0   =     popt_C_I[1]
B0   =     popt_C_I[2]
B0_prime = popt_C_I[3]

pressures_per_F_unit_C_I = P(V_C_I, V0, B0, B0_prime)
output_array_2 = np.vstack((E_C_I, V_C_I, pressures_per_F_unit_C_I)).T
np.savetxt('Volumes_and_pressures_C_I.dat', output_array_2, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

global V0_14, B0_14, B0_prime_14
E0_14   =     popt_14[0] 
V0_14   =     popt_14[1]
B0_14   =     popt_14[2]
B0_prime_14 = popt_14[3]


pressures_per_F_unit_14 = P(V_14, V0_14, B0_14, B0_prime_14)
output_array_2 = np.vstack((E_14, V_14, pressures_per_F_unit_14)).T
np.savetxt('Volumes_and_pressures_14.dat', output_array_2, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

# obtain the Temperature we are working at:
import subprocess
T = subprocess.check_output("basename `pwd`", shell=True)
print 'T = ', T

import string
new_T = string.replace(T, '_', ' ')
print 'new_T = ', new_T
T = T.rstrip()
new_T = new_T.rstrip()

plt.xlabel('$V$ / F.U. (Angstrom$^{3}$)', fontsize=20)
plt.ylabel(r'$(F = E + E_{ZP} + ET - TS)$ / F.U. (a.u.)', fontsize=20)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$ and (0.98 - 1.08)$V_{eq}$. T = %s" %new_T, fontsize=10)

plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot.pdf', bbox_inches='tight')
pl.dump(fig_handle,file('sinus.pickle_calcite_I_and_II_all_2_summary_better_plot','w'))

# Plotting P vs V:
fig_handle = plt.figure()

p2, = plt.plot(V_C_I_lin, P(V_C_I_lin, V0, B0, B0_prime), color='black', label='Cubic fit Calcite I' )
p6, = plt.plot(V_14_lin, P(V_14_lin, V0_14, B0_14, B0_prime_14), 'b', label='Cubic fit Calcite II')


# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, pressures_per_F_unit_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(V_14, pressures_per_F_unit_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
fontP.set_size('13')

plt.legend((p1, p2, p5, p6), ("Calcite I", "Cubic fit Calcite I", "Calcite II", 'Cubic fit Calcite II'), prop=fontP)

plt.xlabel('$V$ / F.U. (Angstrom$^{3}$)', fontsize=20)
plt.ylabel(r'$P = -\frac{\partial (E + E_{ZP} + ET - TS)}{\partial V}$ (GPa)', fontsize=20)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$ and (0.98 - 1.08)$V_{eq}$. T = %s" %new_T, fontsize=10)
plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot_P_vs_V.pdf', bbox_inches='tight')
pl.dump(fig_handle,file('sinus.pickle_calcite_I_and_II_all_2_summary_better_plot_P_vs_V','w'))

plt.figure()
H_C_I = E_plus_E0_C_I + ET_C_I + pressures_per_F_unit_C_I * V_C_I * (2.293710449E+17)*(1E-21) 
H_14 = E_plus_E0_14 + ET_14 + pressures_per_F_unit_14 * V_14 * (2.293710449E+17)*(1E-21)

popt_C_I, pcov_C_I = curve_fit(BM, V_C_I, H_C_I, p0=init_vals)
popt_14, pcov_14 = curve_fit(BM, V_14, H_14, p0=init_vals)

p1 = plt.scatter(V_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(V_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='black', label='Cubic fit Calcite I' )
p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='Cubic fit Calcite II')

fontP = FontProperties()
fontP.set_size('15')
plt.legend((p1, p2, p5, p6), ("Calcite I", "Cubic fit Calcite I", "Calcite II", 'Cubic fit Calcite II'), prop=fontP)

plt.xlabel('$V$ / F.U. (Angstrom$^{3}$)', fontsize=20)
plt.ylabel(r'$(H = E + E_{ZP} + \mathcal{E} + PV)$ / F.U. (a.u.)', fontsize=15)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$ and (0.98 - 1.08)$V_{eq}$. T = %s" %new_T, fontsize=10)
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_ALL_H.pdf', bbox_inches='tight')
pl.dump(fig_handle,file('sinus.pickle_calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_ALL_H','w'))



