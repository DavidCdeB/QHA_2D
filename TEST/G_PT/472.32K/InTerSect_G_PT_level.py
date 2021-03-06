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

filefolder_energetics = 'F_vs_V_472.32K'

# Calcite I (Red triangles): 
V_C_I, E_C_I = np.loadtxt(os.path.join(filefolder_Calcite_I_SG_167, filefolder_energetics, 'F_vs_V_472.32K.dat'), skiprows = 1).T

# 14 (Empty grey triangles):
V_14, E_14 = np.loadtxt(os.path.join(filefolder_Calcite_II_SG_14, filefolder_energetics, 'F_vs_V_472.32K.dat'), skiprows = 1).T


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

#plt.figure()
fig_handle = plt.figure()

# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='black', label='Cubic fit Calcite I' )

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
fontP = FontProperties()
fontP.set_size('15')

plt.legend((p1, p2), ("Calcite I", "Cubic fit Calcite I"), prop=fontP)

plt.xlabel('$V$ / F.U. (Angstrom$^{3}$)', fontsize=20)
plt.ylabel(r'$(F = E + E_{ZP} + ET - TS)$ / F.U. (a.u.)', fontsize=20)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$ and (0.98 - 1.08)$V_{eq}$. T = %s" %new_T, fontsize=10)
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_I_summary_better_plot.pdf', bbox_inches='tight')
pl.dump(fig_handle,file('sinus.pickle_calcite_I_summary_better_plot','w'))


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


#000000000000000000

H_C_I = E_C_I + pressures_per_F_unit_C_I * V_C_I * (2.293710449E+17)*(1E-21) 
H_14 = E_14 + pressures_per_F_unit_14 * V_14 * (2.293710449E+17)*(1E-21)

output_array_3 = np.vstack((E_C_I, V_C_I, pressures_per_F_unit_C_I, H_C_I)).T
np.savetxt('E_V_P_H__C_I.dat', output_array_3, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 

output_array_4 = np.vstack((E_14, V_14, pressures_per_F_unit_14, H_14)).T
np.savetxt('E_V_P_H__14.dat', output_array_4, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 

# Saving into variables:

P_lin_C_I = P(V_C_I_lin, V0, B0, B0_prime)
H_lin_C_I = BM(V_C_I_lin, *popt_C_I) +  P(V_C_I_lin, V0, B0, B0_prime) * V_C_I_lin * (2.293710449E+17)*(1E-21)

P_lin_14 = P(V_14_lin, V0_14, B0_14, B0_prime_14)
H_lin_14 = BM(V_14_lin, *popt_14) + P(V_14_lin, V0_14, B0_14, B0_prime_14) * V_14_lin * (2.293710449E+17)*(1E-21) 


output_array_1 = np.vstack((P_lin_C_I, H_lin_C_I)).T
np.savetxt('P_lin_C_I__H_lin_C_I.dat', output_array_1, header="P(GPa) \t   H per F unit (a.u)", fmt="%0.13f")

output_array_2 = np.vstack((P_lin_14, H_lin_14)).T
np.savetxt('P_lin_14__H_lin_14.dat', output_array_2, header="P(GPa) \t    H per F unit (a.u)", fmt="%0.13f")

init_vals = [E0_init, V0_init, B0_init, B0_prime_init]

popt_HofP_C_I, pcov_HofP_C_I = curve_fit(H, pressures_per_F_unit_C_I, H_C_I, p0=init_vals)
popt_HofP_14, pcov_HofP_14 = curve_fit(H, pressures_per_F_unit_14, H_14, p0=init_vals)


pressures_per_F_unit_14_sorted = np.sort(pressures_per_F_unit_14)
pressures_per_F_unit_14_lin = pressures_per_F_unit_14_sorted

pressures_per_F_unit_C_I_sorted = np.sort(pressures_per_F_unit_C_I)
pressures_per_F_unit_C_I_lin = pressures_per_F_unit_C_I_sorted

P_lin_for_coll = np.linspace(0.5, 16, 10000)#  100000)


# Evaluating:
H_C_I_lin_for_coll = H(P_lin_for_coll, *popt_HofP_C_I)
H_14_lin_for_coll = H(P_lin_for_coll, *popt_HofP_14)

fig_handle = plt.figure()

p1, = plt.plot(P_lin_for_coll, H(P_lin_for_coll, *popt_HofP_C_I), color='black', label='Calcite I cubic fit' )
p2, = plt.plot(P_lin_for_coll, H(P_lin_for_coll, *popt_HofP_14), color='blue', label='Calcite II cubic fit' )

fontP = FontProperties()
fontP.set_size('15')

plt.legend((p1, p2), ("Calcite I cubic fit", "Calcite II cubic Fit"), prop=fontP)

plt.xlabel('$P$ / F.U. (GPa)', fontsize=20)
plt.ylabel(r'$(G = E + E_{ZP} + ET - TS + PV)$ / F.U. (a.u.)', fontsize=15)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$ and (0.98 - 1.08)$V_{eq}$. T = %s" %new_T, fontsize=10)
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_ALL_H.pdf', bbox_inches='tight')
pl.dump(fig_handle,file('sinus.pickle_calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_ALL_H','w'))


# obtain the Temperature we are working at:
import subprocess
T = subprocess.check_output("basename `pwd`", shell=True)

import string
new_T = string.replace(T, '_', ' ')
T = T.rstrip()
new_T = new_T.rstrip()


T_folder = string.replace(new_T, ' ', '')
T_folder = string.replace(new_T, 'K', '')

T_folder_float = float(T_folder)


print 'Performing the Analytic intersection....'

fig = plt.figure()

#  Obtaining a cubic expression for H(P):

# Reminder:
#H_C_I = E_C_I + pressures_per_F_unit_C_I * V_C_I * (2.293710449E+17)*(1E-21) 
#H_14 = E_14 + pressures_per_F_unit_14 * V_14 * (2.293710449E+17)*(1E-21)

init_vals = [E0_init, V0_init, B0_init, B0_prime_init]

popt_HofP_C_I, pcov_HofP_C_I = curve_fit(H, pressures_per_F_unit_C_I, H_C_I, p0=init_vals)
popt_HofP_14, pcov_HofP_14 = curve_fit(H, pressures_per_F_unit_14, H_14, p0=init_vals)


pressures_per_F_unit_14_sorted = np.sort(pressures_per_F_unit_14)
pressures_per_F_unit_14_lin = pressures_per_F_unit_14_sorted


pressures_per_F_unit_C_I_sorted = np.sort(pressures_per_F_unit_C_I)
pressures_per_F_unit_C_I_lin = pressures_per_F_unit_C_I_sorted


# Linspace for plotting the fitting curves:
P_C_I_lin = np.linspace(pressures_per_F_unit_C_I_lin[0], pressures_per_F_unit_C_I_lin[-1], 10000)
P_14_lin = np.linspace(pressures_per_F_unit_14_lin[0], pressures_per_F_unit_14_lin[-1], 10000)

fig_handle = plt.figure()

# Plotting the fitting curves:
p2, = plt.plot(P_C_I_lin, H(P_C_I_lin, *popt_HofP_C_I), color='black', label='Cubic fit Calcite I' )
p6, = plt.plot(P_14_lin, H(P_14_lin, *popt_HofP_14), 'b', label='Cubic fit Calcite II')

# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)


fontP = FontProperties()
fontP.set_size('13')

plt.legend((p1, p2, p5, p6), ("Calcite I", "Cubic fit Calcite I", "Calcite II", 'Cubic fit Calcite II'), prop=fontP)

global a0, a1, a2, a3
a0     =     popt_HofP_C_I[0] 
a1     =     popt_HofP_C_I[1]
a2     =     popt_HofP_C_I[2]
a3     =     popt_HofP_C_I[3]

global a0_s2, a1_s2, a2_s2, a3_s2
a0_s2        =     popt_HofP_14[0] 
a1_s2        =     popt_HofP_14[1]
a2_s2        =     popt_HofP_14[2]
a3_s2        =     popt_HofP_14[3]

print 'a0 = ', a0
print 'a1 = ', a1
print 'a2 = ', a2
print 'a3 = ', a3

print 'a0_s2 = ', a0_s2
print 'a1_s2 = ', a1_s2
print 'a2_s2 = ', a2_s2
print 'a3_s2 = ', a3_s2


print """ 
The equations are the following:
G_I (P) = a0 + a1*P + a2*P**2 + a3*P**3 
G_II (P) = a0_s2 + a1_s2*P + a2_s2*P**2 + a3_s2*P**3 
"""
print('G_I (P) = ({a0}) + ({a1})*P + ({a2})*P**2 + ({a3})*P**3 '.format(a0 = a0, a1 = a1, a2 = a2, a3 = a3, ))

print """
"""
print('G_II (P) = ({a0_s2}) + ({a1_s2})*P + ({a2_s2})*P**2 + ({a3_s2})*P**3 '.format(a0_s2 = a0_s2, a1_s2 = a1_s2, a2_s2 = a2_s2, a3_s2 = a3_s2))

print """
"""

print """
G_I (P) = G_II (P)
"""

# Setting "P" to be symbolic:
P = sym.symbols('P') #, real=True)

def z_I(P):
        return   a0 + a1*P + a2*P**2 + a3*P**3 

def z_II(P):
        return   a0_s2 + a1_s2*P + a2_s2*P**2 + a3_s2*P**3 

# Crude intersection:
sol = sym.solve(z_I(P) - z_II(P) , P)
print 'sol_ H_I(P) - H_II(P)  =', sol

# Transform to complex notation, in order to
# better discard the complex root afterwards.
# Use of evalf to obtain better precision:

def is_close(a,b,tol):
    if abs(a-b)<tol: return True
    else: return False

P = sym.symbols('P')
sol = [complex(x.evalf()) for x in sol]

real_solutions = []
for x in sol:
    print 'x = ', x
    print 'x.real = ', x.real
    print 'abs(x - x.real) = ', abs(x - x.real)

    if is_close(x,x.real,10**(-10)): real_solutions.append(x.real)

print 'real_solutions = ', real_solutions

real_roots = [x for x in real_solutions if 0.1 <= x <= 4] or real_solutions
print 'real_roots = ', real_roots

# Transform each element of the list to a numpy array:
real_roots_zero_to_four = np.array(real_roots)

# Let's grab the root located between 0.1GPa and 4GPa (true for CI-CII phase trans.)
# If no element in the [real_roots] list is between 0.1 and 4.0, then the following line will also return the [real_roots] list.
#real_roots_zero_to_four = [x for x in real_roots if 0.1 <= x <= 4.0] or real_roots
#print 'real_roots_zero_to_four = ', real_roots_zero_to_four

P_real_intersection = real_roots_zero_to_four[0]
H_real_intersection = z_I(real_roots_zero_to_four[0])

output_array_2 = np.vstack((T_folder_float, P_real_intersection, H_real_intersection)).T
np.savetxt('P_H_analytic_intersection_T_%sK.dat' %T_folder , output_array_2, header="Temperature (K) \t Pressure_Intersection (GPa) \t H_Intersection = E + PV (a.u.)", fmt="%0.13f")

plt.xlabel(r'$P$ (GPa)', fontsize=20)
plt.ylabel(r'$(G = E + E_{ZP} + ET - TS + PV)$ / F.U. (a.u.)', fontsize=15)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 0.98)$V_{eq}$ and (0.98 - 1.08)$V_{eq}$ T = %s" %new_T, fontsize=10)
plt.ticklabel_format(useOffset=False)
ax = fig_handle.add_subplot(111)
ax.annotate('Analytic\nIntersection\nP= %g GPa\nG = %g a.u.' %(P_real_intersection, H_real_intersection), xy=(P_real_intersection, H_real_intersection), xytext=(P_real_intersection+2.5767, H_real_intersection-0.05), fontsize=15,
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='purple'),
            )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_with_analytic_intersection.pdf', bbox_inches='tight')
pl.dump(fig_handle,file('sinus.pickle_calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_with_analytic_intersection','w'))

#

#plt.show()


