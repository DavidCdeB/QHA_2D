# 
# This script extracts ET, TS EL E0 from many outputs and prints out:
# V vs F @ each Temperature.
# These files are needed for the last step of the QHA methodology,
# where we need tables of F vs V 
# in order to then compute P=dF/dV and then G(P,T)

import re
import os
import glob
from itertools import islice
import numpy as np
import sys
import shutil
import subprocess
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

os.system("rm -Rf F_vs_V_*")
os.system("rm -Rf F_Anh_vs_V_*")
os.system("rm -Rf F_tot_vs_V_*")

n_volume = []
path='./'
template = os.path.join(path, '*.out')

# Setting the number of formula units as a raw_input:
n_F_u = raw_input("""
Please type as an integer the number of formula units in the primitive cell. 
For example, Calcite I contains 2 formula units in the primitive (rombohedral) cell and 6 formula units in the crystallographic (hexagonal) cell. Thus, the number to be introduced is:   2 <and press ENTER>
""")

n_F_u = float(n_F_u)
n_F_u = int(float(n_F_u))

# Setting the number of atoms per cell as a raw_input:
n_atoms_per_cell = raw_input("""
Please type as an integer the number of atoms per cell in the primitive cell. 
For example, Calcite I contains 2 formula units in the primitive (rombohedral) cell and 6 formula units in the crystallographic (hexagonal) cell. The number of atoms in the primitive cell is 10. Thus, the number to be introduced is:   10 <and press ENTER>
""")

n_atoms_per_cell = float(n_atoms_per_cell)
n_atoms_per_cell = int(float(n_atoms_per_cell))


# Extracting each thermodynamic variable:
ET = []
TS = []
EL = []
E0 = []
VOLUME_EACH = []
T = []
ALL_FREQ = []

for fname in glob.glob(template):
  print fname
  f = open(fname, 'r')
  real_part = False
  n_freqs = []
  Temps = []
  for line in f:

        if re.match(r"^ ET            :", line):
         start = line.find(':') + 8
         end = line.find(':') + 22
         result_ET = line[start:end]
         ET.append(result_ET)

        if re.match(r"^ TS            :", line):
         start = line.find(':') + 8
         end = line.find(':') + 22
         result_TS = line[start:end]
         TS.append(result_TS)

        if re.match(r"^ EL            :", line):
         start = line.find(':') + 4
         end = line.find(':') + 22
         result_EL = line[start:end]
         EL.append(result_EL)

        if re.match(r"^ E0            :", line):
         start = line.find(':') + 8
         end = line.find(':') + 22
         result_E0 = line[start:end]
         E0.append(result_E0)

        if re.match(r"^ AT \(T =", line):
         start = line.find('T =') + 4
         end = line.find('K')
         result_Temperatures = line[start:end]
         T.append(result_Temperatures)
         
         Temps.append(result_Temperatures)

        if 'LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - PRIMITIVE CELL' in line:  
                  print "line 1 = ", line
                  f.next()
                  each_volume_times_4 = []
                  each_volume_times_100 = []
                  
                  parameters = (''.join(islice(f, 1)))
                  columns = parameters.split()
                  each_volume = columns[6]
                  print 'each_volume = ', each_volume

                  VOLUME_EACH.append(each_volume)


        if 'MODES         EIGV          FREQUENCIES     IRREP' in line:

            print line
            print f.next()

            while True:
               target = f.next()
               aux = target.split()
               if not aux:
                  break

               first_No = aux[0]
               second_No = aux[1]
               freq = aux[3]
               print 'freq = ', freq
               print ' first_No original = ', first_No

               first_No = first_No.translate(None, '-')  # remove the '-'

               print ' first_No = ', first_No
               print 'second_No = ', second_No

               factor_freq = abs(int(second_No) - int(first_No)) + 1
               print 'factor_freq = ', factor_freq

               freqs = [freq] * factor_freq
               print 'freqs = ', freqs

               ALL_FREQ.append(freqs)

               n_freqs.append(freqs)
               print 'n_freqs = ', n_freqs

n_freqs = [item for sublist in n_freqs for item in sublist]
number_freqs_per_volume = len(n_freqs)

print 'Temps = ', Temps
print 'type(Temps)', type(Temps)
Temps = np.array(Temps)
print 'Temps[0] = ', Temps[0]
Temps = Temps.astype(np.float)
print 'Temps[0] = ', Temps[0]

#sys.exit()

print 'number_freqs_per_volume = ', number_freqs_per_volume
number_freqs_per_volume_minus_3 = number_freqs_per_volume - 3
print 'number_freqs_per_volume_minus_3 = ', number_freqs_per_volume_minus_3
#sys.exit()

print 'ALL_FREQ[1][:] = ', ALL_FREQ[1][:]
print 'VOLUME_EACH[1] = ', VOLUME_EACH[1]
#sys.exit()
ALL_FREQ = [item for sublist in ALL_FREQ for item in sublist]
print 'ALL_FREQ = ', ALL_FREQ
print 'len(ALL_FREQ) = ', len(ALL_FREQ)
#sys.exit()

thing = '0.0000'
while thing in ALL_FREQ: ALL_FREQ.remove(thing)
print 'ALL_FREQ = ', ALL_FREQ
print 'len(ALL_FREQ) = ', len(ALL_FREQ)

print 'ALL_FREQ[1] = ', ALL_FREQ[1][:]
print 'VOLUME_EACH[1] = ', VOLUME_EACH[1]

#print ' len(ALL_FREQ) = ', len(ALL_FREQ)
#ALL_FREQ = sorted(ALL_FREQ, key=float)  

output_array_2 = np.vstack((ALL_FREQ))#.T
np.savetxt('All_freq.dat', output_array_2, header="FREQS (CM^-1)", fmt="%s")
print 'np.shape(VOLUME_EACH)  = ', np.shape(VOLUME_EACH) 
print 'np.shape(E0)  = ', np.shape(E0) 
print 'np.shape(ET)  = ', np.shape(ET) 
print 'np.shape(TS)  = ', np.shape(TS) 
print 'np.shape(T)  = ', np.shape(T) 
#sys.exit()
# Transform each element of the list from <str> to <float64>:
VOLUME_EACH = [float(i) for i in VOLUME_EACH]
EL = [float(i) for i in EL]
E0 = [float(i) for i in E0]
ET = [float(i) for i in ET]
TS = [float(i) for i in TS]
T = [float(i) for i in T]
ALL_FREQ = [float(i) for i in ALL_FREQ]

# Transform each element of the list to a numpy array:
VOLUME_EACH = np.array(VOLUME_EACH)
EL = np.array(EL)
E0 = np.array(E0)
ET = np.array(ET)
TS = np.array(TS)
T = np.array(T)
ALL_FREQ = np.array(ALL_FREQ)

# Divide per F.U.:
VOLUME_EACH = VOLUME_EACH/n_F_u
EL = EL/n_F_u
E0 = E0/n_F_u
ET = ET/n_F_u
TS = TS/n_F_u

output_array = np.vstack((VOLUME_EACH, EL)).T
np.savetxt('EL_vs_V.dat', output_array, header="Volume           EL", fmt="%0.13f")
os.system("sort -k1 -n EL_vs_V.dat -o EL_vs_V.dat")

EL_plus_E0 = EL + E0

output_array = np.vstack((VOLUME_EACH, EL_plus_E0)).T
np.savetxt('EL_plus_E0_vs_V.dat', output_array, header="Volume           EL+E0", fmt="%0.13f")
os.system("sort -k1 -n EL_plus_E0_vs_V.dat -o EL_plus_E0_vs_V.dat")

n_volume = len(VOLUME_EACH)

n_T = len(T) / n_volume

print '#######'
print 'shape(ET) =  ', np.shape(ET)
print 'type(ET) =  ', type(ET)
print 'shape(ALL_FREQ) =  ', np.shape(ALL_FREQ)
ET = np.reshape(ET, (n_volume, n_T))
TS = np.reshape(TS, (n_volume, n_T))
T  = np.reshape(T, (n_volume, n_T))
ALL_FREQ = np.reshape(ALL_FREQ, (n_volume, number_freqs_per_volume_minus_3))
print 'type(ET) =  ', type(ET)
print 'shape(ET) =  ', np.shape(ET)
print 'shape(ALL_FREQ) =  ', np.shape(ALL_FREQ)

print 'number_freqs_per_volume_minus_3 = ', number_freqs_per_volume_minus_3
number_freqs_per_volume_minus_3_float = float(number_freqs_per_volume_minus_3)
print 'number_freqs_per_volume_minus_3_float = ', number_freqs_per_volume_minus_3_float

# Calculation of F_Anh :######
#############################

# KB = boltmann cte, KB = 1.38064852(79)x10-23 J/K
KB = 1.38064852E-23

# h = plank constant, h = 6.626070040(81)x10-34 J s
h = 6.626070040E-34

# c = speed of light, c = 2.99792458E8 m/s
c = 2.99792458E+8

Theta_each_V = []
rows = ALL_FREQ.shape[0]
print 'rows = ', rows
for i in range(rows):
  print 'ALL_FREQ[i][:] = ', ALL_FREQ[i][:]
  aux = (ALL_FREQ[i][:])**2
  print aux
  suma = sum(aux)
  print suma
  freq_prom = suma/number_freqs_per_volume_minus_3_float
  print freq_prom
  Theta = ((h/(2*np.pi)) / KB ) * (((5./3.) * (freq_prom))**0.5) * 1E+2 * c
  print Theta
  Theta_each_V.append(Theta)
print Theta_each_V
print VOLUME_EACH*2
#sys.exit()
print type(Theta_each_V)
Theta_each_V = np.array(Theta_each_V)
print type(Theta_each_V)

Ln_Theta_each_V = np.log(Theta_each_V)
Ln_VOLUME_EACH = np.log(VOLUME_EACH)

def BM(x, a, b, c, d):
         return  a + b*x + c*x**2 + d*x**3

def dBM(x, b, c, d):
         return  -b -2*c*x - 3*d*x**2


fig = plt.figure()

def autoscale_y(ax,margin=0.25):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>lo) & (xd<hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)

print 'Ln_VOLUME_EACH = ', Ln_VOLUME_EACH
print 'Ln_Theta_each_V = ', Ln_Theta_each_V

# sort the x values, otherwise the fit is not going to work:
p = Ln_VOLUME_EACH.argsort()
Ln_Theta_each_V = Ln_Theta_each_V[p]
Ln_VOLUME_EACH = Ln_VOLUME_EACH[p]

# Plotting the scattered points: 
p1 = plt.scatter(Ln_VOLUME_EACH, Ln_Theta_each_V, color='blue', marker="s", s=100)

# fit to cubic:
popt_CI, pcov_CI = curve_fit(BM, Ln_VOLUME_EACH, Ln_Theta_each_V)

# Plotting the fitting curves:
p2, = plt.plot(Ln_VOLUME_EACH, BM(Ln_VOLUME_EACH, *popt_CI), color='black', label='Cubic fit' )

autoscale_y(plt.gca())

print 'Ln_VOLUME_EACH = ', Ln_VOLUME_EACH
print 'Ln_Theta_each_V = ', Ln_Theta_each_V

plt.xlabel('$\ln V$ [$V =$ Angstrom$^{3}$ / F.U.]', fontsize=20)
plt.ylabel(r'$\ln \theta _{H}$ [$\theta _{H} = $ K]', fontsize=20)
#plt.suptitle("")
plt.title("B3LYP-D3", fontsize=20)
plt.ticklabel_format(useOffset=False)
plt.savefig('Ln_theta_vs_Ln_V.pdf', bbox_inches='tight')

print popt_CI 
print popt_CI[1:] 

Gamma = dBM(Ln_VOLUME_EACH, *popt_CI[1:])
print 'Gamma = ', Gamma

A_2 = ((3*KB)/(Theta_each_V)) * (0.0078*Gamma - 0.0154)
#A_2 = A_2/n_atoms_per_cell
print 'VOLUME_EACH = ', VOLUME_EACH
print 'A_2 = ', A_2 
#sys.exit()

#T = [10.0, 20.0]
#T = np.array(T)
#F_Anh_all
print 'Temps = ', Temps
for Ts in Temps:
 F_Anh = A_2[:] * Ts **2.0 * (1./4.3597482) * 1E+18
 print F_Anh
 output_array = np.vstack((VOLUME_EACH, F_Anh)).T
 output_array_sorted_on_V = output_array[output_array[:,0].argsort()]
 np.savetxt('F_Anh_vs_V_%0.2fK.dat'  %Ts, output_array_sorted_on_V, header="Volume           F Anh at %0.2fK" %Ts, fmt="%0.13f")
 os.makedirs('F_Anh_vs_V_%0.2fK' %Ts)
 shutil.move("./F_Anh_vs_V_%0.2fK.dat" %Ts, "./F_Anh_vs_V_%0.2fK" %Ts)

 V, F = np.loadtxt('../G_PT/F_vs_V_%0.2fK/F_vs_V_%0.2fK.dat' %(Ts, Ts), skiprows = 1).T
 F_tot = F + F_Anh
 output_array_2 = np.vstack((VOLUME_EACH, F_tot)).T

 np.savetxt('F_tot_vs_V_%0.2fK.dat'  %Ts, output_array_2, header="Volume           F + F Anh at %0.2fK" %Ts, fmt="%0.13f")
 os.makedirs('F_tot_vs_V_%0.2fK' %Ts)
 shutil.move("./F_tot_vs_V_%0.2fK.dat" %Ts, "./F_tot_vs_V_%0.2fK" %Ts)



os.system('rm -Rf EL_vs_V')
os.system('rm -Rf EL_plus_E0_vs_V')
os.system('rm -Rf G_PT')

#os.system('mv InTerSect_EL_level.py ./EL_vs_V')


