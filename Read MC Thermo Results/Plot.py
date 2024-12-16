from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import math
import os

os.chdir("D:/3A/CF2/Results/Thermo/Test 25 - Phycalc2 - Low temperature behaviour at the transition between 1st and Kasteleyn transition")
plt.rcParams['text.usetex'] = True
plt.style.use('classic')

SystemSize=['18', '24', '36', '48', '60', '72', '96']

# What do you want to plot?
# <E>       5
# Cb        7
# r3r3      9
# χ         11
# |m_z|     13
# χm_z      14
# nA        28
# nB        29
# nC        30

Data=13
if Data == 5:
    title='E'
    DataName='<E>'

elif Data == 7:
    title='Cb'
    DataName='Cb'

elif Data == 9:
    title='BinderCumulent'
    DataName=r'$V_L$'

elif Data == 10:
    title='m_x'
    DataName=r'$m_x$'

elif Data == 11:
    title='m_y'
    DataName=r'$m_y$'

elif Data == 13:
    title='m'
    DataName=r'$<m>$'

elif Data == 15:
    title='χ'
    DataName=r'$\chi$'

elif Data == 16:
    title='Tχ'
    DataName='Tχ'

elif Data == 17:
    title='r3r3'
    DataName=r'$\sqrt{3}\times\sqrt{3}$'

elif Data == 19:
    title='χr3r3'
    DataName=r'$\chi_{\sqrt{3}\times\sqrt{3}}$'

elif Data == 20:
    title='r3r3BinderCumulent'
    DataName=r'$V_{\sqrt{3}\times\sqrt{3}}$'


for Length in SystemSize:

    File_Name = 'L=' + Length
    print(File_Name)
    x = np.loadtxt(File_Name + '.txt', skiprows=1, usecols=(2))
    y = np.loadtxt(File_Name + '.txt', skiprows=1, usecols=(Data))

    plt.plot(x, y, '-+' , label=File_Name)


plt.legend(loc='lower left')

plt.xlabel('Temperature (K)')
plt.ylabel(DataName)
#plt.xlim(0,2.5)
plt.rcParams['text.usetex'] = False

plt.savefig(title + '.png')
plt.clf()
plt.close('all')