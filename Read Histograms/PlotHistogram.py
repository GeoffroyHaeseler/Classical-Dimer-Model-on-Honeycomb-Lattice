import matplotlib.pyplot as plt
import numpy as np
import os

os.chdir("D:/3A/CF2/Results/OP Histogram/Test 4 - PSMN - BKT to 1st order transition point")
#plt.rcParams['text.usetex'] = True
#plt.style.use('classic')


magneticfield = '0.025'
Temperature = ['1.2','1.3','1.4']

bins = np.linspace(0, 0.5, 100)

for temperature in Temperature:
    FileName='H=' + magneticfield + '_L=90_T=' + temperature
    print(FileName)

    x = np.loadtxt(FileName + '.txt', usecols=(0))
    x = x[0:6500000]

    LegendName=' T=' + temperature
    plt.hist(x,bins, density=False, alpha=0.5, label=LegendName, edgecolor = 'k')

    locs, _ = plt.yticks()
    plt.yticks(locs,np.round(locs/len(x),3))


plt.xlabel(r'$\sqrt{3}\times\sqrt{3}$')
plt.xlim((0, 0.5))
plt.ylabel(r'$p(\sqrt{3}\times\sqrt{3})$')

plt.title('H='+magneticfield)

plt.legend(loc='upper right')

print('saving figure')

plt.savefig('H=' + magneticfield + '.png');
plt.clf();



#plt.rcParams['text.usetex'] = False
plt.clf();




