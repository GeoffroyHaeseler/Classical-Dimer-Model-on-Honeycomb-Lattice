import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
import os
import matplotlib.colors as mcolors

os.chdir("D:/3A/CF2")
plt.rcParams['text.usetex'] = True
plt.style.use('classic')



LegendList=[]

File_Name = 'evolution'
print(File_Name)
x1 = np.loadtxt(File_Name + '.txt', skiprows=1, usecols=(0))
y1 = np.loadtxt(File_Name + '.txt', skiprows=1, usecols=(7))
y1=y1/x1
LegendName = 'new PBC L=9 Heated AllTS'
LegendList.append(LegendName)
print(trapz(y1,x=x1))

Y1=[]
for i in range (len(x1)):
    Y1.append(trapz(y1[0:i],x=x1[0:i]))
plt.plot(x1, Y1)


#plt.ylim(0,0.2)

plt.xlabel("T")
plt.ylabel(r'$S$')

plt.legend(LegendList, prop={'size': 10})

plt.rcParams['text.usetex'] = False




print('saving figure')
plt.savefig('test.png');
plt.clf();



