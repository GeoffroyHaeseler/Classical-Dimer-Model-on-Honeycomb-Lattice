from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import math
import os

os.chdir("D:/3A/CF2/Results/line/Test 1 - PSMN - B=0=V")
plt.rcParams['text.usetex'] = True
plt.style.use('classic')

FileName = 'line.txt'
print(FileName)

NSP = np.loadtxt(FileName)

NSP[:, 0] = NSP[:, 0] / (2 * math.pi)
NSP[:, 1] = NSP[:, 1] / (2 * math.pi)
NSP[:, 4] = NSP[:, 4]

PlottedData = np.argwhere(np.sqrt(NSP[:,0]==0))
x0=NSP[PlottedData, 1] / math.sqrt(3)
y0=NSP[PlottedData, 4]


plt.plot(x0,y0)



plt.xlabel('k k -2k')
plt.ylabel('I')
plt.xlim(-math.sqrt(3),math.sqrt(3))
#plt.xlim(0,0.5)

plt.savefig('line.png')
plt.clf()
plt.close('all')