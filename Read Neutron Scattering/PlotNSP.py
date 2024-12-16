from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import math
import os

os.chdir("Read Neutron Scattering")
plt.rcParams['text.usetex'] = True
plt.style.use('classic')

color=1/np.logspace(0.0005,math.e,256)
transparency=np.zeros(256)
customcolor=np.vstack((color, color, color, transparency+1))
customcolor=np.transpose(customcolor)
newcmp = ListedColormap(customcolor)

cutoff = 0.5

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

images =[]

FileName = 'T=0.1'
print(FileName)

NSP = np.loadtxt(FileName + '.txt')

NSP[:, 0] = NSP[:, 0] / (2 * math.pi)
NSP[:, 1] = NSP[:, 1] / (2 * math.pi)
NSP[:, 4] = NSP[:, 4]

PlottedData = np.argwhere(np.sqrt(NSP[:,0]**2+NSP[:,1]**2) < 3)
x0=NSP[PlottedData, 0]
y0=NSP[PlottedData, 1] / math.sqrt(3)
z0=NSP[PlottedData, 4]

FirstScale  = np.argwhere( z0[:] < cutoff )
x1=x0[FirstScale]
y1=y0[FirstScale]
z1=z0[FirstScale]

SecondScale = np.argwhere( z0[:] >= cutoff )
x2=x0[SecondScale]
y2=y0[SecondScale]
z2=z0[SecondScale]

PCMa = axes[0].scatter(x1, y1, c=z1, vmin=0     , vmax=cutoff, s=1, cmap="CMRmap", edgecolor="none")
PCMb = axes[0].scatter(x2, y2, c=z2, vmin=cutoff, vmax=1000  , s=1, cmap=newcmp, edgecolor="none")
#images.append([PCMa, PCMb ])
images.append([PCMa ])

axes[0].set_xlim((-3, 3))
axes[0].set_ylim((-math.sqrt(3), math.sqrt(3)))
axes[0].set_ylabel(r'$(k,k,2\bar{kh})$')
axes[0].set_xlabel(r'$(h,\bar{h},0)$')
axes[0].set_title('NSP')

fig.colorbar(       PCMb, ax=[axes[1]], location = 'left')
cbar = fig.colorbar(PCMa, ax=[axes[1]], location = 'right')
cbar.formatter.set_powerlimits((0, 0))

axes[1].remove()

plt.subplots_adjust(wspace=0, hspace=0)

plt.rcParams['text.usetex'] = False

plt.savefig(FileName + '.png')
plt.clf()
plt.close('all')
