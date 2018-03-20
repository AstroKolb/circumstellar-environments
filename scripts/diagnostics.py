import numpy as np
import matplotlib.pylab as plt


path = '../output/Boots_diagnostics.dat'

data = np.loadtxt(path)


time = data[:,0]
mass = data[:,1]


plt.plot(time, mass)
plt.xlabel('time')
plt.ylabel('mass')
plt.show()