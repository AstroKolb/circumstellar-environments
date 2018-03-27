import numpy as np
import matplotlib.pylab as plt


path = '../output/launch4_diagnostics.dat'

data = np.loadtxt(path)


time  = data[:,0]
mass1 = data[:,1]
mass2 = data[:,2]
mass3 = data[:,3]
mass4 = data[:,4]
flux1 = data[:,5] * 3.15e7/2e33
flux2 = data[:,6] * 3.15e7/2e33
flux3 = data[:,7] * 3.15e7/2e33
flux4 = data[:,8] * 3.15e7/2e33

plt.plot(time, mass1+mass2+mass3+mass4)
plt.xlabel('time')
plt.ylabel('mass')
plt.show()


plt.plot(time, mass1/mass1[-1], label='mass1')
plt.plot(time, mass2/mass2[-1], label='mass2')
plt.plot(time, mass3/mass3[-1], label='mass3')
plt.plot(time, mass4/mass4[-1], label='mass4')
plt.xlabel('time')
plt.ylabel('mass')
plt.legend()
plt.show()

plt.plot(time, flux1, label='flux1')
plt.plot(time, flux2, label='flux2')
plt.plot(time, flux3, label='flux3')
plt.plot(time, flux4, label='flux4')
plt.xlabel('time')
plt.ylabel('flux')
plt.legend()
plt.show()
