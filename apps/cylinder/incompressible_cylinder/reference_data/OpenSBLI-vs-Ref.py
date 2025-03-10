# import numpy 
import numpy as np
# import matplotlib for plotting figures
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams['figure.figsize'] = 32, 20
plt.rcParams.update({'font.size': 40})
plt.rcParams.update({"text.usetex": True,"font.family": "serif","font.serif": ["Palatino"]})

mark_size = 20

#----------------------------------------------------------------------#
# read the probe file, made of float values
op = np.loadtxt("OpenSBLI_Re-Lsep-NLSt-NLGR.dat", dtype='f', skiprows=1)
gl = np.loadtxt("GL-Re-Lsep-NLSt-LGR.dat", dtype='f')
b  = np.loadtxt("B_Re-NLSt.dat", dtype='f')
#----------------------------------------------------------------------#
# OpenSBLI
reop   = op[:,0]
lsop   = op[:,1]
nlstop = op[:,2]
nlgrop = op[:,3]
# Giannetti and Luchini
regl   = gl[:,0]
lsgl   = gl[:,1]
lstgl  = gl[:,2]
lgrgl  = gl[:,3]
# Barkley and Henderson
reb    = b[:,0]
nlstb  = b[:,1]
#----------------------------------------------------------------------#
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(regl,lsgl,'kX',markersize=mark_size*2, label='Giannetti and Luchini (JFM, 2007)')
plt.plot(reop,lsop,'ro',markersize=mark_size, label='OpenSBLI')
plt.xlabel("Reynolds")
plt.ylabel("Base Flow Separation Length")
plt.legend(loc='best')
plt.grid(True)

plt.subplot(2,2,3)
plt.plot(reb, nlstb,'kX',markersize=int(mark_size*1.5), label='Barkley and Henderson (JFM, 1996)')
plt.plot(reop[1:],nlstop[1:],'ro',markersize=mark_size, label='OpenSBLI')
plt.xlabel("Reynolds")
plt.ylabel("Strouhal Number")
plt.legend(loc='best')
plt.grid(True)

plt.subplot(2,2,4)
plt.plot(regl[1:], lgrgl[1:],'kX',markersize=int(mark_size*1.5), label='Giannetti and Luchini (JFM, 2007)')
plt.plot(reop[1:],nlgrop[1:],'ro',markersize=mark_size, label='OpenSBLI')
plt.xlabel("Reynolds")
plt.ylabel("Linear growth Rate")
plt.legend(loc='best')
plt.grid(True)

plt.tight_layout()
plt.savefig("./cylinder-validation.png",dpi=300)
# show the plot
plt.show()
