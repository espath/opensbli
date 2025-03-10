# Read a mixing layer dataset and make a contour plot
# Runs on spyder (p updateython 3.8), one deprecated feature to uodate (dataset)
# NDS 21/1/2021
#
import numpy as np
import h5py
import matplotlib.pyplot as plt
gama=1.4
Mref=0.4  # needs changing here whenever M changes - better to read it from the file
print('Check: Mach number = ',Mref)

def read_dataset(file, dataset):
    group = file["opensbliblock00"]
    d_m = group["%s" % (dataset)].attrs['d_m']
    size = group["%s" % (dataset)].shape
    start=[abs(d+2) for d in d_m]
    end=[s-abs(d+2) for d, s in zip(d_m, size)]
    read_data=group["%s" % (dataset)][start[0]:end[0],start[1]:end[1]]
    return read_data

print('Reading data')
#fname='opensbli_output.h5'
fname='opensbli_output_001000.h5'
ff=h5py.File(fname, 'r')
x0dum=read_dataset(ff,'x0_B0')
x1dum=read_dataset(ff,'x1_B0')
rdum=read_dataset(ff, 'rho_B0')
rudum=read_dataset(ff, 'rhou0_B0')
rvdum=read_dataset(ff, 'rhou1_B0')
rEdum=read_dataset(ff, 'rhoE_B0')
rfdum=read_dataset(ff, 'rhof_B0')

print('Setting flowfield arrays')
x=x0dum[2:-2,2:-2]
y=x1dum[2:-2,2:-2]
rho=rdum[2:-2,2:-2]
rhou=rudum[2:-2,2:-2]
rhov=rvdum[2:-2,2:-2]
rhoE=rEdum[2:-2,2:-2]
rhof=rfdum[2:-2,2:-2]
u=rhou/rho
v=rhov/rho
e=rhoE/rho-0.5*(u**2+v**2)
p=(gama-1.0)*rho*e
T=e*Mref**2*gama*(gama-1.0)
f=rhof/rho

# note that the first array index is y, the second x (python convention)
print('Array size for plotting',rho.shape)
print('x range',x[0,0],x[0,-1])
print('y range',y[0,0],y[-1,0])
print('u range',np.amax(u),np.amin(u))
print('v range',np.amax(v),np.amin(v))
print('f range',np.amax(f),np.amin(f))

fig1,ax=plt.subplots()
CS=ax.contourf(x,y,f,20)
ax.set_aspect('equal')
cbar=fig1.colorbar(CS)
ax.set_ylim([-7.5,7.5])
plt.savefig('mixlay_plot.pdf')
