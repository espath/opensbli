from opensbli.postprocess.plot_functions import *
from opensbli.core.block import SimulationBlock
import multiprocessing
import numpy as np
import numexpr as ne
import h5py
import os, psutil
import matplotlib.style
import matplotlib as mpl
import time
global fontsize

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

## Useful functions, turn into a class later
def instantaneous_q(base_dir, time_group, block_number, constants):
    # Extract surface slice
    gamma = constants['gama']
    Minf = constants['Minf']
    Re = constants['Re']
    SuthT = constants['SuthT']
    RefT = constants['RefT']

    start = time.time()
    f = h5py.File(base_dir + 'J0.h5', 'r')
    # Read Q vector
    rho = f['opensbliblock0%d' % block_number][time_group]["rho_B%d" % block_number][()]
    rhou = f['opensbliblock0%d' % block_number][time_group]["rhou0_B%d" % block_number][()]
    rhov = f['opensbliblock0%d' % block_number][time_group]["rhou1_B%d" % block_number][()]
    rhow = f['opensbliblock0%d' % block_number][time_group]["rhou2_B%d" % block_number][()]
    rhoE = f['opensbliblock0%d' % block_number][time_group]["rhoE_B%d" % block_number][()]
    u = ne.evaluate('rhou/rho')
    v = ne.evaluate('rhov/rho')
    w = ne.evaluate('rhow/rho')
    p = ne.evaluate('(rhoE - 0.5*rho*(u**2 + v**2 + w**2))*(gamma-1)')

    end = time.time()
    #print("Time taken for reading Q vector: {:.2f}s".format(end - start))
    f.close()
    start = time.time()
    a = np.sqrt(gamma*p/rho)
    M = np.sqrt(u**2 + v**2 + w**2)/a
    T = ne.evaluate('gamma*(Minf**2)*p/rho')
    mu = ne.evaluate("T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)")
    # mu = T[:]
    end = time.time()
    #print("Time taken for calculating quantities: {:.2f}s".format(end - start))
    return rho, u, v, w, rhoE, p, T, M, mu, constants

def calculate_angles(dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, y, ordering):
    if ordering == 'clockwise':
        ordering = 1
    else:
        ordering = -1
    """ Calculated on a single y plane."""
    Nx = np.size(dxdxi)
    a = detJ[0]*dz
    detJ = (dxdxi*dydeta-dxdeta*dydxi)*dz ## Markus scaling
    b = detJ[0]
    # print(a / b)
    aidet = 1.0/detJ
    dydx = dydxi / dxdxi
    ss, dss, th_, S_, S2_ = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)
    for i in range(0, Nx):
        # th_ = np.arctan(dydx)*180.0 / np.pi
        th_ = np.arctan(dydx)
    ## Opposite of Lloyd Jones / Markus, because xi is in the other direction
    # S = 1 is pressure side, S = -1 is suction side
    for i in range(0,Nx):        
        if(dxdxi[i]>0): # dxdxi # suction side, anti-clockwise xi ordering
            S_[i] = ordering*1
        else:
            S_[i] = ordering*-1
        if(dydxi[i]>0):
            S2_[i]= ordering*1
        else:
            S2_[i]= ordering*-1
        # Calculate the arc lengths at each x[i] location
        if(i>0):
            ss[i] = ss[i-1]+((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)**0.5
            dss[i] = np.fabs(ss[i]-ss[i-1])
    return ss, dss, th_, S_, S2_, aidet, dydx

def compute_wall_normal_derivative(variable, Ny):
    Nz, Nx = variable.shape[0], variable.shape[-1]
    delta = 1.0/(Ny-1.0) ## no length scaling
    var = variable[:, 0:6, :]
    coeffs = np.array([-1.83333333333334, 3.00000000000002, -1.50000000000003, 0.333333333333356, -8.34617916606957e-15, 1.06910884386911e-15]) # want (6, Nz) shape
    coeffs = np.reshape(coeffs, (6, 1))
    dudy = np.zeros((Nz, Nx))
    init = var[:,:,:]*coeffs
    dudy[:,:] = np.sum(init, axis=1)/delta
    return dudy

def compute_skin_friction(tau_wall, Re):
    Cf = ne.evaluate('tau_wall/(0.5*Re)')
    return Cf

def surface_derivatives_span(base_dir, time_group, block_number, dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, Ny, Nz):
    """ Uses the metric relations to find the derivatives on the airfoil surface."""
    start = time.time()
    Nx = np.size(dxdxi)
    # Surface plane
    rho, u, v = np.zeros((Nz, 6, Nx)), np.zeros((Nz, 6, Nx)), np.zeros((Nz, 6, Nx))
    # Read the averaged values over the first 6 points normal to the wall
    for j in range(0, 6):
        f = h5py.File(base_dir + 'J%d.h5' % j, 'r')
        rho[:,j,startx:] = f['opensbliblock0%d' % block_number][time_group]["rho_B%d" % block_number][()]
        u[:,j,startx:] = f['opensbliblock0%d' % block_number][time_group]["rhou0_B%d" % block_number][()]/rho[:,j,startx:]
        v[:,j,startx:] = f['opensbliblock0%d' % block_number][time_group]["rhou1_B%d" % block_number][()]/rho[:,j,startx:]
        f.close()

    dudeta, dvdeta = np.zeros((Nz,Nx)), np.zeros((Nz,Nx))
    # Uniform derivative
    # for k in range(Nz):
    dudeta[:,:] = compute_wall_normal_derivative(u[:,:,:], Ny)
    dvdeta[:,:] = compute_wall_normal_derivative(v[:,:,:], Ny)

    dudeta = ne.evaluate('dudeta / (Ny - 1.)')
    dvdeta = ne.evaluate('dvdeta / (Ny - 1.)')
    # Recalculate detJ
    detJ = ne.evaluate('(dxdxi*dydeta-dxdeta*dydxi)*dz') ## Markus scaling
    # # Transform the derivative
    dudy = np.zeros((Nz,Nx)) # Nz, Nx ordering
    dvdx = np.zeros((Nz,Nx))

    # for k in range(Nz):
    dudy[:,:] = ne.evaluate('(1./detJ)*(dudeta*dxdxi)*dz')
    dvdx[:,:] = ne.evaluate('(1./detJ)*(-dvdeta*dydxi)*dz')
    end = time.time()
    #print("Time taken for reading derivative slices: {:.2f}s".format(end - start))
    return dudy, dvdx


def aerodynamic_coefficients_span(dudy, dvdx, p, mu, constants, ss, dss, th_, S_, S2_, aidet, dydx, Nz):
    start = time.time()
    Nx = p.shape[-1]
    gamma, Minf, Re = constants['gama'], constants['Minf'], constants['Re']
    # Freestream pressure
    pinf = 1.0 / (gamma*Minf*Minf)
    tau_wall, Cp = np.zeros((Nz,Nx)), np.zeros((Nz,Nx))
    Cl, Cdp, Cdf = np.zeros((Nz, Nx)), np.zeros((Nz,Nx)), np.zeros((Nz, Nx))
    # First entry at i=0 location
    i = 0
    tau_wall[:,i] = S_[i]*mu[:,i]*(dudy[:,i]*np.abs(np.cos(th_[i]))-dvdx[:,i]*np.abs(np.sin(th_[i])))
    Cp[:,i] = (p[:,i]-pinf) / (gamma*0.5*pinf*Minf**2)

    # Temporary integration arrays
    fa, fb = np.zeros(Nz), np.zeros(Nz)
    for i in range(1,Nx):
        # print(i)
        # Size of the arc length, ds
        dlts = ss[i]-ss[i-1]
        # Lift coefficient
        fa[:] = -S_[i]*(p[:,i]-pinf)*np.abs(np.cos(th_[i]))
        fb[:] = -S_[i-1]*(p[:, i-1]-pinf)*np.abs(np.cos(th_[i-1])) 
        Cl[:,i] = Cl[:,i] + 0.5*(fa+fb)*dlts  / (0.5*1*1**2)
        # Pressure drag coefficient 
        fa[:] = S2_[i]*p[:,i]*np.abs(np.sin(th_[i]))
        fb[:] = S2_[i-1]*p[:,i-1]*np.abs(np.sin(th_[i-1]))
        Cdp[:,i] = Cdp[:,i] + 0.5*(fa+fb)*dlts / (0.5*1*1**2)
        # # Skin-friction drag coefficient
        fa[:] = S_[i]*mu[:,i]*(dudy[:,i]*np.abs(np.cos(th_[i]))-dvdx[:,i]*np.abs(np.sin(th_[i]))) / Re
        fb[:] = S_[i-1]*mu[:,i-1]*(dudy[:,i-1]*np.abs(np.cos(th_[i-1]))-dvdx[:,i-1]*np.abs(np.sin(th_[i-1]))) / Re
        Cdf[:,i] = Cdf[:,i] + 0.5*(fa+fb)*dlts / (0.5*1*1**2)
        # Wall shear stress for skin friction
        tau_wall[:,i] = S_[i]*mu[:,i]*(dudy[:,i]*np.abs(np.cos(th_[i]))-dvdx[:,i]*np.abs(np.sin(th_[i])))
        # Pressure coefficient
        Cp[:,i] = (p[:,i]-pinf) / (gamma*0.5*pinf*Minf**2)
    # Sum the 1D quantities over the x coordinate
    Cl = np.sum(Cl, axis=1)
    Cdp = np.sum(Cdp, axis=1)
    Cdf = np.sum(Cdf, axis=1)
    # Skin friction scaling
    Cf = compute_skin_friction(tau_wall, Re)
    end = time.time()
    #print("Time taken for aerodynamic coefficient calculation: {:.2f}s".format(end - start))
    return Cl, Cdp, Cdf, Cp, Cf


## Main calling function
def main(cases):
    """ Main processing loop, do each case in parallel."""
    ordering = 'clockwise'

    Lz = cases[1]
    case = cases[0]

    base_dir = '../../%s/' % case

    Lx, Ly = 2.0461756979465546, 22.5

    mblock = True
    outPath = base_dir + 'time_history/'
    global startx
    startx = 0
    try:
        os.mkdir(outPath)
    except FileExistsError:
        pass

    grid_file = base_dir + '/data.h5'

    if mblock:
        block_number = 1
    else:
        block_number = 0
    # Read the grid
    grid_file, block_name, dsets, shape, constants = PP.read_block(grid_file, block_number)
    halos = np.abs(grid_file[block_name][dsets[0]].attrs['d_p'])
    global nhalo 
    nhalo = halos[0]
    # Get the coordinates on the surface
    yloc = 0
    surface_line = np.s_[10, yloc + nhalo,startx:]
    x = PP.read_full_dset(grid_file, block_name, 'x0_B%d' % block_number, partial_slice=surface_line)
    y = PP.read_full_dset(grid_file, block_name, 'x1_B%d' % block_number, partial_slice=surface_line)
    z = PP.read_full_dset(grid_file, block_name, 'x2_B%d' % block_number, partial_slice=surface_line)

    # Check airfoil length
    dist_array = (x[:-1]-x[1:])**2 + (y[:-1]-y[1:])**2
    print("Airfoil length: {}".format(np.sum(np.sqrt(dist_array))))
    # Number of grid points per direction (without halos) on this surface_line slice
    Nx, Ny, Nz = shape[-1], shape[1], shape[0]
    grid_file.close()
        
    zloc = int(Nz/2.0) + nhalo
    # Grid spacing
    dx = 0
    dy = 0
    dz = Lz / float(Nz)
    # Load metrics on the surface_line
    # Scaling of the derivatives used in the OpenSBLI code
    metrics_file = base_dir  + '/metrics.h5'
    metrics_file, block_name, dsets, shape, constants = PP.read_block(metrics_file, block_number)
    x_fact, y_fact = constants['Delta0block%d' % block_number], constants['Delta1block%d' % block_number]

    dxdxi = PP.read_full_dset(metrics_file, block_name, 'wk0_B%d' % block_number, partial_slice=surface_line) * x_fact
    dxdeta = PP.read_full_dset(metrics_file, block_name, 'wk1_B%d' % block_number, partial_slice=surface_line) * y_fact
    dydxi = PP.read_full_dset(metrics_file, block_name, 'wk3_B%d' % block_number, partial_slice=surface_line) * x_fact
    dydeta = PP.read_full_dset(metrics_file, block_name, 'wk4_B%d' % block_number, partial_slice=surface_line) * y_fact
    detJ = PP.read_full_dset(metrics_file, block_name, 'detJ_B%d' % block_number, partial_slice=surface_line)

    # Calculate angles around the surface
    ss, dss, th_, S_, S2_, aidet, dydx = calculate_angles(dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, y, ordering=ordering)

    # Get the unsteady input files
    input_slice_file = base_dir + 'J0.h5'
    # Find the time_snapshots
    f0 = h5py.File(input_slice_file, 'r')
    time_snapshots = [x for x in f0['opensbliblock0%d' % block_number].keys() if x != "0"]
    # Sort based on index number
    time_snapshots.sort(key=natural_keys)

    f0.close()

    dt = constants['dt']
    simulation_times = [dt*float(x) for x in time_snapshots]  
    num_to_process = -1

    time_snapshots = time_snapshots[0:num_to_process] # skip most recent one?
    nsamples = len(time_snapshots)

    # Output log
    log_file = outPath + 'aero_coefficients2.dat'
    # Cp curve 1D arrays, in time
    Cp_profiles = outPath + 'Cp_history'
    Cf_profiles = outPath + 'Cf_history'
    x_coords = outPath + 'x_coords'
    # Columns to save
    nvars = 6

    # Resume processing or start from beginning if the file doesn't exist
    if os.path.exists(log_file) and os.path.exists(Cp_profiles + '.npy'):
        load_coefficients = np.genfromtxt(log_file)
        load_Cp = np.load(Cp_profiles + '.npy')
        load_Cf = np.load(Cf_profiles + '.npy')
        last_processed = np.shape(load_coefficients)[0]
        store_coefficients = np.zeros((nsamples, nvars))
        store_Cp = np.zeros((nsamples, Nx))
        store_Cf = np.zeros((nsamples, Nx))
        # Append previously processed to the start of the new one
        store_coefficients[0:last_processed,:] = load_coefficients[:,:]
        store_Cp[0:last_processed,:] = load_Cp[:,:]
        store_Cf[0:last_processed,:] = load_Cf[:,:]

    else:
        print("Fresh start.")
        store_coefficients = np.zeros((nsamples, nvars))
        store_Cp = np.zeros((nsamples, Nx))
        store_Cf = np.zeros((nsamples, Nx))
        last_processed = -1

    # Begin looping over the output files, skip any that are already processed
    print("Last processed: ", last_processed)
    for idx, time_group in enumerate(time_snapshots[0:nsamples]):
        if idx > (last_processed - 1):
            if idx % 10 == 0:
                print("{} files remaining to process.".format(nsamples - idx))
                print("Memory usage:", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2) 
            rho, u, v, w, E, p, T, M, mu, constants = instantaneous_q(base_dir, time_group, block_number, constants)
            # simulation_time = constants['simulation_time']
            simulation_time = simulation_times[idx]
            # Get velocity derivatives at the surface
            dudy, dvdx = surface_derivatives_span(base_dir, time_group, block_number, dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, Ny, Nz)
            # Calculate aerodynamic coefficients
            Cl, Cdp, Cdf, Cp, Cf = aerodynamic_coefficients_span(dudy, dvdx, p, mu, constants, ss, dss, th_, S_, S2_, aidet, dydx, Nz)
            print("Case: {:}, Output: {:}, Time: {:.5f}, Cl: {:.5f}, Cdp: {:.5f}, Cdf: {:.5f}, Cd: {:.5f}".format(case, time_group, simulation_time, np.mean(Cl), np.mean(Cdp), np.mean(Cdf), np.mean(Cdp) + np.mean(Cdf)))
            # Average over the span
            Cp = np.mean(Cp, axis=0)
            Cf = np.mean(Cf, axis=0)

            # Save the data
            store_coefficients[idx,0] = int(time_snapshots[idx])
            store_coefficients[idx,1] = simulation_time
            store_coefficients[idx,2] = np.mean(Cl)
            store_coefficients[idx,3] = np.mean(Cdp)
            store_coefficients[idx,4] = np.mean(Cdf)
            store_coefficients[idx,5] = np.mean(Cdp) + np.mean(Cdf)
            store_Cp[idx, :] = Cp
            store_Cf[idx, :] = Cf

    # Output updated values to disk
    np.savetxt(log_file, store_coefficients)
    # Save Cp profiles
    np.save(Cp_profiles, store_Cp)
    # Save Cf profiles
    np.save(Cf_profiles, store_Cf)
    # Save x coordinates
    np.save(x_coords, x)

cases = [('/3D/', 0.05)]

grid_sizes = [pc[1] for pc in cases]
AoA = [5 for _ in range(len(cases))]

labels = [r'$AR = %.3f$, $\alpha = %d^{\circ}$' % (i, b) for (i,b) in zip(grid_sizes, AoA)]
colors = ['k', 'b' 'r', 'g', 'm', 'saddlebrown', 'c', 'y']

def select_data(case):
    # Plot up to this time, excluding the transient
    transient_time = 0
    final_time = 100
    base_dir = '../../%s/' % case
    outPath = base_dir + 'time_history/'
    log_file = outPath + 'aero_coefficients2.dat'
    data = np.genfromtxt(log_file)
    times = data[:,1]
    # Remove transient if there is one
    try:
        start_index = np.where(times > transient_time)[0][0]
    except:
        start_index = 0
    try:
        end_index = np.where((times - times[0]) > final_time)[0][0] - 1 + start_index
    except:
        end_index = times.shape[-1] + 2
    if start_index == 0:
        end_index = end_index - 1
    times = times[start_index:end_index]
    times = times - times[0]
    # Reset start time to zero
    print("For case: {}, using {} samples. Start/End Index: {}, {}. Start Time: {}, End Time: {}".format(case, end_index-start_index, start_index, end_index, times[0], times[-1]))
    return start_index, end_index, times, final_time

# Utility function
PP = OpenSBLIPreProcess()

skip_process = False

#  Multi-processing of the cases
if not skip_process:
    with multiprocessing.Pool() as pool:
        pool.map(main, cases)

    # Wait for all tasks to finish processing
    barrier = multiprocessing.Barrier(10)

## Plot all on same plot
plt.clf()

fontsize = 13
plt.rcParams.update({'font.size': fontsize})
plt.rcParams.update({"text.usetex": True,"font.family": "serif","font.serif": ["Palatino"]})

# fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(6.4*2, 4.8*2))
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 6.5))
ax0, ax1, ax2, ax3 = axes[0,0], axes[0,1], axes[1,0], axes[1,1]


for i, case in enumerate(cases):
    case = case[0]
    base_dir = '../../%s/' % case
    outPath = base_dir + 'time_history/'

    log_file = outPath + 'aero_coefficients2.dat'
    data = np.genfromtxt(log_file)

    times = data[:,1]
    Cl = data[:,2]
    Cd = data[:,5]
    Cdp = data[:,3]
    Cdf = data[:,4]
    nsamples = Cl.shape[0]
    start_index, end_index, times, final_time = select_data(case)
    Cl, Cd, Cdp, Cdf = Cl[start_index:end_index], Cd[start_index:end_index], Cdp[start_index:end_index], Cdf[start_index:end_index]
    print("Case: {}, mean Cl: {:.3f}, mean Cl/Cd: {:.3f}, mean Cd: {:.4f}, mean Cdp: {:.4f}, mean Cdf: {:.4f}.".format(case, np.mean(Cl), np.mean(Cl) / np.mean(Cd), np.mean(Cd), np.mean(Cdp), np.mean(Cdf)))
    # Add mean line
    ax0.plot(times, Cl, label='%s : $\overline{C_L} : %.3f$' % (labels[i], np.mean(Cl)), color=colors[i])
    ax0.axhline(y=np.mean(Cl), ls='--', color=colors[i])


# ax0.legend(prop={'size': fontsize-3}, ncol=1, loc='best')
ax0.set_xlabel(r'$t - t_0$')
ax0.set_ylabel(r'$C_L$')
ax0.set_xlim([0, final_time])
# ax0.set_ylim([0.825, 1.17])
ax0.grid()
# These should be read from the constants in the output file
Minf = 0.72
gamma = 1.4

# Time average and plot Cp
for i, case in enumerate(cases):
    case = case[0]
    base_dir = '../../%s/' % case
    outPath = base_dir + 'time_history/'
    # Load x profile
    x = np.load(outPath + 'x_coords.npy')
    # Load Cp data
    Cp_data = np.load(outPath  + 'Cp_history.npy')
    start_index, end_index, times, final_time = select_data(case)
    Cp = np.mean(Cp_data[start_index:end_index], axis=0)
    connection = (Cp[0]+Cp[-1])/2.0
    Cp[0], Cp[-1] = connection, connection

    print("Case: %s ---" % case, np.min(Cp), np.max(Cp))

    if i == 0:
         Cp_crit = (2./(gamma*Minf*Minf)) * ( ((2. + (gamma - 1.)*Minf*Minf)/(gamma + 1.))**(gamma / (gamma - 1.)) - 1.) ## WARNING, Reload constants here
         ax2.axhline(y=Cp_crit, color='k', linestyle='--')

    ax2.plot(x, Cp, label='%s' % labels[i], color=colors[i])

# Add critical Cp value
ax2.invert_yaxis()
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$C_p$')
ax2.set_xlim([0,1])
ax2.grid()

# Time average and plot Cf
LE = 499
for i, case in enumerate(cases):
    case = case[0]
    base_dir = '../../%s/' % case
    outPath = base_dir + 'time_history/'
    # Load x profile
    x = np.load(outPath + 'x_coords.npy')
    # Load Cp data
    Cf_data = np.load(outPath  + 'Cf_history.npy')
    nsamples = Cf_data.shape[0]
    start_index, end_index, times, final_time = select_data(case)
    Cf = np.mean(Cf_data[start_index:end_index], axis=0)
    if i == 0:
        plt.axhline(y=0.0, color='k', linestyle='--')

    xp, xs = x[0:LE], x[LE-1:]
    xp = xp*-1
    ax3.plot(xp, Cf[0:LE], label='%s' % labels[i], color=colors[i])
    ax3.plot(xs, Cf[LE-1:], color=colors[i])
    # ax3.legend(ncol=2, prop={'size': fontsize-3}, loc='lower center')

ax3.set_xlabel(r'$x^{\prime}$')
ax3.set_ylabel(r'$C_f$')
ax3.set_xlim([-0.99, 0.98])
ax3.set_ylim([-0.001, 0.01])
ax3.set_xticks([-1, 0, 1])

ax3.text(0.325, 0.9, 'Pressure Side',
    verticalalignment='bottom', horizontalalignment='right',
    transform=ax3.transAxes,
    color='black')

ax3.text(0.95, 0.9, 'Suction Side',
    verticalalignment='bottom', horizontalalignment='right',
    transform=ax3.transAxes,
    color='black')

plt.grid()

import numpy as np
import matplotlib.pylab as plt
from scipy import signal

def compute_PSD(x, q):

  q_prime = q - np.mean(q)
  xs      = 1/(x[1]-x[0])
  x, PSD  = signal.periodogram(q_prime[0:-1], xs, window='hann', scaling='spectrum')
  ## use the following to get the right amplitudes
  #PSD     = 2 * np.sqrt(PSD / 2)
  # use the following to get PSD
  # note 1: the scaling='spectrum' does not normalize with xs
  # note 2: the multiplication by 2 accounts for the one-sided spectrum
  PSD     = 2 * PSD / xs
  max_idx = np.argmax(PSD)
  return x, PSD, max_idx

def compute_PSD_Welch(x, q, segment=2, overlap=0.9):

  nPseg   = int(np.ceil(np.size(q)/segment))
 
  q_prime = q - np.mean(q)
  xs      = 1/(x[1]-x[0])
  x, PSD  = signal.welch(q_prime[0:-1], xs, window='hann', nperseg=nPseg, noverlap=nPseg*overlap,scaling='spectrum')
  ## use the following to get the right amplitudes
  #PSD     = 2 * np.sqrt(PSD / 2)
  # use the following to get PSD
  # note 1: the scaling='spectrum' does not normalize with xs
  # note 2: the multiplication by 2 accounts for the one-sided spectrum
  PSD     = 2 * PSD / xs
  max_idx = np.argmax(PSD)
  return x, PSD, max_idx


for i, case in enumerate(cases):
    case = case[0]
    base_dir = '../../%s/' % case
    # OpenSBLI data
    outPath = base_dir + 'time_history/'
    log_file = outPath + 'aero_coefficients2.dat'
    data = np.genfromtxt(log_file)
    Cl = data[:,2]
    Cd = data[:,5]
    Cdp = data[:,3]
    Cdf = data[:,4]
    nsamples = Cl.shape[0]
    start_index, end_index, times, final_time = select_data(case)
    Cl, Cd, Cdp, Cdf = Cl[start_index:end_index], Cd[start_index:end_index], Cdp[start_index:end_index], Cdf[start_index:end_index]
    # Time period considered
    sf = 1 / (times[1] - times[0])
    # print("Available samples: {}".format(nsamples))
    # Mean Cl over the interval
    mean_Cl  = np.mean(Cl)
    # Mean Cd over the interval
    mean_Cd  = np.mean(Cd)
    print("Case: {}, Mean CL: {:.3f}, Mean CD: {:.3f}, Sampling Frequency: {:.3f}".format(case, mean_Cl, mean_Cd, sf))
    # Save the FFT data
    ff, PP, max_idx = compute_PSD(times, Cl)
    ff_w, PP_w, max_idx_w = compute_PSD_Welch(times, Cl)
    print('OpenSBLI (ILES) - St: {:.4f}, PSD: {:}.'.format(ff_w[max_idx_w], PP_w[max_idx_w]))
    # Make the plot
    #plt.loglog(ff, PP, label=labels[i] + ' : ' r'$S_t = %s$' % str(np.around(ff[max_idx], 4)), color=colors[i])
    ax1.loglog(ff,   PP,   linestyle='solid',  label=labels[i] , color=colors[i])

ax1.grid(True, which="both")
ax1.set_xlabel(r'$S_t$')
ax1.set_ylabel(r'$PSD\left(C_{L}^{\prime}\right)$')
ax1.set_xlim([1/(times[-1]-times[0]), 0.5*sf])

# Add letters
letters = ['a', 'b', 'c', 'd']
offsets = [-0.13, -0.13, -0.13, -0.13]
axes = [ax0, ax1, ax2, ax3]
for i, my_ax in enumerate(axes):
    my_ax.text(offsets[i], 0.95, '(%s)' % letters[i],
        verticalalignment='bottom', horizontalalignment='right',
        transform=my_ax.transAxes,
        color='black')

plt.tight_layout()
ax0.legend(prop={'size': fontsize}, 
           loc='upper center', bbox_to_anchor=(1.1, 1.2),frameon=False, ncol=4)
plt.subplots_adjust(hspace=0.25)
plt.savefig('Airfoil_coefficients.pdf', bbox_inches='tight')


