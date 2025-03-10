""" Post-processing of airfoil coefficients."""
from opensbli.postprocess.plot_functions import *
from opensbli.core.block import SimulationBlock
import multiprocessing
import numpy as np
# import numexpr as ne
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
    f, block_name, dsets, shape, constants = PP.read_block(time_group, block_number)
    # Extract surface slice
    gamma = constants['gama']
    Minf = constants['Minf']
    Re = constants['Re']
    SuthT = constants['SuthT']
    RefT = constants['RefT']
    simulation_time = constants['simulation_time']

    start = time.time()
    # Read Q vector
    rho = PP.read_full_dset(f, block_name, "rho_B%d" % block_number, verbose=False)
    rhou = PP.read_full_dset(f, block_name, "rhou0_B%d" % block_number, verbose=False)
    rhov = PP.read_full_dset(f, block_name, "rhou1_B%d" % block_number, verbose=False)
    rhoE = PP.read_full_dset(f, block_name, "rhoE_B%d" % block_number, verbose=False)

    u = rhou/rho
    v = rhov/rho
    p = (rhoE - 0.5*rho*(u**2 + v**2))*(gamma-1)

    end = time.time()
    #print("Time taken for reading Q vector: {:.2f}s".format(end - start))
    f.close()
    start = time.time()
    a = np.sqrt(gamma*p/rho)
    M = np.sqrt(u**2 + v**2)/a
    T = gamma*(Minf**2)*p/rho
    mu = T**(1.5)*(1.0+SuthT/RefT)/(T+SuthT/RefT)
    end = time.time()
    #print("Time taken for calculating quantities: {:.2f}s".format(end - start))
    return rho, u, v, rhoE, p, T, M, mu, constants

def calculate_angles(dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, y, ordering):
    if ordering == 'clockwise':
        ordering = 1
    else:
        ordering = -1
    """ Calculated on a single y plane."""
    Nx = np.size(dxdxi)
    a = detJ[0]*dz
    detJ = (dxdxi*dydeta-dxdeta*dydxi)*dz
    b = detJ[0]
    # print(a / b)
    aidet = 1.0/detJ
    dydx = dydxi / dxdxi
    ss, dss, th_, S_, S2_ = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)
    for i in range(0, Nx):
        th_ = np.arctan(dydx)
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
    var = variable[0:6, :]
    coeffs = np.array([-1.83333333333334, 3.00000000000002, -1.50000000000003, 0.333333333333356, -8.34617916606957e-15, 1.06910884386911e-15]) # want (6, Nz) shape
    coeffs = np.reshape(coeffs, (6, 1))
    dudy = np.zeros(Nx)
    init = var[:,:]*coeffs
    dudy[:] = np.sum(init, axis=0)/delta
    return dudy

def compute_skin_friction(tau_wall, Re):
    Cf = tau_wall/(0.5*Re)
    return Cf

def surface_derivatives_span(base_dir, time_group, block_number, dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, Ny, rho, u, v):
    """ Uses the metric relations to find the derivatives on the airfoil surface."""
    start = time.time()
    Nx = np.size(dxdxi)
    dudeta, dvdeta = np.zeros(Nx), np.zeros(Nx)
    # Take the points near the wall
    rho, u, v = rho[0:6,:], u[0:6,:], v[0:6,:]
    # Uniform derivative
    # for k in range(Nz):
    dudeta = compute_wall_normal_derivative(u[:,:], Ny)
    dvdeta = compute_wall_normal_derivative(v[:,:], Ny)

    dudeta = dudeta / (Ny - 1.)
    dvdeta = dvdeta / (Ny - 1.)
    # Recalculate detJ
    detJ = (dxdxi*dydeta-dxdeta*dydxi)*dz
    # # Transform the derivative
    dudy = np.zeros(Nx) # Nz, Nx ordering
    dvdx = np.zeros(Nx)

    dudy[:] = (1./detJ)*(dudeta*dxdxi)*dz
    dvdx[:] = (1./detJ)*(-dvdeta*dydxi)*dz
    end = time.time()
    #print("Time taken for reading derivative slices: {:.2f}s".format(end - start))
    return dudy, dvdx


def aerodynamic_coefficients_span(dudy, dvdx, p, mu, constants, ss, dss, th_, S_, S2_, aidet, dydx):
    start = time.time()
    Nx = p.shape[-1]
    gamma, Minf, Re = constants['gama'], constants['Minf'], constants['Re']
    # Freestream pressure
    pinf = 1.0 / (gamma*Minf*Minf)
    tau_wall, Cp = np.zeros(Nx), np.zeros(Nx)
    Cl, Cdp, Cdf = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)
    # Take wall values
    mu = mu[0,:]
    p = p[0,:]
    # First entry at i=0 location
    i = 0
    tau_wall[i] = S_[i]*mu[i]*(dudy[i]*np.abs(np.cos(th_[i]))-dvdx[i]*np.abs(np.sin(th_[i])))

    Cp[i] = (p[i]-pinf) / (gamma*0.5*pinf*Minf**2)

    # Temporary integration arrays
    for i in range(1,Nx):
        # print(i)
        # Size of the arc length, ds
        dlts = ss[i]-ss[i-1]
        # Lift coefficient
        fa = -S_[i]*(p[i]-pinf)*np.abs(np.cos(th_[i]))
        fb = -S_[i-1]*(p[i-1]-pinf)*np.abs(np.cos(th_[i-1])) 
        Cl[i] = Cl[i] + 0.5*(fa+fb)*dlts  / (0.5*1*1**2)
        # Pressure drag coefficient 
        fa = S2_[i]*p[i]*np.abs(np.sin(th_[i]))
        fb = S2_[i-1]*p[i-1]*np.abs(np.sin(th_[i-1]))
        Cdp[i] = Cdp[i] + 0.5*(fa+fb)*dlts / (0.5*1*1**2)
        # # Skin-friction drag coefficient
        fa = S_[i]*mu[i]*(dudy[i]*np.abs(np.cos(th_[i]))-dvdx[i]*np.abs(np.sin(th_[i]))) / Re
        fb = S_[i-1]*mu[i-1]*(dudy[i-1]*np.abs(np.cos(th_[i-1]))-dvdx[i-1]*np.abs(np.sin(th_[i-1]))) / Re
        Cdf[i] = Cdf[i] + 0.5*(fa+fb)*dlts / (0.5*1*1**2)
        # Wall shear stress for skin friction
        tau_wall[i] = S_[i]*mu[i]*(dudy[i]*np.abs(np.cos(th_[i]))-dvdx[i]*np.abs(np.sin(th_[i])))
        # Pressure coefficient
        Cp[i] = (p[i]-pinf) / (gamma*0.5*pinf*Minf**2)
    # Sum the 1D quantities over the x coordinate
    Cl = np.sum(Cl, axis=0)
    Cdp = np.sum(Cdp, axis=0)
    Cdf = np.sum(Cdf, axis=0)
    # Skin friction scaling
    Cf = compute_skin_friction(tau_wall, Re)
    end = time.time()
    #print("Time taken for aerodynamic coefficient calculation: {:.2f}s".format(end - start))
    return Cl, Cdp, Cdf, Cp, Cf


## Main calling function
def main(cases):
    """ Main processing loop, do each case in parallel."""
    # ordering = 'anti-clockwise'
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
    block_number = 1

    # Read the grid
    PP = OpenSBLIPreProcess()
    grid_file, block_name, dsets, shape, constants = PP.read_block(grid_file, block_number)
    halos = np.abs(grid_file[block_name][dsets[0]].attrs['d_p'])
    global nhalo 
    nhalo = halos[0]
    # Get the coordinates on the surface
    yloc = 0
    surface_line = np.s_[yloc + nhalo,startx:]
    x = PP.read_full_dset(grid_file, block_name, 'x0_B%d' % block_number, partial_slice=surface_line)
    y = PP.read_full_dset(grid_file, block_name, 'x1_B%d' % block_number, partial_slice=surface_line)
    # Number of grid points per direction (without halos) on this surface_line slice
    Nx, Ny = shape[-1], shape[0]
    grid_file.close()
        
    # Grid spacing
    dx = 0
    dy = 0
    dz = 1
    # Load metrics on the surface_line
    # Scaling of the derivatives used in the OpenSBLI code
    metrics_file = base_dir  + '/metrics.h5'
    metrics_file, block_name, dsets, shape, constants = PP.read_block(metrics_file, block_number)
    x_fact, y_fact = constants['Delta0block%d' % block_number], constants['Delta1block%d' % block_number]


    dxdxi = PP.read_full_dset(metrics_file, block_name, 'wk0_B%d' % block_number, partial_slice=surface_line) * x_fact
    dxdeta = PP.read_full_dset(metrics_file, block_name, 'wk1_B%d' % block_number, partial_slice=surface_line) * y_fact
    dydxi = PP.read_full_dset(metrics_file, block_name, 'wk2_B%d' % block_number, partial_slice=surface_line) * x_fact
    dydeta = PP.read_full_dset(metrics_file, block_name, 'wk3_B%d' % block_number, partial_slice=surface_line) * y_fact
    detJ = PP.read_full_dset(metrics_file, block_name, 'detJ_B%d' % block_number, partial_slice=surface_line)

    # Calculate angles around the surface
    ss, dss, th_, S_, S2_, aidet, dydx = calculate_angles(dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, y, ordering=ordering)

    # Get the unsteady input files
    # Extract data
    time_snapshots, iters = PP.find_files(base_dir)
    time_snapshots.sort(key=natural_keys)
    if len(time_snapshots) < 3:
        raise ValueError("Please post-process the results when there are more than 3 restart files available..")
    # Time-step
    dt = constants['dt']
    num_to_process = 10000
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
            rho, u, v, E, p, T, M, mu, constants = instantaneous_q(base_dir, time_group, block_number, constants)
            simulation_time = constants['simulation_time']
            # Get velocity derivatives at the surface
            dudy, dvdx = surface_derivatives_span(base_dir, time_group, block_number, dxdxi, dxdeta, dydxi, dydeta, detJ, dz, x, Ny, rho, u, v)
            # Calculate aerodynamic coefficients
            Cl, Cdp, Cdf, Cp, Cf = aerodynamic_coefficients_span(dudy, dvdx, p, mu, constants, ss, dss, th_, S_, S2_, aidet, dydx)
            print("Case: {:}, Output: {:}, Time: {:.5f}, Cl: {:.5f}, Cdp: {:.5f}, Cdf: {:.5f}, Cd: {:.5f}".format(case, time_group, simulation_time, np.mean(Cl), np.mean(Cdp), np.mean(Cdf), np.mean(Cdp) + np.mean(Cdf)))

            # Save the data
            store_coefficients[idx,0] = int(iters[idx])
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

cases = [('2D', 0.05)]
Mach_numbers = [0.2]
markers = ['^', 'o', 's']
new_labels = [r'$M=%.1f$' % M for M in Mach_numbers]
colors = ['b', 'g', 'r', 'k', 'm', 'c', 'y', '0.25', '0.75']

# Utility function
PP = OpenSBLIPreProcess()
skip_process = False
# Take only 1 spanwise plane, don't span average
#  Multi-processing of the cases
if not skip_process:
    with multiprocessing.Pool() as pool:
        pool.map(main, cases)

# Wait for all tasks to finish processing
    barrier = multiprocessing.Barrier(10)
## Plot all on same plot
plt.clf()

fontsize = 13

start_index = 0
end_index = -1
fig, ax1 = plt.subplots()

# Read these based on the case later
Minf = 0.6
gamma = 1.4
Cp_start, Cp_end = 0, -1

linewidth = 2

fig, ax1 = plt.subplots()

# Time average and plot Cp
for i, case in enumerate(cases):
    case = case[0]
    base_dir = '../../%s/' % case
    outPath = base_dir + 'time_history/'
    # Load x profile
    x = np.load(outPath + 'x_coords.npy')
    # Load Cp data
    Cp_data = np.load(outPath  + 'Cp_history.npy')
    # Cp = np.mean(Cp_data[start_index:end_index], axis=0)
    try:
        Cp = np.mean(Cp_data[Cp_start:Cp_end], axis=0)
    except:
        Cp = np.mean(Cp_data[0:], axis=0)

    print("Case: %s ---" % case, np.min(Cp), np.max(Cp))
    plt.plot(x, Cp, label='%s' % new_labels[i], linewidth=linewidth, color=colors[i])
    # Add reference data
    LJ = np.genfromtxt('../reference_data/LJ_reference_CP_%s.csv' % 'M02', delimiter=',')
    LJ_x, LJ_Cp = LJ[:,0], LJ[:,1]
    print(np.min(LJ_x), np.max(LJ_x))
    print(np.min(LJ_Cp), np.max(LJ_Cp))
    plt.plot(LJ_x, LJ_Cp, color='k', marker=markers[i], markerfacecolor='none', markeredgewidth=1, ls='', ms=9, markevery=3)

plt.legend(loc='lower right', prop={'size': 13})
# Add critical Cp value
plt.gca().invert_yaxis()
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlabel(r'$x$', fontsize=fontsize+2)
plt.ylabel(r'$C_p$', fontsize=fontsize+2)
plt.xlim([0,1])
plt.ylim([1,-0.8])
plt.grid()
ax1.text(-0.15, 0.95, '(a)',
    verticalalignment='bottom', horizontalalignment='right',
    transform=ax1.transAxes,
    color='black', fontsize=13)
plt.savefig('LJ_Re10000_validation_Cp.pdf', bbox_inches='tight')
plt.clf()


fig, ax1 = plt.subplots()

# Time average and plot Cf
# Separate the two sides
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
    try:
        Cf = np.mean(Cf_data[Cp_start:Cp_end], axis=0)
    except:
        Cf = np.mean(Cf_data[0:], axis=0)

    if i == 0:
        plt.axhline(y=0.0, color='k', linestyle='--')

    xp, xs = x[0:LE], x[LE-1:]
    xp = xp*-1
    plt.plot(xs, Cf[LE-1:], label='%s' % new_labels[i] ,color=colors[i], linewidth=linewidth)
    # Add reference data
    LJ = np.genfromtxt('../reference_data/LJ_reference_CF_%s.csv' % 'M02', delimiter=',')
    LJ_x, LJ_Cf = LJ[:,0], LJ[:,1]
    plt.plot(LJ_x, LJ_Cf, color='k', marker=markers[i], markerfacecolor='none', markeredgewidth=1, ls='', ms=9, markevery=3)
    plt.legend(loc='best', prop={'size': 13})

plt.xlabel(r'$x$', fontsize=fontsize+2)
plt.ylabel(r'$C_f$', fontsize=fontsize+2)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.xlim([-0.001, 0.99])

ax1.text(-0.15, 0.95, '(b)',
    verticalalignment='bottom', horizontalalignment='right',
    transform=ax1.transAxes,
    color='black', fontsize=13)

plt.grid()
plt.savefig('LJ_Re10000_validation_Cf.pdf', bbox_inches='tight')
plt.clf()


