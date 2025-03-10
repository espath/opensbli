import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.pyplot as plt
import h5py
import os.path
import matplotlib.cm as cm
import os
from scipy import signal
from opensbli.postprocess.plot_functions import *

from pyspod.spod.standard  import Standard  as spod_standard
from pyspod.spod.streaming import Streaming as spod_streaming
import pyspod.spod.utils     as utils_spod
import pyspod.utils.weights  as utils_weights
import pyspod.utils.errors   as utils_errors
import pyspod.utils.io       as utils_io
import pyspod.utils.postproc as post

plt.style.use('classic')

## -------------------------------------------------------------------
## initialize MPI
## -------------------------------------------------------------------
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.rank
except:
    comm = None
    rank = 0
## -------------------------------------------------------------------

# select number of snapshots for SPOD
# change this to 256 when re-run. Change also the n_dft to 128 later in the input file
it_start = -200

save_dir = './simulation_plots/'

# Extract data
def atoi(text):
        return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def plot_PSD(t,q,f,psd,ts,qs,fs,psds,peaks):

  # plot probe fluctuation history
  plt.figure()
  plt.subplot(2,1,1)
  plt.plot(t, q ,'b-',alpha=0.5,label='Probe sampling')
  plt.plot(ts,qs,'bo',          label='POD snapshot sampling')
  plt.xlabel("Time")
  plt.grid(True)
  # plt.legend()
  
  # plot the probe PSD distributions
  plt.subplot(2,1,2)
  plt.loglog(f, psd,               'b-',alpha=0.5,label='Probe sampling')
  plt.loglog(fs,psds,              'b-',          label='POD snapshot sampling')
  plt.loglog(fs[peaks],psds[peaks],'ro',          label='Selected POD frequencies')
  plt.xlabel("Frequency")
  plt.ylabel("PSD")
  # skips the first mode, whose amplitude is corrupted by the hanning function
  plt.xlim([2/(t[-1]-t[0]),0.5/(t[1]-t[0])])
  plt.grid(True)
  # plt.legend()
  plt.savefig('./simulation_plots/frequency-peaks_plot.pdf', bbox_inches='tight')

def calculate_PSD(delta_iteration_snapshots):

  probe = np.loadtxt("cylinder_probes.log", dtype='f', delimiter=',', skiprows=1)
  delta_iteration_probe = probe[1,0] - probe[0,0]
  skip  = int(delta_iteration_snapshots / delta_iteration_probe) 
  probe = probe[it_start*skip:,:]
  #----------------------------------------------------------------------#
  time    = probe[::,1] - probe[0,1] # iteration x dt_code => nondimensional
  q_prime = probe[::,2] - np.mean(probe[::,2])
  fs      = 1/(time[1]-time[0])
  f, PSD  = signal.periodogram(q_prime[0:-1], fs, window='hanning', scaling='spectrum')
  PSD     = 2 * (PSD / fs) / np.sqrt(2)
  #----------------------------------------------------------------------#
  time_s     = probe[::skip,1] - probe[0,1] # iteration x dt_code => nondimensional
  q_s_prime  = probe[::skip,2] - np.mean(probe[::skip,2])
  fs_s       = 1/(time_s[1]-time_s[0])
  f_s, PSD_s = signal.periodogram(q_s_prime[0:-1], fs_s, window='hanning', scaling='spectrum')
  PSD_s      = 2 * (PSD_s / fs_s) / np.sqrt(2)
  peaks,_    = signal.find_peaks(PSD_s,height=1E-5) 
  #----------------------------------------------------------------------#
  plot_PSD(time,q_prime,f,PSD,time_s,q_s_prime,f_s,PSD_s,peaks)  

  return time_s[1]-time_s[0], f_s[peaks]

def plot_mode(x,y,mode,f_idx,f,nvar):

    if nvar == "rho":
        label = "Density"
    elif nvar == "vortz":
        label = "Vorticity"
    elif nvar == "p":
        label = 'Pressure'
    elif nvar == 'u':
        label = 'u Velocity'
    elif nvar == 'v':
        label = 'v Velocity'

    # just for visualization purposes, as vorticity is not transferred at the boundaries
    if nvar == "vortz":
        mode[:,-1] = mode[:,-2]

    n_levels = 256
    min_val = np.min(np.real(mode))
    max_val  = np.max(np.real(mode))
    if nvar == 'v':
        levels0  = np.linspace(min_val, -min_val, n_levels)
    else:
        levels0  = np.linspace(-max_val, max_val, n_levels)

    ax, fig = plt.subplots()
    plt.contourf(x, y, mode, levels=levels0, cmap=cm.seismic)
    plt.gca().set_aspect('equal')
    # plt.title("Real part; Frequency =" + str(f))
    plt.xlabel(r"$x$", fontsize=18)
    plt.ylabel(r"$y$", fontsize=18)
    plt.xlim([-1,15])
    plt.ylim([-3.5,3.5])
    plt.plot(x[0,:], y[0,:], color='k', lw=1)
    ax = plt.gca()
    ax.text(0.025, 0.875, "%s" % (label),
    verticalalignment='bottom', horizontalalignment='left',
    transform=ax.transAxes,
    color='black', fontsize=14)
    # Save data
    np.save('./simulation_plots/x_OpenSBLI', x)
    np.save('./simulation_plots/y_OpenSBLI', y)
    np.save('./simulation_plots/' + nvar, mode)
    # Add GSA to some modes
    GSA_path = '/data/D/DA201G24/time-capsule/CYL-GSA-MODE/'
    if nvar == 'u' or nvar == 'v':
        x, z = np.load(GSA_path + 'GSA_x.npy'), np.load(GSA_path + 'GSA_z.npy')

        if nvar == 'u':
            GSA_mode = np.load(GSA_path + 'GSA_u.npy')
        elif nvar == 'v':
            GSA_mode = np.load(GSA_path + 'GSA_w.npy')
        plt.contour(x, z, GSA_mode, color='k')

    plt.savefig('./simulation_plots/POD_mode_'+str(nvar)+'_'+str(f_idx)+'.png', dpi=300, bbox_inches='tight')
    plt.clf()

class plotFunctions(object):
    def __init__(self):
        return

    def read_file(self, fname):
        f = h5py.File(fname, 'r')
        group = f["opensbliblock00"]
        return f, group

    def read_dataset(self, group, dataset):
        d_m = group["%s" % (dataset)].attrs['d_m']
        size = group["%s" % (dataset)].shape
        read_start = [abs(d) for d in d_m]
        read_end = [s-abs(d) for d, s in zip(d_m, size)] # include the interface line in the plot
        if len(read_end) == 2:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1]+1]
        else:
            raise NotImplementedError("")
        return read_data


class Plot(plotFunctions):
    def __init__(self):
        self.Minf = 0.1
        return

    def extract_coordinates(self):
        fname = 'data.h5'
        f, group1 = self.read_file(fname)
        x = self.read_dataset(group1, "x0_B0")
        y = self.read_dataset(group1, "x1_B0")
        return x, y

    def extract_flow_variables(self, group, Minf, gama):
        rho   = self.read_dataset(group, "rho_B0")
        rhou  = self.read_dataset(group, "rhou0_B0")
        rhov  = self.read_dataset(group, "rhou1_B0")
        rhoE  = self.read_dataset(group, "rhoE_B0")
        vortz = self.read_dataset(group, "wz_B0")
        u = rhou/rho
        v = rhov/rho
        p = (rhoE - 0.5*rho*(u**2 + v**2))*(gama-1)
        return rho, u, v, p, vortz

    def get_OpenSBLI_files(self):
        # Extract data
        PP = OpenSBLIPreProcess()
        fnames, iters = PP.find_files('./')
        fnames.sort(key=natural_keys)
        fnames = fnames[it_start:]
        iters  = iters[it_start:]
        dt = 0.0001
        times = [dt*float(x) for x in iters]
        iters = [int(x) for x in iters]
        return fnames, iters, times

    def spod(self, fname, n_levels,nvar):
        # Get 2D OpenSBLI output snapshot file names
        if not restart:
            fnames, iters, times = self.get_OpenSBLI_files()
            # Get x, y coordinates
            self.x, self.y = self.extract_coordinates()
            Nx, Ny = self.x.shape[-1], self.x.shape[0]
            nsamples = len(iters)
            OpenSBLI_data = []
            data = np.zeros((nsamples, Ny, Nx))

            for i, fname in enumerate(fnames):
                f, group1 = self.read_file(fname)
                Minf, gama = f['Minf'][()][0], f['gama'][()][0]
                rho, u, v, p, vortz = self.extract_flow_variables(group1, Minf, gama)
                if nvar == "rho":
                  data[i,:,:] = rho
                elif nvar == "u":
                  data[i,:,:] = u
                elif nvar == "v":
                  data[i,:,:] = v
                elif nvar == "p":
                  data[i,:,:] = p
                elif nvar == "vortz":
                  data[i,:,:] = vortz

            f.close()

            names = nvar
            # perform PSD
            delta_iterations_snapshots = iters[1]-iters[0]
            dt, f_peaks = calculate_PSD(delta_iterations_snapshots)
            # start SPOD
            config_file = 'input_cylinder.yaml'
            params = utils_io.read_config(config_file)
            params['time_step'] = dt

            # Call SPOD
            ## -------------------------------------------------------------------
            ## compute spod modes and check orthogonality
            ## -------------------------------------------------------------------
            standard  = spod_standard (params=params, comm=comm)
            streaming = spod_streaming(params=params, comm=comm)
            # if not restart:
            spod = standard.fit(data_list=data)
            # spod = streaming.fit(data_list=data)
            results_dir = spod.savedir_sim
            flag, ortho = utils_spod.check_orthogonality(
                results_dir=results_dir, mode_idx1=[1],
                mode_idx2=[0], freq_idx=[5], dtype='double',
                comm=comm)
            print(f'flag = {flag},  ortho = {ortho}')
        ## -------------------------------------------------------------------
            ## compute coefficients
            ## -------------------------------------------------------------------
            file_coeffs, coeffs_dir = utils_spod.compute_coeffs_op(
                data=data, results_dir=results_dir, comm=comm)
            ## -------------------------------------------------------------------
            ## only rank 0
            if rank == 0:
                ## identify frequency of interest and plot corresponding POD modes
                for i in f_peaks:
                  f, f_idx = spod.find_nearest_freq(freq_req=i, freq=spod.freq)
                  modes = post.get_modes_at_freq(results_path=results_dir, freq_idx=f_idx)
                  plot_mode(self.x,self.y,modes[:,:,0,0],f_idx,f,nvar)
        else:
            print("Loading previous data for variable: {}.".format(nvar))
            # Load previous data
            x, y = np.load('./simulation_plots/x_OpenSBLI.npy'), np.load('./simulation_plots/y_OpenSBLI.npy')
            dset = np.load('./simulation_plots/%s.npy' % nvar)
            f_idx, f = 3, 1
            plot_mode(x, y, dset,f_idx,f,nvar)

fname = "opensbli_output.h5"
n_contour_levels = 256
directory = './simulation_plots/'

restart = True

if not os.path.exists(directory):
    os.makedirs(directory)

KP = Plot()
variables=["rho","u","v","p","vortz"]
for nvar in variables:
  KP.spod(fname, n_contour_levels,nvar)
