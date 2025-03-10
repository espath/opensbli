import numpy
import matplotlib.pyplot as plt
import h5py
import os.path
import os

plt.style.use('classic')


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
        read_end = [s-abs(d) for d, s in zip(d_m, size)]
        if len(read_end) == 1:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0]]
        elif len(read_end) == 2:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1]]
        elif len(read_end) == 3:
            read_data = group["%s" % (dataset)][read_start[0]:read_end[0], read_start[1]:read_end[1], read_start[2]:read_end[2]]
        else:
            raise NotImplementedError("")
        return read_data


class Plot(plotFunctions):
    def __init__(self):
        return

    def line_graphs(self, x, variable, name, xref, ref, xref2=None, ref2=None):
        if ref is not 0:
            plt.plot(xref, ref, color='b', label='OldLLF')
            #plt.plot(xref2, ref2, color='k', label='Fine_Mesh')
        if name == "rho":
            plt.semilogy(x, variable, color='r', label='NewLLF')
            plt.ylim([0.0008, 1.1])
        elif name == "P":
            plt.semilogy(x, variable, color='r', label='NewLLF')
            # plt.ylim([0.0008, 1.1])           
        else:
            plt.plot(x, variable, color='r', label='NewLLF')
        plt.xlabel(r'$x_0$', fontsize=20)
        plt.ylabel(r'$%s$' % name, fontsize=20)
        plt.legend(loc="best")
        plt.savefig(directory + "output_%s.pdf" % name, bbox_inches='tight')
        plt.clf()
        return

    def extract_flow_variables(self, group):
        rho = self.read_dataset(group, "rho_B0")
        rhou = self.read_dataset(group, "rhou0_B0")
        rhoE = self.read_dataset(group, "rhoE_B0")
        u = rhou/rho
        p = (0.4)*(rhoE - 0.5*(u**2)*rho)
        #kappa = self.read_dataset(group, "kappa_B0")
        #q0 = self.read_dataset(group, "q0_B0")
        #q1 = self.read_dataset(group, "q1_B0")
        #q2 = self.read_dataset(group, "q2_B0")
        return rho, u, rhoE, p

    def save_data(self, fname, x, rho, u, P):
        numpy.savetxt(fname, numpy.c_[x, rho, u, P])
        return

    def main_plot(self, fname, n_levels):
        f, group1 = self.read_file(fname)
        rho, u, rhoE, p  = self.extract_flow_variables(group1)
        variables = [rho, u, p]
        names = ["rho", "u", "P"]
        x = numpy.linspace(0, 10, rho.size)
        save = True

        if save:
                self.save_data('NewLF_WENO7Z.txt', x, rho, u, p)

        # Load reference data
        #data = numpy.loadtxt('OldLF_WENO7Z.txt')
        #xref, rhoref, uref, Pref = data[:,0], data[:,1], data[:,2], data[:,3]
        #ref = [rhoref, uref, Pref, 0, 0, 0, 0]

        # Load reference data
        data = numpy.loadtxt('TENO6_reference.txt')
        xref2, rhoref2, uref2, Pref2 = data[:,0], data[:,1], data[:,2], data[:,3]
        ref2 = [rhoref2, uref2, Pref2, 0, 0, 0, 0]

        for i, (var, name) in enumerate(zip(variables, names)):
            self.line_graphs(x, var, name, xref2, ref2[i])#, xref2, ref2[i])
            f.close()


fname = "opensbli_output.h5"
n_contour_levels = 25
directory = './simulation_plots/'

if not os.path.exists(directory):
    os.makedirs(directory)

KP = Plot()
KP.main_plot(fname, n_contour_levels)
