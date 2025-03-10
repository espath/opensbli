# import numpy 
import numpy as np
# import matplotlib for plotting figures
import matplotlib.pyplot as plt
# import scipy for the PSD calculation
from scipy import signal
import os

np.seterr(divide="ignore")

#----------------------------------------------------------------------#
# read the probe file, made of float values
probe = np.loadtxt("../block2_airfoil_output.log", dtype='f', delimiter=',', skiprows=1)
#----------------------------------------------------------------------#
# pick only second probe location
i = 5

# entire signal
q_prime_all  = probe[:,i] - np.mean(probe[:,i]) 
time_all     = probe[:,1] - probe[0,1] # iteration x dt_code => nondimensional
print("Total collection time: {}".format(time_all[-1]))

# saturated signal only
start_index  = -3900
start_index = 0
time_sat     = probe[start_index:,1] - probe[start_index,1] # iteration x dt_code => nondimensional
q_prime_sat  = probe[start_index:,i] - np.mean(probe[start_index:,i])

# # initial transient only
# end_index    = 3500
# time_intrans = probe[:end_index,1] - probe[0,1]
# q_intrans    = np.abs(probe[:end_index,i] - probe[0,i])

#----------------------------------------------------------------------#
# PSD on saturated signal #

# get the sampling frequency
fs = 1/(time_sat[1]-time_sat[0])

#calculate the PSD for the lift coefficients
f, PSD = signal.periodogram(q_prime_sat[0:-1], fs, window='hann', scaling='spectrum')
PSD = 2 * (PSD / fs) / np.sqrt(2)

#----------------------------------------------------------------------#
# find peaks in the initial transient #

# peaks,_ = signal.find_peaks(np.log(q_intrans))

#----------------------------------------------------------------------#
# Get growth rate and frequency #

# growth_rate = (np.log(q_intrans[peaks[-5]])-np.log(q_intrans[peaks[-10]]))/(time_intrans[peaks[-5]]-time_intrans[peaks[-10]])
# print("Growth rate = ", growth_rate)

max_idx = np.argmax(PSD)
print("Peak Frequency",f[max_idx])

#----------------------------------------------------------------------#
# make plot #
linewidth = 2
# plot the entire time evolution
plot_range = 1
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(time_all[-plot_range:] - time_all[-plot_range],q_prime_all[-plot_range:],color="black",linewidth=linewidth, label='$10\%$ of the data collection period')
# plt.axvline(x = time_all[start_index], color='r', linestyle='dashed',linewidth=linewidth,label="Sample collection for PSD")
plt.xlabel(r"$t - t_0$")
plt.ylabel(r"$v$")
# plt.xlim([0,)
plt.legend(loc="upper right")
plt.grid(True)

# plot the initial transient and linear growth rate
# plt.subplot(2,2,3)
# plt.plot(time_all,np.log(q_prime_all),color="black",linewidth=linewidth)
# plt.plot(time_intrans[peaks[5:]],np.log(q_intrans[peaks[5:]]),color="red", linestyle='dashed',linewidth=linewidth,label="Growth rate")
# plt.xlabel(r"$t$")
# plt.ylabel(r"$log\left|v\prime - v\prime\left(t=0\right)\right|$")
# plt.xlim([0,time_all[-1]])
# plt.legend(loc="lower right")
# plt.grid(True)

# plot the PSD distribution
plt.subplot(2,1,2)
# Add LJ reference line at f = 2.52
plt.axvline(x=2.52, color='red', linewidth=linewidth, ls='--', label='Reference, Jones (2008)')
plt.loglog(f,PSD,color="black",linewidth=linewidth, label='OpenSBLI')

plt.legend()
plt.xlabel("Frequency")
plt.ylabel("PSD")
print(2*f[1])
plt.xlim([2*f[1],f[-1]])
# plt.xlim([0.2, f[-1]])
plt.grid(True)

plt.tight_layout()

directory = './simulation_plots/'
if not os.path.exists(directory):
    os.makedirs(directory)

# show the plot
plt.savefig('./simulation_plots/LJ_growth-rate_frequency_plot.pdf', bbox_inches='tight')
