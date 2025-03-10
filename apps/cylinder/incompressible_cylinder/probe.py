# import numpy 
import numpy as np
# import matplotlib for plotting figures
import matplotlib.pyplot as plt
# import scipy for the PSD calculation
from scipy import signal

np.seterr(divide="ignore")

#----------------------------------------------------------------------#
# read the probe file, made of float values
probe = np.loadtxt("block0_cylinder_probes.log", dtype='f', delimiter=',', skiprows=1)

#----------------------------------------------------------------------#
# construct the iterations array
start_index = -10000
iterations  = probe[start_index:,0]-probe[start_index,0] # iterations
time        = probe[start_index:,1]-probe[start_index,1] # iteration x dt_code => nondimensional

end_index    = 3500
time_intrans = probe[:end_index,1]-probe[0,1]

for i in range(2,3):

  q1       = probe[start_index:,i]
  q1_mean  = np.mean(q1)
  q1_prime = q1 - q1_mean

  q1_intrans = np.abs(probe[:end_index,i]-probe[0,i])

  #----------------------------------------------------------------------#
  # get the sampling frequency
  fs = 1/(time[1]-time[0])
  
  #calculate the PSD for the lift coefficients
  f, PSD = signal.periodogram(q1_prime[0:-1], fs, window='hamming', scaling='spectrum') # window=hanning/blackman
  PSD = 2 * (PSD / fs) / np.sqrt(2)
  max_idx = np.argmax(PSD)

  print("Peak Frequency",f[max_idx])

  peaks,_ = signal.find_peaks(np.log(q1_intrans))

#----------------------------------------------------------------------#
# plot the drag and the lift coefficients
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(time,q1)
plt.xlabel("Time")
plt.grid(True)

# plot the lift PSD distributions
plt.subplot(2,1,2)
plt.loglog(f,PSD)
plt.xlabel("Frequency")
plt.ylabel("PSD")
plt.grid(True)
plt.savefig('./simulation_plots/frequency_plot.pdf', bbox_inches='tight')
  
plt.figure(2)
plt.plot(time_intrans,np.log(q1_intrans))
plt.plot(time_intrans[peaks[5:]],np.log(q1_intrans[peaks[5:]]))

growth_rate = (np.log(q1_intrans[peaks[-5]])-np.log(q1_intrans[peaks[-10]]))/(time_intrans[peaks[-5]]-time_intrans[peaks[-10]])
print("Growth rate = ", growth_rate)

# show the plot
plt.savefig('./simulation_plots/growth_rate_plot.pdf', bbox_inches='tight')
