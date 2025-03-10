# import numpy 
import numpy as np
# import matplotlib for plotting figures
import matplotlib.pyplot as plt

schemes = ['Divergence', 'Blaisdell', 'Jameson', 'Kok', 'Feiereisen', 'KGP', 'KEEP']
fname = 'block0_TGV.log'

# 256^3 runtimes for reference
run_times = [53.7508, 59.4306, 63.2548, 63.3555, 64.715, 68.4992, 67.3297]
run_times = [x / run_times[0] for x in run_times]
colors = ['b', 'g', 'r', 'k', 'm', 'c', 'y'][::-1]

# plot the input_data and forces iteration evolution 
fontsize = 16
plt.rcParams.update({'font.size': fontsize})
plt.rcParams.update({"text.usetex": True,"font.family": "serif","font.serif": ["Palatino"]})

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6.4*2, 4.8))

ax0, ax1 = axes[0], axes[1]

for i, scheme in enumerate(schemes):
    print("Plotting case: {}".format(scheme))

    # read the input_data file
    input_data  = np.genfromtxt('./%s/' % scheme + '/' + fname, delimiter=',', skip_header=1)

    # construct the arrays
    iterations = input_data[:,0]
    times = input_data[:,1]
    KE = input_data[:,2]

    ax0.plot(times[2:], KE[2:]/KE[2], label=scheme, color=colors[i])
    ax0.set_ylim([0, 1.5])
    ax0.set_xlabel(r"Time")
    ax0.set_ylabel(r"$KE / KE_0$")
    ax0.grid(True)
    ax0.legend(ncol=2, loc='best')

ax1.grid(zorder=10000)
ax1.bar(schemes, run_times, label=schemes, color=colors, zorder=1000)
plt.xticks(rotation='vertical')
ax1.set_ylabel('Computational cost')

# show the plot
# plt.show()
fig.tight_layout()
plt.savefig('SplitSchemes_TGV.png', dpi=400)
plt.clf()
