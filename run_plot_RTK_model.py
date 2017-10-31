import pandas as pd
import plotting_functions as pf
import os
import matplotlib.pyplot as plt
from decimal import Decimal
outdir = "/Users/priyaveeraraghavan/Documents/SorgerRotation/energy_model/RTK_model/ode_plots"
h5_storage = 'RTK_simulation_output.h5'
store = pd.HDFStore(h5_storage)

ode_sim_df = store['ode_sim_out']
"""
if not os.path.exists(outdir):
    os.mkdir(outdir)
else:
    print "Had to make a new directory!!!"
    os.mkdir(outdir.strip("/") + "_new")
    outdir = outdir.strip("/") + "_new"
"""

for a_b_tuple, a_b_subset in ode_sim_df.groupby(["a", "b"]):
    a_b_tuple = tuple([str(x).replace(".", "_") for x in a_b_tuple])
    pf.plot_3d(a_b_subset, "conc_RTK_0", "conc_GF_0", "time")
    plt.savefig(os.path.join(outdir, "a_%s_b_%s_time.pdf" % a_b_tuple), bbox_inches='tight')
    plt.close()

max_time = max(ode_sim_df.time)
ss_df = ode_sim_df[ode_sim_df.time == max_time]
for species in set(ode_sim_df.species):
    pf.plot_heatmap_array(ss_df, species, "conc_RTK_0", "conc_GF_0", "a", "b")
    plt.savefig(os.path.join(outdir, "%s_heatmap.png" % species), bbox_inches='tight')
    plt.close()

# plot the distribution of activated receptors for cells
# this will be a proxy for distribution over many cells, due to time invariance
for a_b_tuple, a_b_subset in sim_df.groupby(["a", "b"]):
    a_b_tuple = tuple([str(x).replace(".", "_") for x in a_b_tuple])
    ss_timepoint = store['ss_timepoints']
    pf.plot_steady_state_distribution_generator(a_b_subset, ss_timepoint)
    plt.savefig(os.path.join(outdir, "%s_activated_histogram.png" % species), bbox_inches='tight')
    plt.close()
