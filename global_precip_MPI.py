# Calculate global precipitation due to GLD data
# (in parallel, on nansen)

import os
import commands
import datetime
import random
import matplotlib
matplotlib.use("Agg")

import time
import numpy as np 
import itertools
import subprocess
from partition import partition

import pickle

from mpi4py import MPI 
from calc_global_precip import calc_global_precip
from precip_model import precip_model
from GLD_file_tools import GLD_file_tools
from plotting import plot_flux_basemap


import matplotlib.pyplot as plt


# Grids for output arrays
map_lats = np.arange(-90, 90, 1)
map_lons = np.arange(-180, 180, 1)

# Grids for the precalculated precipitation model
in_lat_grid  = np.arange( 15, 55, 1)
out_lat_grid = np.arange(-60, 60, 1);
out_lon_grid = np.arange(-10, 10, 1);
out_time_grid= np.arange(  0, 20, 0.1);

window_time = 60  # Window length for time integration
task_spacing = datetime.timedelta(minutes=30) # Spacing between tasks (synoptic period)


# Dates to search thru
sim_start = datetime.datetime(2015, 02, 1, 0, 0, 0)
sim_stop  = datetime.datetime(2015, 03, 1, 0, 0, 0)

root_dir ='/shared/users/asousa/WIPP/global_precip_agu2016/' 
# db_file  = os.path.join(root_dir,'db_agu2016_kp0_full.pkl')
db_file  = os.path.join(root_dir,'banded_4_mode2_diff/db_band_3_10Mev_100Mev.pkl')

out_dir = os.path.join(root_dir, 'banded_4_mode2_diff')
fig_dir = os.path.join(out_dir, 'figures_band_3')
dump_dir= os.path.join(out_dir, 'saves_band_3')



# Initialize MPI:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")

# GLD mount:
gld = GLD_file_tools(os.path.join(os.path.expanduser("~"),'GLD_mount'), prefix='GLD')


# --------------------- Prep jobs to scatter ------------------------
if rank==0:
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)
    if not os.path.exists(dump_dir):
        os.mkdir(dump_dir)

    print "available nodes: ",comm.Get_size()
    print "Setting up parallel runs..."

    run_starttime = time.time()
    print "start time: %s"%run_starttime
    sss = sim_start

    # generate tasklist
    tasklist = []
    while sss<= sim_stop:
        tasklist.append(sss)
        sss += task_spacing

    nTasks = 1.0*len(tasklist)
    nProcs = 1.0*comm.Get_size()
    nSteps = np.ceil(nTasks/nProcs).astype(int)

    chunks = partition(tasklist, nProcs)

    print "%d total chunks"%(len(chunks))
else:
    chunks = None
    tasklist = None

# -------------------- Prep model bits -----------------------------
if rank==0:
    print "Setting up model stuff..."
    p = precip_model(database=db_file, cumsum=True)
    p.sc.I0 = -10000 # Setting here because it's not in the consts file anymore
    p.precalculate_gridded_values(in_lat_grid, out_lat_grid, out_lon_grid, out_time_grid)
else:
    p = None

# Broadcast from root to nodes:
p = comm.bcast(p, root=0)
chunks = comm.bcast(chunks, root=0)

# Send iiiiit
chunk = comm.scatter(chunks)
print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunk))

# vv----------------- Do stuff here --------------------------------
for in_time in chunk:
    print "[%s/%d] starting at %s"%(host, rank, in_time)
    try:
        flux, flashes = calc_global_precip(p, gld, in_time, window_time, map_lats, map_lons)

        if len(flashes) > 0:
            if len(flashes) < 5000:    # There's a bug in the GLD file loader, I guess

                # fig = plot_flux_basemap(flux/window_time, map_lats, map_lons, flashes,
                #                         plottime=in_time, logscale=True, clims=[-3,5],
                #                         num_contours=20, mode='counts')


                # fig = plot_flux_basemap(flux/window_time, map_lats, map_lons, flashes,
                #                         plottime=in_time, logscale=True, clims=[-4,1],
                #                         num_contours=20, mode='energy')

                # plt.savefig('%s/%s.png'%(fig_dir, in_time.strftime('%Y-%m-%d_%H-%M-%S')),bb_inches='tight')
                # plt.close('all')

                file = open(os.path.join(dump_dir,'%s.pkl'%in_time),'wb')
                pickle.dump([[in_time.isoformat()],flux],file)
                file.close()
            else:
                print "[%s/%d] Too many flashes -- skipping at %s"%(host, rank, in_time)


    except:
        print "[%s/%d] Something went wrong at %s"%(host, rank, in_time)


# ^^----------------- Do stuff here --------------------------------
