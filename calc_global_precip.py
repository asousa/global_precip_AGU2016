import numpy as np
from GLD_file_tools import GLD_file_tools
from precip_model import precip_model
import datetime as dt
from coordinate_structure import transform_coords
from daynite_scaler import daynite_scaler
import matplotlib.pyplot as plt

# from mpl_toolkits.basemap.solar import daynight_grid
# from mpl_toolkits.basemap.solar import daynight_terminator
# ---------------------------------------------------------
# calc_global_precip.py
# ---------------------------------------------------------
#
# Calculate a heatmap of global precipitation:
#   p:          precip_model object to use
#   gld:        GLD_file_tools object to retrieve flash data
#   in_time:    time to plot at (datetime object)
#   window_time:time previous to look for flashes at (timedelta object)
#   grid_lats:  Output latitudes (geomagnetic)
#   grid_lons:  Output longitudes (geomagnetic)
#
#   V1.0 ~ 6.2016 APS

def calc_global_precip(p, gld, in_time, window_time, grid_lats, grid_lons):

    # Instantiate objects:
    # p = precip_model(database="db_test.pkl", cumsum=True)
    # gld = GLD_file_tools('GLD_mount',prefix='GLD')

    # How far back to load flashes from:
    # (window time plus the full length of the model - this way we 
    #  account for all possible flux)
    lookback_time = dt.timedelta(seconds=p.t[-1] + window_time)
    # lookback_time = dt.timedelta(seconds=window_time)
    in_lat_grid  = np.arange(-60, 60, step=0.1)


    # print "Precalculating..."
    # p.precalculate_gridded_values(in_lat_grid, grid_lats, p.t)

    map_lat_inds = nearest_index(grid_lats, p.pc_out_lats)

    # Get Day-night terminator map:
    dn_s = daynite_scaler(in_time)

    # dn_lon = dn_lon_grid[0,:]
    # dn_lat = dn_lat_grid[:,0]
    print "starting at %s"%in_time

    lat_ind = 7
    lon_ind = 8
    mag_ind = 9

    flux = np.zeros([len(grid_lats), len(grid_lons)], dtype='d')

    # print "Loading flashes..."
    # (Loads flashes between (in_time - lookback_time) and (in_time), including current second)
    flashes, flash_times = gld.load_flashes(in_time, lookback_time)

    if flashes is not None:
        flashes = flashes[:,(lat_ind, lon_ind, mag_ind, mag_ind)]
        flash_coords = transform_coords(flashes[:,0], flashes[:,1], np.zeros_like(flashes[:,0]),
                                        'geographic', 'geomagnetic')
        flashes[:,:2] = flash_coords[:,:2]
        flashes[:,3] = [(in_time - s).microseconds*1e-6 + (in_time - s).seconds for s in flash_times]


        # Mask out flashes outside the range of the interpolator:
        mask = (  (np.abs(flashes[:,0]) >= 15) 
                & (np.abs(flashes[:,0]) <= 55)) 
                # (flashes[:,2] > 50))



        print "%g flashes (post-filter)" % np.sum(mask)

        # masked_flashes = flashes[mask, :]
        flashes = flashes[mask, :]

        for ind, f in enumerate(flashes):
            t_end   = np.min([p.pc_t[-1], f[3]])   # Min to account for the interpolator returning 0 when outside of range
            t_start = np.max([0, t_end - window_time]) # np.max([0,f[3] - p.t[-1]])


            in_lat_ind = nearest_index(p.pc_in_lats, [f[0]])
            t_end_ind  = nearest_index(p.pc_t, [f[3]])
            t_start_ind = nearest_index(p.pc_t, [f[3] - window_time])

            # Get longitude indexes for output array
            map_lons = p.pc_out_lons + f[1]
            map_lons[map_lons < grid_lons[0]]  = grid_lons[-1] - map_lons[map_lons < grid_lons[0]] 
            map_lons[map_lons > grid_lons[-1]] -= grid_lons[-1]
            map_lon_inds = nearest_index(grid_lons, map_lons)

            # Get latitude indexes for output array
            # print np.shape(map_lat_inds)

            # print "flash coords: ", f
            # print "in_lat_ind:", in_lat_ind
            # print "t_end_ind:", t_end_ind
            # print "t_start_ind:", t_start_ind

            # print "Scalefactor: ", dn_s.scaling_factor_at(f[0], f[1])

            lv = (p.precalculated[in_lat_ind, :, :, t_end_ind] - p.precalculated[in_lat_ind,:, :, t_start_ind]).squeeze()

            power_scaling_factor = np.sqrt(np.abs(f[2]/(p.sc.I0/1000.))) 
            print "Flash: ", f[2], "pwr scaling: ", power_scaling_factor
            lv *= dn_s.scaling_factor_at(f[0], f[1])    # day-nite scaling
            lv *= power_scaling_factor                # Peak current scaling
            # # Freshly interpolated:
            # lv = (p.get_multiple_precip_at([f[0]], grid_lats, [t_end]) - p.get_multiple_precip_at([f[0]], grid_lats, [t_start])).squeeze()
            # print "flux shape: ", np.shape(flux)
            # print "flux_section shape:", np.shape(flux[map_lat_inds[:,np.newaxis],map_lon_inds[np.newaxis,:]])
            # print "new section shape:", np.shape(lv)

            # flux[:, map_lon_inds] += lv.T         
            flux[map_lat_inds[:,np.newaxis],map_lon_inds[np.newaxis,:]] += lv
            # flux += scalefactor*(np.outer(lv, daynite_vector))

    else:
        print "No flashes found at ", in_time


    print "finished %s"%in_time

    return flux/window_time, flashes





def nearest_index(grid, values):
    # Find closest index of a value in an array (i.e., quick quantize to grid value)
    idx = np.searchsorted(grid, values, side="left")
    idx = np.clip(idx, 0, len(grid) - 1)
    idx_l = np.clip(idx - 1, 0, len(grid) - 1)

    idx[abs(values - grid[idx_l]) < abs(values - grid[idx])] -= 1
    return idx