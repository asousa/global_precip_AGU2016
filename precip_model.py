# precip_model.py
# AGU2016 version -- interpolating over inlat, outlat, outlon, time
import numpy as np
from scipy import interpolate
import pickle




class precip_model(object):
    def __init__(self,database="database.pkl", cumsum = False):

        self.D2R = np.pi/180.0

        with open(database,'r') as file:
            self.db = pickle.load(file)

        self.in_lats = self.db['in_lats']
        self.out_lats= self.db['out_lats']
        self.out_lons= self.db['out_lons']
        self.t = self.db['t']
        self.sc = self.db['consts']


        # Hack to get it to interpolate when there's only one longitude:
        # (sorry)
        if len(self.out_lons) == 1:
            self.out_lons = np.arange(-10, 11, 1)
            weights = np.exp(-1.0*np.abs(self.out_lons)/2.1)
            print weights

            self.db['N'] = self.db['N']*weights[np.newaxis, np.newaxis, :, np.newaxis]
            self.db['S'] = self.db['S']*weights[np.newaxis, np.newaxis, :, np.newaxis]
            # self.db['N'] = np.tile(self.db['N'], [1,1,2,1])
            # self.db['S'] = np.tile(self.db['S'], [1,1,2,1])
            # self.out_lons = np.array([-0.5, 0.5]) + self.out_lons[0]

        if cumsum:
            self.N = np.cumsum(np.maximum(0, self.db['N']),axis=3)*self.sc.T_STEP
            self.S = np.cumsum(np.maximum(0, self.db['S']),axis=3)*self.sc.T_STEP
        else:
            self.N = np.maximum(0, self.db['N'])
            self.S = np.maximum(0, self.db['S'])

        # Interpolating objects:
        self.N_interp = interpolate.RegularGridInterpolator((self.in_lats, self.out_lats, self.out_lons, self.t), self.N, fill_value=0, bounds_error=False)
        self.S_interp = interpolate.RegularGridInterpolator((self.in_lats, self.out_lats, self.out_lons, self.t), self.S, fill_value=0, bounds_error=False)



        # Initialize any other parameters we might store:
        self.precalculated = None
        self.pc_in_lats = None
        self.pc_out_lats = None
        self.pc_t = None



    def get_precip_at(self, in_lat, out_lat, out_lon, t):
        ''' in_lat:  Flash latitude (degrees)
            out_lat: Satellite latitude (degrees)
            out_lon: Satellite longitude (degrees)
            t:       Time elapsed from flash (seconds)
            '''
        # Model is symmetric around northern / southern hemispheres (mag. dipole coordinates):
        # If in = N, out = N  --> Northern hemisphere
        #    in = N, out = S  --> Southern hemisphere
        #    in = S, out = N  --> Southern hemisphere
        #    in = S, out = S  --> Northern hemisphere

        # assert len(in_lat) == len(out_lat) == len(t) == len(out_lon), "Length mismatch!"

        use_southern_hemi = np.array(((in_lat > 0) ^ (out_lat > 0)).ravel())
        
        keys = np.array([np.abs(in_lat.ravel()),np.abs(out_lat.ravel()), np.abs(out_lon.ravel()), t.ravel()]).T

        keys_N = keys[~use_southern_hemi,:]
        keys_S = keys[ use_southern_hemi,:]

        out_data = np.zeros(len(keys))
        out_data[ use_southern_hemi] = self.S_interp(keys_S)
        out_data[~use_southern_hemi] = self.N_interp(keys_N)
        # out_data = self.N_interp(keys)

        return out_data

    def get_multiple_precip_at(self, in_lats, out_lats, out_lons, t):
        '''A vectorized version of get_precip_at(). This should be faster!
            in_lats, out_lats, t are all vectors.
            Returns: A numpy array of dimension [in_lats x out_lats x out_lons x t]
        '''

        tx, ty, tw, tz  = np.meshgrid(in_lats, out_lats, out_lons, t)

        keys =  np.array([np.abs(tx.ravel()),np.abs(ty.ravel()), np.abs(tw.ravel()), tz.ravel()]).T
        # print np.shape(tx), np.shape(ty), np.shape(tz)
        # print keys
        # keys = cartesian([in_lats, out_lats, t])

        
        # Model is symmetric around northern / southern hemispheres (mag. dipole coordinates):
        # If in = N, out = N  --> Northern hemisphere
        #    in = N, out = S  --> Southern hemisphere
        #    in = S, out = N  --> Southern hemisphere
        #    in = S, out = S  --> Northern hemisphere
        use_southern_hemi = np.array(((tx > 0) ^ (ty > 0)).ravel())
        keys = np.abs(keys)
        # use_southern_hemi = np.array(((keys[:,0] > 0) ^ (keys[:,1] > 0)).ravel())
        # keys = np.abs(keys)

        keys_N = keys[~use_southern_hemi,:]
        keys_S = keys[ use_southern_hemi,:]

        out_data = np.zeros(len(keys))


        out_data[ use_southern_hemi] = np.maximum(0, self.S_interp(keys_S))
        out_data[~use_southern_hemi] = np.maximum(0, self.N_interp(keys_N))

        # out_data = self.get_precip_at(tx, ty, tz)


        # Indexes need to be swapped because of how meshgrid arranges things. Works tho!
        # return out_data.reshape(len(in_lats),len(out_lats),len(t))
        return out_data.reshape(len(out_lats), len(in_lats),  len(out_lons), len(t)).swapaxes(0,1)

        # return out_data.reshape(len(out_lats),len(in_lats),len(t)).swapaxes(0,1)

    def precalculate_gridded_values(self, in_lats, out_lats, out_lons, t):
        print "Precalculating..."
        
        self.precalculated = self.get_multiple_precip_at(in_lats, out_lats, out_lons, t)
        self.pc_in_lats = in_lats
        self.pc_out_lats = out_lats
        self.pc_out_lons = out_lons
        self.pc_t = t




    def power_scaling(self, I0):

        return np.sqrt(np.abs(I0/self.sc.I0))



if __name__=='__main__':
    p = precip_model(cumsum=True)
    print p.get_multiple_precip_at(np.array([30, -30]),
                          np.array([40, 50]),
                          np.array([0,  0]),
                          np.array([0.1, 20]))

