import re
import numpy as np
import pickle
import os
import itertools
from load_sim_constants import load_sim_constants
from load_phi_files import load_phi_files_latlon as load_phi_files

# Build a database of fluxes to interpolate over

def particle_flux(phi, sc):
    '''Convert a phi file into particle fluxes.'''
    ev2joule = (1.60217657)*1e-19 # Joules / ev
    joule2millierg = 10*1e10
    
    # Energy vector, in ev
    E_EXP = sc.E_EXP_BOT + np.linspace(1,sc.NUM_E,sc.NUM_E)*sc.DE_EXP
    E = np.power(10, E_EXP)
    
    #  Energy differential term dE, in kev
    DE = np.exp(np.linspace(1,sc.NUM_E,sc.NUM_E)*sc.DE_EXP/np.log10(np.e))
    DE = DE*sc.DE_EXP/np.log10(np.e)
    DE = DE*(1e-3)*np.power(10, sc.E_EXP_BOT + sc.DE_EXP/2.)
    
    
    # Particle flux: integrate (phi*dE) over each bin.
    # Output is (particles/sec) within each energy bin. Total particle flux ~ sum over bins in range.
    N = phi*DE[:, np.newaxis, np.newaxis]

#     totals = np.sum(phi, axis=1)
    return N

def energy_flux(phi, sc):    
    '''Convert a phi file into energy fluxes.'''
    ev2joule = (1.60217657)*1e-19 # Joules / ev
    joule2millierg = 10*1e10
    
    # Energy vector, in ev
    E_EXP = sc.E_EXP_BOT + np.linspace(1,sc.NUM_E,sc.NUM_E)*sc.DE_EXP
    E = np.power(10, E_EXP)
    
    #  Energy differential term dE, in kev
    DE = np.exp(np.linspace(1,sc.NUM_E,sc.NUM_E)*sc.DE_EXP/np.log10(np.e))
    DE = DE*sc.DE_EXP/np.log10(np.e)
    DE = DE*(1e-3)*np.power(10, sc.E_EXP_BOT + sc.DE_EXP/2.)
    
    

    # Output is (mErgs/sec) within each energy bin. Total particle flux ~ sum over bins in range.    
    # Energy flux: Integrate E*phi*dE
    Q = phi*(E*DE)[:,np.newaxis, np.newaxis]*ev2joule*joule2millierg

    return Q
    

def build_database(input_dir='outputs', output_filename='database.pkl'):
    # Load constants:
    sc = load_sim_constants(os.path.join(input_dir, 'consts.h'))



    p = re.compile("\d+")  #regular expression number finder magic
    inlats = sorted([int(p.findall(s)[0]) for s in os.listdir(input_dir) if s.startswith('in_')])
    print "input latitudes", inlats
    # inlats = [20, 30, 40]


    # output latitudes (assume constant for all input latitudes)
    outlats = sorted([int(p.findall(s)[0]) for s in os.listdir(os.path.join(input_dir,'in_%d'%inlats[0])) if s.startswith('lat_')])
    print "output latitudes", outlats

    # # output longitudes (assume constant for all in lats + out lats)
    outlons = sorted(np.unique([int(p.findall(s)[1]) for s in os.listdir(os.path.join(input_dir,'in_%d'%inlats[0], 'lat_%d'%outlats[0]))]))
    print "output longitudes", outlons

    # outlons = [0,2]


    obj = dict()
    obj['consts'] = sc
    obj['out_lats'] = outlats
    obj['in_lats'] = inlats
    obj['out_lons'] = outlons
    obj['t'] = np.linspace(sc.T_STEP,sc.T_MAX,sc.NUM_STEPS)



                #  input lats,  output lats,  output lons,   time
    N = np.zeros([len(inlats), len(outlats), len(outlons), sc.NUM_STEPS])
    S = np.zeros([len(inlats), len(outlats), len(outlons), sc.NUM_STEPS])

    for ilat_ind, inlat in enumerate(inlats):
        for olon_ind, outlon in enumerate(outlons):
            curdir = os.path.join(input_dir,'in_%d/phi_mode2'%inlat)
            print "loading ", inlat, outlon

            phiN, phiS, latvec = load_phi_files(curdir, outlon, sc, zipped=True)
            
            # Get total particle fluxes (integrate over energy bins)
            N[ilat_ind, :, olon_ind, :] = np.sum(particle_flux(phiN, sc),axis=0).T
            S[ilat_ind, :, olon_ind, :] = np.sum(particle_flux(phiS, sc),axis=0).T

    obj['N'] = N
    obj['S'] = S

    print "Saving database"
    with open(output_filename,'wb') as f:
        pickle.dump(obj,f,pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    build_database(input_dir="/shared/users/asousa/WIPP/WIPPv4/outputs/agu2016_kp0_v3/pwr_-10000",
                    output_filename="db_agu2016_mode2.pkl")
