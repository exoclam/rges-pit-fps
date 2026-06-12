import numpy as np
import pandas as pd
import os
import astropy
from astropy.io import fits
from glob import glob
#import VBMicrolensing
import RTModel
import math

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 1})
pylab_params = {'legend.fontsize': 'large',
    'axes.labelsize': 'x-large',
    'axes.titlesize': 'x-large',
    'xtick.labelsize': 'large',
    'ytick.labelsize': 'large'}
pylab.rcParams.update(pylab_params)

path = "/Users/lamchris/Desktop/wg6/wg6/"

times = path + 'data/times/'
f146_times = np.load(times + 'roman_times_shortcadence.npy')
f087_times = np.load(times + 'roman_times_longcadence.npy')
f213_times = np.load(times + 'roman_times_longcadence2.npy')
f146_times = f146_times - 2450000
f087_times = f087_times - 2450000
f213_times = f213_times - 2450000

t2cep_path = path + 'data/RGES_filters_T2CEP_lightcurves/'
t2cep_path_processed = path + 'data/processed/RGES_filters_T2CEP_lightcurves/'
t2cep_files = sorted(glob(t2cep_path+'*'))

fl_path = path + 'data/RGES_filters_FL_lightcurves/RGES_filters_FL_lightcurves/'
fl_path_processed = path + 'data/processed/RGES_filters_FL_lightcurves/RGES_filters_FL_lightcurves/'
fl_files = sorted(glob(fl_path+'*'))

for i in range(len(fl_files)):
    fits_image_filename = fl_files[i]
    # strip preamble and suffix from file name
    stripped_name = fits_image_filename.removeprefix(fl_path+'RGES_filters_')
    stripped_name = stripped_name.removesuffix('_lightcurves.fits')

    event_path = fl_path_processed+stripped_name+'/'
    event_data_path = fl_path_processed+stripped_name+'/Data/'
    try:
        os.mkdir(event_path)
        print(f"Directory '{event_path}' created")
    except:
        print(f"Directory '{event_path}' already exists")
    try:
        os.mkdir(event_data_path)
        print(f"Directory '{event_data_path}' created")
    except:
        print(f"Directory '{event_data_path}' already exists")

    # open fits data
    hdul = fits.open(fits_image_filename)
    # print(hdul[0].header)
    # print(hdul[1].header)
    # print(hdul[2].header)
    # print(hdul[3].header)
    # print(hdul[1].data)
    # print(hdul[2].data)
    # print(hdul[3].data)
    f087_y, f087_yerr = zip(*hdul[1].data)
    f146_y, f146_yerr = zip(*hdul[2].data)
    f213_y, f213_yerr = zip(*hdul[3].data)

    df_f087 = pd.DataFrame({'# Mag': f087_y, 'err': f087_yerr, 'HJD-2450000': f087_times})
    df_f146 = pd.DataFrame({'# Mag': f146_y, 'err': f146_yerr, 'HJD-2450000': f146_times})
    df_f213 = pd.DataFrame({'# Mag': f213_y, 'err': f213_yerr, 'HJD-2450000': f213_times})
    df_f087.to_csv(event_data_path+'w087sat3.dat', sep=' ', index=False)
    df_f146.to_csv(event_data_path+'w146sat1.dat', sep=' ', index=False)
    df_f213.to_csv(event_data_path+'w213sat2.dat', sep=' ', index=False)

    rtm = RTModel.RTModel(event_path)
    rtm.set_satellite_dir('/satellitedir')
    rtm.config_InitCond(npeaks=5, nostatic = False, usesatellite = 1)
    rtm.run()
    quit()

    plt.errorbar(f146_times, f146_y, f146_yerr, label='F146')
    plt.errorbar(f087_times, f087_y, f087_yerr, label='F087')
    plt.errorbar(f213_times, f213_y, f213_yerr, label='F213')
    plt.xlabel('time [BJD-2450000]')
    plt.ylabel('magnitude')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    quit()

    ### fit binary lens event
    #VBM = VBMicrolensing.VBMicrolensing()
    rtm = RTModel.RTModel()
    s = 0.9       # Separation between the lenses
    q = 0.1       # Mass ratio
    u0 = 0.0       # Impact parameter with respect to center of mass
    alpha = 1.0       # Angle of the source trajectory
    rho = 0.01       # Source radius
    tE = 1.0      # Einstein time in days
    t0 = 1663.45      # Time of closest approach to center of mass

    # Array of parameters. Note that s, q, rho and tE are in log-scale
    pr = [math.log(s), math.log(q), u0, alpha, math.log(rho), math.log(tE), t0]

    #t = np.linspace(t0-tE, t0+tE, 300) # Array of times
    t = f146_times

    magnifications, y1, y2 = VBM.BinaryLightCurve(pr,t)      # Calculation of binary-lens light curve
    plt.plot(t, magnifications)
    plt.show()
    quit()

for i in range(len(t2cep_files)):
    fits_image_filename = t2cep_files[i] 
    hdul = fits.open(fits_image_filename)
    # print(hdul[0].header)
    # print(hdul[1].header)
    # print(hdul[2].header)
    # print(hdul[3].header)
    # print(hdul[1].data)
    # print(hdul[2].data)
    # print(hdul[3].data)
    
    f087_y, f087_yerr = zip(*hdul[1].data)
    f146_y, f146_yerr = zip(*hdul[2].data)
    f213_y, f213_yerr = zip(*hdul[3].data)

    plt.errorbar(f146_times, f146_y, f146_yerr, label='F146')
    plt.errorbar(f087_times, f087_y, f087_yerr, label='F087')
    plt.errorbar(f213_times, f213_y, f213_yerr, label='F213')
    plt.xlabel('time [BJD-]')
    plt.ylabel('magnitude')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()


    quit()