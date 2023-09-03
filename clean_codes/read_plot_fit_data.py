import os
os.environ["CDF_LIB"] = "~/CDF/lib"
from spacepy import pycdf
from matplotlib import pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
from read_rpw_data import read_tnr_autoch_full
import matplotlib.colors as colors
import pandas as pd
import glob
import math
from mpfit import *
from utils import *
from utils_wind_marsis import get_wind_data 




#defining some matplot stuff for the standard font style
rfont = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rfont)
plt.rcParams["legend.labelspacing"] = 0.001
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

myFmt = mdates.DateFormatter('%H:%M')   #this is needed to define the format of the time axis in spectrograms
                                        #this way you only print the hour:minute 





if __name__ == "__main__":
        instr = 'STEREO' #choose the instrument: WIND / RPW / STEREO / SOLO

        
        path = '/Volumes/L7aler_HD/Final_master_project/real_observations/stereo_files/' 
        #path = '/Volumes/L7aler_HD/Final_master_project/real_observations/psp_files/L3/cdf/'
        #path = '/Volumes/L7aler_HD/Final_master_project/real_observations/wind_files/'
        #files = glob.glob(path)

        file = path + 'wi_wa_rad1_l3_df_20200711_v01.cdf'
        #file = path + 'psp_vk_fld_lfr_20200711_v01.cdf'
        file = path + 'sta_l3_wav_hfr_20200711_v01.cdf'

        burst_list = pd.read_csv('intervals/WIND_burst_intervals.csv', index_col = 0)  #you can change the file depening which instrument you use

        ti, tf = burst_list.B1i.iloc[0], burst_list.B1f.iloc[0] # approximate start and end time of the burst
        ti, tf = convert_string_datetime(ti), convert_string_datetime(tf)
        datum = burst_list.date.iloc[0]
        
    
        #the following opens the cdf files for the different spacecraft - gives Intensity, Time, and Frequency
        if instr == 'RPW':
                cdf = read_tnr_autoch_full(file)
                I = cdf['voltage'].T   
                T, F = cdf['time'][:], cdf['frequency'][:] * 1E3

                fchannels  = np.array([np.argmin((abs(F - 411E3))), np.argmin((abs(F - 430E3))),np.argmin((abs(F - 511E3))) \
                      , np.argmin((abs(F - 635E3))), np.argmin((abs(F - 663E3))), np.argmin((abs(F - 755E3))) \
                      , np.argmin((abs(F - 788E3))), np.argmin((abs(F - 979E3)))])
                fchannels = np.unique(fchannels)
        elif instr == 'WIND':
                cdf = pycdf.CDF(file)
                F, I, T = get_wind_data(cdf)
                F, I, T = F * 1E3, I.T, T.T   #freq units in Hz (kHz actually)

                #need to do this again because wind data is inconsistent with the ordering of the frequency... nervt mech!
                fchannels  = np.array([np.argmin((abs(F[:, 0] - 411E3))), np.argmin((abs(F[:, 0] - 430E3))),np.argmin((abs(F[:, 0] - 511E3))) \
                      , np.argmin((abs(F[:, 0] - 635E3))), np.argmin((abs(F[:, 0] - 663E3))), np.argmin((abs(F[:, 0] - 755E3))) \
                      , np.argmin((abs(F[:, 0] - 788E3))), np.argmin((abs(F[:, 0] - 979E3)))])
                fchannels = np.unique(fchannels)
                F = F.T
        else:
                cdf = pycdf.CDF(file)
                T, I, F = cdf['Epoch'][:], cdf['PSD_FLUX'][:] * 1E22, cdf['FREQUENCY'][:] #intensity units are in sfu, freq units are in Hz
                fchannels  = np.array([np.argmin((abs(F - 411E3))), np.argmin((abs(F - 430E3))),np.argmin((abs(F - 511E3))) \
                      , np.argmin((abs(F - 635E3))), np.argmin((abs(F - 663E3))), np.argmin((abs(F - 755E3))) \
                      , np.argmin((abs(F - 788E3))), np.argmin((abs(F - 979E3)))])
                fchannels = np.unique(fchannels)


        


        #the follwoing visulaizes the data, either the time profile (with fitting), or the dynamic soectrum
        vis = input('What do you want to see? Spectograph (s) / Time Profile (p)')

        if vis == 'p':
                df_name = 'test.csv'
                df = pd.DataFrame(columns = ['date', 'ti', 'tf', 'decay_time', 'error_low', 'error_up', 'continuum_level', 'frequency(Hz)', 'fpeak_sfu', 'burst_number', 'fit_method'])
                df.to_csv(df_name)

                _ = fchannel_loop_decay_time(fchannels, ti, tf, I, T, F, df_name, datum, 1, instr)
        else:
                if instr == 'WIND':
                        cdf = pycdf.CDF(file)
                        T, I, F = cdf['Epoch'][:], cdf['STOKES_I'][:] *1E22, cdf['FREQUENCY'][:] * 1E3
                        F = F[:64]  #there are 64 frequency channels (though some channels like 1040 and 940 are there multiple times for some reason)
                        I_new = np.zeros((len(I) // 64, 64))
                        F = F[:64]  #there are 64 frequency channels (though some channels like 1040 and 940 are there multiple times for some reason)
                        I_new = np.zeros((len(I) // 64, 64))
                        try:
                                for i in range(64):
                                        I_new[:, i] = I[i::64]
                                I = I_new
                                T = T[::64]
                                plotting_dynamic_spectrum(ti, tf, datum, I, T, F, instr)
                        except:
                                print('Wind structure data is problematic')
                else:
                        plotting_dynamic_spectrum(ti, tf, datum, I, T, F, instr)

                



      

       
















