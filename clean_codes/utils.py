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
import sys

from mpfit import *


#defining some matplot stuff for the standard font style

rfont = {"font.family" : "serif", 
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rfont)
plt.rcParams["legend.labelspacing"] = 0.001
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

myFmt = mdates.DateFormatter('%H:%M')   #this is needed to define the format of the time axis in spectrograms
                                        #this way you only print the hour:minute 

xlim_coord, fpeak_xcoord = [], 0
ylim_coord = []


def on_press(event):
    global xlim_coord, ylim_coord, fpeak_xcoord
    if event.key == 'x':
        print('Data recorded:', event.xdata)
        xlim_coord.append(event.xdata)
        ylim_coord.append(event.ydata)
    if event.key == 'd':
        print('Delete last record')
        xlim_coord = xlim_coord[:-1]
        ylim_coord = ylim_coord[:-1]
    if event.key == 'a':
        print('New_peak')
        fpeak_xcoord = event.xdata
        




def get_hour_minute(time):
    #function which returns the hour and minute of the burst in integer format (rather than string)
    #time: string, gives the approximate time (hour:min) when the burst occured 
    t = time.split(':')
    return int(t[0]), int(t[1])

def convert_string_datetime(t):
    dt_tuple=tuple([int(x) for x in t[:10].split('-')])+tuple([int(x) for x in t[11:19].split(':')])   #had 11:18 before
    datetimeobj=datetime.datetime(*dt_tuple)
    #print(type(datetimeobj))
    return datetimeobj
    


def fchannel_loop(fchannels, ti, tf, I, T, F, df_name, datum):
    #function which loops through the frequency channels and then updates the dataframe which returns the peak flux and corresponding time of the burst
    #fchannels: list of indexes of frequency channels to analyse
    #ti and tf: datetime variables, marking the range where the burst approximately takes place
    #I: array all the flux intensity data (in sfu)
    #T: array of the epochs/time (datetime format)
    #F: array of all the frequency data (in Hz)
    #df_name: name of the data frame where we want to save the data
    #datum: string, the date when the burst happended (format day/month/year)
    
    for i in fchannels:
        df = pd.read_csv(df_name, index_col = 0)
        
        flux = I[:, i]   #flux in frequency channel I
        flux_max, i_max = get_peak_flux(ti, tf, T, flux)
        df_new = pd.DataFrame(columns = ['date', 'ti', 'tf', 'fpeak', 'tmax', 'frequency(Hz)'])
        df_new.date, df_new.ti, df_new.tf, df_new.fpeak, df_new.tmax, df_new['frequency(Hz)'] = [datum], [ti], [tf], [flux_max], [T[i_max][0]], [F[i]]
        df = pd.concat([df, df_new])
        df.to_csv(df_name)
        
    return _


def get_peak_flux(ti, tf, date, flux):
    #function which returns the peak flux and corresponding time index of a burst
    #ti and tf: datetime variables, marking the range where the burst approximately takes place
    #date: array of the epochs/time (datetime format)
    #flux: flux array withing a specific frequency channel for which we want to find the peak
    i = np.where((date > ti) & (date < tf))[0] 

    f_peak = np.amax(flux[i])
    idx = np.where(flux == f_peak)
    return f_peak, idx





def integer_to_datetime(t_int, t_datetime):
    year, month, day = t_datetime.year, t_datetime.month, t_datetime.day
    hour = t_int // 3600
    minute = (t_int - hour * 3600) // 60
    second = math.floor(t_int - hour * 3600 - minute * 60)
    microsecond = (t_int - hour * 3600 - minute * 60 - seconds) * 1E-6



def to_integer(time):
    #gives the time over which the bursts occurs in integer form

    ti, tf = time[0], time[-1]
    dt = (tf - ti).seconds
    #t = np.linspace(0, dt, len(time))

    t = np.array([t.hour * 3600 + t.minute * 60 + t.second for t in time])
    return t

def comp_bkgr_err(I):
    #measuring the background continuum level and the associated error
    m = np.median(I) #m = np.median(I[i])
    i = np.where(I - 2*m < 0)   #experiment with this a litttle
    err = np.std(I[i])
    return m, err

def get_time_idx(ti, tf, date):
    #function which returns time indices around which the burst occurs
    #ti and tf: datetime variables, marking the range where the burst approximately takes place
    #date: array of the epochs/time (datetime format)
    
    i = np.where((date > ti) & (date < tf))[0] 
    return i



def background_substr(I, i):
    if i[0] < 30:
        return np.median(I[:i[0]])
    else:
        return np.median(I[i[0] - 30:i[0]])




def fchannel_loop_decay_time(fchannels, ti, tf, I, T, F, df_name, datum, burst_num, instr = 'None'):
    #ti and tf: datetime variables, marking the range where the burst approximately takes place
    global xlim_coord, ylim_coord, fpeak_xcoord
    print('============================================')
    print('Burst at date: ', datum)
    print('Around time:', tf)

    F_copy, I_copy, T_copy = np.copy(F), np.copy(I), np.copy(T)

    for j in fchannels:
        accord = False
        while not accord:
            #try:
                xlim_coord, ylim_coord = [], []
                
                if instr == 'WIND':
                    nz = np.where(F_copy[:, j] != 0)[0]
                    F, flux, T = F_copy[nz, j], I_copy[nz, j], T_copy[nz, j]
                else:
                    flux = I[:, j]   #flux in frequency channel I

                

                df = pd.read_csv(df_name, index_col = 0)

                i = get_time_idx(ti, tf, T)
                

                t = to_integer(T[i])  #time in integer format (in seconds only)
                if instr == 'WIND':

                    Iplot = flux[i] - background_substr(flux, i)
                    flux_max, _ = get_peak_flux(ti, tf, T, flux - background_substr(flux, i) )
                else:
                    Iplot = I[i, j] - background_substr(I[:, j], i) 
                    flux_max, _ = get_peak_flux(ti, tf, T, I[:, j] - background_substr(I[:, j], i))    #background subtraction accounted for

                
                
                i_fp = np.where(Iplot == flux_max) #index of the 'peak'
                fpeak_xcoord = t[i_fp][0]
                


                m, error = comp_bkgr_err(I[i, j])

                fig, ax = plt.subplots(1, 1, figsize = (8, 6))

                f = fig.canvas.mpl_connect('key_press_event', on_press)
                
                ax.plot(t, Iplot, color = 'black', lw = 1.2)
                ax.scatter(t, Iplot, color = 'black', lw = 1.2)
                ax.scatter(t[i_fp], Iplot[i_fp], color = 'magenta')
                ax.set_yscale('log')

                plt.show()

                i = [np.argmin(abs(t - xlim_coord[0])), np.argmin(abs(t - xlim_coord[1]))]
                
                #caused problems with WIND data for some reason
                #d = np.sqrt(t**2 + Iplot**2)
                #i = [np.argmin(abs(d - np.sqrt(xlim_coord[0]**2 + ylim_coord[0]**2))), np.argmin(abs(d - np.sqrt(xlim_coord[1]**2 + ylim_coord[1]**2)))]
                i_fp = np.argmin(abs(t - fpeak_xcoord))

                
                

                fig, ax = plt.subplots(1, 1, figsize = (8, 6))
                
                ax.plot(t, Iplot, color = 'black', lw = 1.2)
                ax.scatter(t, Iplot, color = 'black', lw = 1.2)
                ax.scatter(t[i[0]], Iplot[i[0]], color = 'red')
                ax.scatter(t[i[-1]], Iplot[i[-1]], color = 'red')
                ax.scatter(t[i_fp], Iplot[i_fp], color = 'magenta')

                ax.set_yscale('log')
                plt.show()


                #fitting
                Ifit, tfit = Iplot[i[0]:i[-1] + 1], t[i[0]:i[-1] + 1]
                #Ipeak, imax = np.amax(Ifit), np.argmax(Ifit)
                Ipeak, imax = Ifit[0], 0
                tpeak = tfit[imax]


                err = np.zeros(len(tfit)) + error
                values = [100] #the initial guess for the decay time

                fit_method = input('What do you want to minimize? Absolute error (enter a)? Relative error (enter r)? Compare the two methods (c)?')
                if fit_method == 'r':
                    tau_decay = fitting(tfit, Ifit, err, values, efitter_rel)
                elif (fit_method != 'r') or (fit_method != 'c'):
                    tau_decay = fitting(tfit, Ifit, err, values, efitter_abs)
                if fit_method == 'c':
                    tau_decay_rel = fitting(tfit, Ifit, err, values, efitter_rel)
                    tau_decay_abs = fitting(tfit, Ifit, err, values, efitter_abs)
                    cont = np.amin(Ifit)
                    quick_plot_fit(Iplot, t, cont, tau_decay_rel, tau_decay_abs, Ipeak, tpeak, tfit, datum, F[j])

                if fit_method != 'c':
                    #calculate the uncertainties on the fit
                    err = tau_decay.perror[0]
                    if fit_method == 'r':
                        err_up, err_low = compute_error_v2(Iplot, t, i, tau_decay.params[0], efitter_rel)
                    else:
                        err_up, err_low = compute_error_v2(Iplot, t, i, tau_decay.params[0], efitter_abs)
                    #print('Fitting error:', err)
                    tot_err_up = np.sqrt(err_up**2 + err**2)
                    tot_err_low = np.sqrt(err_low**2 + err**2)
                    print('Fitted decay time:', tau_decay.params[0])
                    print('Uncertainties are +-', tot_err_up, tot_err_low)


                    cont = np.amin(Ifit)
                    plotting_fit(Iplot, t, cont, tau_decay, Ipeak, tpeak, tfit, datum, F[j], burst_num, save = False)

                    decision = input('Are you happy with the fit (if you enter ''y'', it means you agree to continue)')
                    if decision == 'y':
                        accord = True
                        plotting_fit(Iplot, t, cont, tau_decay, Ipeak, tpeak, tfit, datum, F[j], burst_num, save = True)

                        t_all = to_integer(T)   #the entire day in integer format (time in integer format, i.e in seconds)
                        i, q = np.argmin(abs(t_all - tfit[0])), np.argmin(abs(t_all - tfit[-1]))
                        tfit_i, tfit_f = T[i], T[q]   #interval where we performed the fit

                        df_new = pd.DataFrame(columns = ['date', 'ti', 'tf', 'decay_time', 'error_low', 'error_up', 'continuum_level', 'frequency(Hz)', 'fpeak_sfu', 'burst_number', 'fit_method'])
                        df_new.date, df_new.ti, df_new.tf, df_new.decay_time, df_new.error_low, df_new.error_up, df_new.continuum_level, df_new['frequency(Hz)'], df_new.fpeak_sfu, \
                        df_new.burst_number, df_new.fit_method = [datum], [tfit_i], [tfit_f], [tau_decay.params[0]], tot_err_low, tot_err_up, [cont], [F[j]], Iplot[i_fp], burst_num, fit_method
                        df = pd.concat([df, df_new])
                        df.to_csv(df_name)
            #except:
            #    if len(xlim_coord) < 2:
            #        print('Error: You have to choose a time interval for the fit!')
            #        stop = input('Do you want to stop the code (answer with y if you want to stop?; enter s if you want to skip to the next frequency channel)')
            #        if stop == 'y':
            #            sys.exit()
            #        elif stop == 's':
            #            accord = True
            
        
    return 'we are done here'








def plotting_fit(I, t, cont, tau_decay, Ipeak, tpeak, tfit, datum, F, burst_num, save = False):
    fig, ax = plt.subplots(1, 1, figsize = (8, 6))
    
    ax.scatter(t, I, color = 'black', lw = 1.2)
    ax.plot(t, I, color = 'black', lw = 1.2, label = 'Signal')
    ax.plot(t, np.zeros(len(t)) + cont, color = 'grey', lw = 0.9, label = 'Background')

    y = eplot(tau_decay.params, tfit, Ipeak, tpeak)
    err = tau_decay.perror[0]
    tau = str(np.round(tau_decay.params[0], 2)) + 's'
    ax.plot(tfit, y, ls = '--', color = 'orange', lw = 1.2, label = r'$\tau$ = ' + tau)
    ax.set_title('Date: ' + datum + ', frequency ' + str(F) + 'Hz')

    ax.legend(frameon = False, fontsize = 12)
    ax.minorticks_on()
    ax.tick_params(which = 'both', bottom = True, top = True, left = True, right = True)
    ax.tick_params(which = 'major', length = 10, direction = 'in', labelsize = 10)
    ax.tick_params(which = 'minor', length = 5, direction = 'in', labelsize = 10)
    ax.set_yscale('log')
    ax.set_xlabel('Time (s)', fontsize = 15)
    ax.set_ylabel(r'Flux Density (sfu)', fontsize = 15)

    if not save:
        plt.show()
    if save:
        path = '/Users/louissiebenaler/Dropbox/Private_7aler/Louis/University Material/Leiden University/Final_Research_Project/real_observations/analysis/'
        folder = path + 'WIND_Decay_Time_Fits/'
        date = datum.replace("/", "_" )
        name = date + '__' + str(F) + 'Hz' + '__' + str(burst_num) + '.png'
        plt.savefig(folder + name)
        plt.close(fig)


def quick_plot_fit(I, t, cont, tau_decay_rel, tau_decay_abs, Ipeak, tpeak, tfit, datum, F):
    fig, ax = plt.subplots(1, 1, figsize = (8, 6))
    
    ax.scatter(t, I, color = 'black', lw = 1.2)
    ax.plot(t, I, color = 'black', lw = 1.2, label = 'Signal')
    ax.plot(t, np.zeros(len(t)) + cont, color = 'grey', lw = 0.9, label = 'Background')

    y = eplot(tau_decay_rel.params, tfit, Ipeak, tpeak)
    tau = str(np.round(tau_decay_rel.params[0], 2)) + 's'
    ax.plot(tfit, y, ls = '--', color = 'orange', lw = 1.2, label = r'Rel. error, $\tau$ = ' + tau)

    y = eplot(tau_decay_abs.params, tfit, Ipeak, tpeak)
    err = tau_decay_abs.perror[0]
    tau = str(np.round(tau_decay_abs.params[0], 2)) + 's'
    ax.plot(tfit, y, ls = '--', color = 'green', lw = 1.2, label = r'Abs. error, $\tau$ = ' + tau)

    ax.set_title('Date: ' + datum + ', frequency ' + str(F) + 'Hz')

    ax.legend(frameon = False, fontsize = 12)
    ax.minorticks_on()
    ax.tick_params(which = 'both', bottom = True, top = True, left = True, right = True)
    ax.tick_params(which = 'major', length = 10, direction = 'in', labelsize = 10)
    ax.tick_params(which = 'minor', length = 5, direction = 'in', labelsize = 10)
    ax.set_yscale('log')
    ax.set_xlabel('Time (s)', fontsize = 15)
    ax.set_ylabel(r'Flux Density (sfu)', fontsize = 15)

    plt.show()
    
    
    

#======================================================================
#here are the functions important for the fitting, which require the mpfit package


def fitting(tfit, Ifit, err, values, efit):
    fa = {'t':tfit, 'I':Ifit, 'tpeak':tfit[0], 'Ipeak':Ifit[0], 'err':err}
    parinfo = [{'value':1., 'fixed':0, 'limited':[0, 0], 'limits':[0.,0.], 'tied':''} for i in range(1)]

    parinfo[0]['value'] = values[0]
    
    return mpfit(efit, functkw=fa, parinfo = parinfo, maxiter = 400, quiet = True, ftol = 1e-15)

def efitter_abs(p, fjac = None, t = None, I = None, Ipeak = None, tpeak = None, err = None):
    model = Ipeak * np.exp(- (t - tpeak) / p[0])
    status = 0
    return [status, (I - model)/ err]

def efitter_rel(p, fjac = None, t = None, I = None, Ipeak = None, tpeak = None, err = None):
    model = Ipeak * np.exp(- (t - tpeak) / p[0])
    status = 0
    return [status, abs(I - model)/ (I*err)]


def eplot(p, t, Ipeak, tpeak):
    return Ipeak * np.exp(- (t - tpeak) / p[0])



#======================================================================
#functions to estimate the error on the fitted decay time

def fitting_for_error(Iplot, t, i, jl, ju, efitter):
    Ifit, tfit = Iplot[i[0] + jl:i[-1] + ju], t[i[0] + jl:i[-1] + ju]
    if len(tfit) > 2:
        #Ipeak, imax = np.amax(Ifit), np.argmax(Ifit)
        Ipeak, imax = Ifit[0], 0
        tpeak = tfit[imax]
        err = np.zeros(len(tfit)) + 1
        values = [50] #the initial guess for the decay time
        tau_decay = fitting(tfit, Ifit, err, values, efitter)
        return tau_decay.params[0]
    else:
        return -1


def compute_error_v2(Iplot, t, i, tau, efitter):
    tau_decay_list = []

    for j in range(1, 5):
        x = fitting_for_error(Iplot, t, i, 0, j, efitter)
        if (abs(tau - x) / tau < 0.2) and (x > 0):
            tau_decay_list.append(x)

        x = fitting_for_error(Iplot, t, i, 0, -j, efitter)
        if abs(tau - x) / tau < 0.2 and (x > 0): 
            tau_decay_list.append(x)

        x = fitting_for_error(Iplot, t, i, j, 1, efitter)
        if (abs(tau - x) / tau < 0.2) and (x > 0):
            tau_decay_list.append(x)

        x = fitting_for_error(Iplot, t, i, -j, 1, efitter)
        if (abs(tau - x) / tau < 0.2) and (x > 0):
            tau_decay_list.append(x)

        x = fitting_for_error(Iplot, t, i, j, j, efitter)
        if (abs(tau - x) / tau < 0.2) and (x > 0):
            tau_decay_list.append(x)

        x = fitting_for_error(Iplot, t, i, -j, -j, efitter)
        if (abs(tau - x) / tau < 0.2) and (x > 0):
            tau_decay_list.append(x)

        x = fitting_for_error(Iplot, t, i, j, -j, efitter)
        if (abs(tau - x) / tau < 0.2) and (x > 0):
            tau_decay_list.append(x)

        x = fitting_for_error(Iplot, t, i, -j, j, efitter)
        if (abs(tau - x) / tau < 0.2) and (x > 0):
            tau_decay_list.append(x)


    tau_list = np.array(tau_decay_list)
    #print('Tau_list:', tau_list)
    err_up, err_low = abs(tau - np.amax(tau_list)), abs(tau - np.amin(tau_list))
    return err_up, err_low



#function for plotting the dynamic spectrum of the burts

def plotting_dynamic_spectrum(ti, tf, datum, I, T, F, instr):
    ti-= datetime.timedelta(minutes=30)
    tf+= datetime.timedelta(minutes=30)
    
    fig, ax = plt.subplots(1, 1, figsize = (12, 5))
    i = get_time_idx(ti, tf, T)
    y, x = np.meshgrid(F, T[i])  

    I_min, I_max = -np.abs(I[i]).max(), np.abs(I[i]).max()   #not used

    c = ax.pcolormesh(x, y, I[i], cmap='rainbow', norm=colors.LogNorm(vmin=(1E2), vmax=1E7))

    ax.xaxis.set_major_formatter(myFmt)
    ax.set_xlabel('Time (hour:minute)', fontsize = 12)
    ax.set_ylabel('Frequency (Hz)', fontsize = 12)
    ax.axis([x.min(), x.max(), y.min(), 1])
    cb = plt.colorbar(c)  #producing the colorbar s
    cb.set_label(label= '{} (sfu)'.format(instr), size = 12)
    #plt.savefig('Spectrograms/Sophie_Paper/Dynamic_spectrum_SolO_processed.pdf', bbox_inches='tight')
    ax.set_ylim(100E3, 970E3)
    



    plt.show()






    