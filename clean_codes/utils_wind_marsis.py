import numpy as np
import pandas as pd
import re
from datetime import datetime



def convert_string_freq(freq_string):
    # Split the input string by whitespace to get a list of words
    words = freq_string[0].split()

    # Find the word that contains the floating point number
    number_word = [word for word in words if '.' in word][0]

    # Convert the number word to a float
    freq = float(number_word)

    return freq

def convert_string_flux(flux_string):
    # Split the input string by whitespace to get a list of string integers
    flux_floats = flux_string[0].split()

    # Use a list comprehension to convert each string integer into an actual integer
    flux = [float(i) for i in flux_floats]
    
    return flux

def convert_string_time(time_string):
    # Define a regex pattern to match the time in the string
    pattern = r'\d{4}-\d{3}T\d{2}:\d{2}:\d{2}\.\d{3}'

    # Use the regex `search()` function to find the first match of the pattern in the string
    match = re.search(pattern, time_string)

    # Extract the matched time string
    time_string = match.group()
    
    # Define a datetime format string that matches the input time string
    format_string = '%Y-%jT%H:%M:%S.%f'

    # Use the `strptime()` function from the datetime module to parse the time string
    time_obj = datetime.strptime(time_string, format_string)

    return time_obj


def get_marsis_data(df):
    i = 0
    for i_flux, i_freq, i_time in zip(np.arange(5, len(df), 6), np.arange(1, len(df), 6), np.arange(0, len(df), 6)):
        flux_string, freq_string = np.array(df.iloc[i_flux]), np.array(df.iloc[i_freq]) 
        flux, freq = convert_string_flux(flux_string), convert_string_freq(freq_string)

        time_string = np.array(df.iloc[i_time])
        time = convert_string_time(time_string[0])


        if i == 0:
            flux_array = flux
            freq_array = np.array([freq])

            time_array = np.array([time])
        elif i == 1:
            flux_array = np.concatenate(([flux_array], [flux]), axis=0)
            freq_array = np.concatenate((freq_array, [freq]))

            time_array = np.concatenate((time_array, [time]))
        else:
            flux_array = np.concatenate((flux_array, [flux]), axis=0)
            freq_array = np.concatenate((freq_array, [freq]))

            time_array = np.concatenate((time_array, [time]))
        i+= 1
        
        
        
    #the instrument has 160 frequency channels
    F = freq_array[:160]  #array of unique frequency channels
    T = time_array[::160] #array of observation time
    I = np.reshape(np.median(flux_array, axis = 1), (len(T), len(F)))  #intensity matrix
    #-> took the median of all the intensities for a given frequency-time instance
    return F, I, T



def get_wind_data(cdf):
    T = cdf['Epoch'][:]
    I = cdf['STOKES_I'][:] *1E22
    F = cdf['FREQUENCY'][:]

    #sort the arrays such that it goes from earlies time to latest
    i = np.argsort(T)
    T = T[i]
    I = I[i]
    F = F[i]



    #produce the Intensity, time, and frequency grids to produce lightcurves
    count = np.zeros(len(np.unique(F)))
    for i, f in enumerate(np.unique(F)):
        count[i] = len(np.where(F == f)[0])

    I_grid = np.zeros((len(np.unique(F)), int(np.amax(count))))
    F_grid = np.copy(I_grid)

    TMP = T[:int(np.amax(count))] - T[:int(np.amax(count))] 
    for i in range(len(np.unique(F) - 1)):
        tmp = T[:int(np.amax(count))] - T[:int(np.amax(count))] 
        TMP = np.concatenate((TMP, tmp), axis = 0)

    q = int(len(np.unique(F)) * np.amax(count))
    T_grid = np.reshape(TMP[:q], (len(np.unique(F)), int(np.amax(count))))

    for i, f in enumerate(np.unique(F)):
        j = np.where(F == f)[0]
        I_grid[i, :len(j)] = I[j]
        F_grid[i, :len(j)] = F[j]
        T_grid[i, :len(j)] = T[j]

    

    return F_grid, I_grid, T_grid

    
