# -*- coding: utf-8 -*-
## +
# NAME:
#   read_TNR_data
#
# PURPOSE:
#   READ L2 data from TNR
#
#
# INPUTS:
#   filepath - STRING OF THE FILENAME IN .CDF CONTAINING THE L2 DATA 
#
#   sensor : TNR SENSOR TO BE READ
#             	1: V1
#             	2:V2
#				3:V3
#				4:V1-V2
#				5:V2-V3
#				6:V3-V1
#				7: B
#
#   indstart: starting time index
#
#   indend:  end time index. indend=-99 means the last point of the timeserie.
#
# OUTPUTS:
#    V: array (time,frequency) of the measured signals
#    time: time of each measurement in Julian days
#    freq_tnr: frequency in kHz
#
# OPTIONAL OUTPUTS:
#   None.
#
#NOTES:
#  The script check if there are data from the two channel and put them together.
#
# MODIFICATION HISTORY:
#   Written by A.VECCHIO (LESIA, CNRS - RRL Radboud University): February 2021
#
def read_tnr_autoch_full(
    filepath, sensor=4, start_index=0, end_index=-99, data_index=0):
    from spacepy import pycdf
    import numpy as np
    import datetime

    with pycdf.CDF(filepath) as data_L2:

        freq_tnr1 = np.append(
            data_L2['TNR_BAND_FREQ'][0, :], data_L2['TNR_BAND_FREQ'][1, :]
        )
        freq_tnr2 = np.append(
            data_L2['TNR_BAND_FREQ'][2, :], data_L2['TNR_BAND_FREQ'][3, :]
        )
        freq_tnr = np.append(freq_tnr1, freq_tnr2)
        freq_tnr = freq_tnr / 1000.0  # frequency in kHz
        nn = np.size(data_L2['Epoch'][:])
        if end_index == -99:
            end_index = nn
        epochdata = data_L2['Epoch'][start_index:end_index]
        sensor_config = np.transpose(
            data_L2['SENSOR_CONFIG'][start_index:end_index, :]
        )
        auto1_data = np.transpose(data_L2['AUTO1'][start_index:end_index, :])
        auto2_data = np.transpose(data_L2['AUTO2'][start_index:end_index, :])
        sweep_num = data_L2['SWEEP_NUM'][start_index:end_index]
        bande = data_L2['TNR_BAND'][start_index:end_index]
        if sensor == 7:
            auto1_data = np.transpose(
                #data_L2['MAGNETIC_SPECTRAL_POWER1'][start_index:end_index, :]   #this is to get voltage
                data_L2['FLUX_DENSITY1'][start_index:end_index, :]               #this is to get flux
            )
            auto2_data = np.transpose(
                #data_L2['MAGNETIC_SPECTRAL_POWER2'][start_index:end_index, :]   #this is to get voltage
                data_L2['FLUX_DENSITY2'][start_index:end_index, :]               #this is to get flux
            )
        puntical = (data_L2['FRONT_END'][start_index:end_index] == 1).nonzero()
    epochdata = epochdata[puntical[0]]
    sensor_config = sensor_config[:, puntical[0]]
    auto1_data = auto1_data[:, puntical[0]]
    auto2_data = auto2_data[:, puntical[0]]
    sweep_numo = sweep_num[puntical[0]]
    bande = bande[puntical[0]]
    sweep_num = sweep_numo
    timet=epochdata
    #deltasw = sweep_numo[ 1:: ] - sweep_numo[ 0:np.size ( sweep_numo ) - 1 ]
    deltasw = abs (np.double(sweep_numo[ 1::]) - np.double(sweep_numo[ 0:np.size(sweep_numo)-1 ]))
    xdeltasw = np.where ( deltasw > 100 )
    xdsw = np.size ( xdeltasw )
    if xdsw > 0:
        xdeltasw = np.append ( xdeltasw, np.size ( sweep_numo ) - 1 )
        nxdeltasw = np.size ( xdeltasw )
        for inswn in range ( 0, nxdeltasw - 1 ):
            #sweep_num[ xdeltasw[ inswn ] + 1:xdeltasw[ inswn + 1 ] ] = sweep_num[
            #                                                           xdeltasw[ inswn ] + 1:xdeltasw[ inswn + 1 ] ] + \
            #                                                           sweep_num[ xdeltasw[ inswn ] ]
            sweep_num[ xdeltasw[ inswn ] + 1:xdeltasw[ inswn + 1 ] + 1 ] = sweep_num[
                                                                           xdeltasw[ inswn ] + 1:xdeltasw[
                                                                                                     inswn + 1 ] + 1 ] + \
                                                                           sweep_numo[ xdeltasw[ inswn ] ]
    sens0 = (sensor_config[0, :] == sensor).nonzero()[0]
    sens1 = (sensor_config[1, :] == sensor).nonzero()[0]
    psens0 = np.size(sens0)
    psens1 = np.size(sens1)

    if (np.size(sens0) > 0 and np.size(sens1) >0):
        auto_calib = np.hstack ((auto1_data[ :, sens0 ], auto2_data[:, sens1 ]))
        sens = np.append(sens0, sens1)
        timet_ici = np.append(timet[ sens0 ], timet[ sens1 ] )
    else:
        if (np.size(sens0) > 0):
            auto_calib = auto1_data[ :, sens0 ]
            sens = sens0
            timet_ici = timet[ sens0 ]
        if  (np.size(sens1) > 0):
            auto_calib = auto2_data[ :, sens1 ]
            sens = sens1
            timet_ici = timet[ sens1 ]
        if (np.size(sens0) == 0 and np.size(sens1) == 0):
            print('no data at all ?!?')
            V = (128, 128)
            V = np.zeros ( V ) + 1.0
            time = np.zeros ( 128 )
            sweepn_TNR = 0.0
            return {
                'voltage': V,
                'time': time,
                'frequency': freq_tnr,
                'sweep': sweepn_TNR,
                'sensor': sensor,
            }
    ord_time = np.argsort ( timet_ici )
    timerr = timet_ici[ ord_time ]
    sens = sens[ ord_time ]
    bandee = bande[ sens ]
    auto_calib=auto_calib[:,ord_time]
    maxsweep = max ( sweep_num[ sens ] )
    minsweep = min ( sweep_num[ sens ] )
    sweep_num = sweep_num[ sens ]
    V1 = np.zeros(128)
    V = np.zeros(128)
    time = 0.0
    sweepn_TNR = 0.0
    for ind_sweep in range ( minsweep, maxsweep + 1 ):
        ppunt = (sweep_num == ind_sweep).nonzero ()[ 0 ]
        xm = np.size ( ppunt )
        if xm > 0:
            for indband in range ( 0, xm ):
                V1[
                32
                * bandee[ ppunt[ indband ] ]: 32
                                              * bandee[ ppunt[ indband ] ]
                                              + 32
                ] = np.squeeze ( auto_calib[ :, [ ppunt[ indband ] ] ] )

        if np.sum ( V1 ) > 0.0:
            V = np.vstack ( (V, V1) )
            sweepn_TNR = np.append ( sweepn_TNR, sweep_num[  ppunt[ 0 ] ]  )
        V1 = np.zeros ( 128 )
        if xm > 0:
            time = np.append ( time, timerr[ min ( ppunt ) ] )
    V = np.transpose ( V[ 1::, : ] )
    time = time[ 1:: ]
    sweepn_TNR = sweepn_TNR[ 1:: ]
    return {
    'voltage': V,
    'time': time,
    'frequency': freq_tnr,
    'sweep': sweepn_TNR,
    'sensor': sensor,}