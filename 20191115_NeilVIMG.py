# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 23:49:42 2019

@author: KEVIN-DORGANS
"""
import numpy as np
import pandas as pd
import scipy as sp
import igor.binarywave as bw
import easygui
import os
from matplotlib import pyplot as plt
from neo import io
from scipy import signal

def load_directory_content__(): #   This lists waves from selected directory
    directory = easygui.diropenbox(default=r'Z:\\')
    paths__ = []
    for i in range(len(os.listdir(directory))):
        paths__.append(directory+"\\"+os.listdir(directory)[i])
    return (paths__, directory)

def load_IgorWave__(path): #    This reads Igor wave
    r = io.IgorIO(filename=path)
    data = bw.load(path)
    timestamp = int(data["wave"]["wave_header"]["modDate"])
    analogsignal_number = int(data["wave"]['wave_header']['nDim'][1])
    sampling__ = float(r.read_analogsignal().sampling_rate)
    sweep_position__ = int(r.filename.split('\\R')[1].split('_')[0])
    if analogsignal_number > 1:
        sweep_lag = np.array(str(data["wave"]['note']).split('%SweepStartTimes^>')[1].split('|')[0:analogsignal_number], dtype=np.float64)
    else:
        sweep_lag = [0]
    return r, sampling__, sweep_position__, timestamp, sweep_lag


__IMAGING_DATA_LIST__ = []
__RAW_ePhy_LIST__ = []
SWEEP_TIME = []

ID__LIST = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5] #The cell_ID per trial list
Dye1_CELL_ID = [4, 5] #Cell IDs for condition1
Dye4_CELL_ID = [1, 2, 4] #Cell IDs for condition2

DOWNSAMPLING__ = 1
VOLTAGE_STEP_NUMBER = 15
ePhy_Imaging_DELAY__ = 2 #in secs.
Camera_SAMPLING_INTERVAL = 0.025497 #in sec.
SAMPLE_NUMBER_PER_TRIAL = 400
PROTOCOL_START_TIME_ePhy = 5.01 #in sec.
PROTOCOL_START_POINT_imaging = 120
STEP_LENGTH_time = 0.2 #in sec.
STEP_LENGTH_imaging = STEP_LENGTH_time/Camera_SAMPLING_INTERVAL #pulses in current step protocol
WINDOW_LENGTH_imaging = int(0.05/Camera_SAMPLING_INTERVAL) #calculation time window

#Parameters for final graph
NBINS__ = 8
MIN_RANGE = -150
MAX_RANGE = 150

#STEP1 - IMPORT IMAGING DATA from clipboard
temp_list = np.array(pd.read_clipboard())
temp_list = np.transpose(temp_list)
for i in range(len(temp_list)):
    __IMAGING_DATA_LIST__.append(temp_list[i])

#STEP1 - IMPORT ePhy DATA from folder (order has to match with IMAGING DATA)
PATHS__, DIRECTORY__ = load_directory_content__()
for i in range(len(PATHS__)):
    r, sampling__, sweep_position__, start_time, sweep_lag = load_IgorWave__(PATHS__[i])
    signal = r.read_analogsignal().transpose()
    print('LOADING -- '+PATHS__[i])
    for j in range(len(signal)):
        RAW = np.array(signal[j], dtype=np.float64)*1000
        __RAW_ePhy_LIST__.append(sp.signal.resample(RAW, int(len(RAW)/DOWNSAMPLING__)))
        SWEEP_TIME.append(start_time+sweep_lag[j])

#STEP2 - CALCULATE AVERAGE VALUES PER CELL FOR IMAGING and ePhy
fig = plt.figure(figsize=(3, 4))
ax = plt.subplot(211)
AVERAGE_VOLTAGE_STEP_VALUES = []
for i in range(len(__RAW_ePhy_LIST__)):
    ax.plot(np.linspace(-ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])/sampling__-ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])), __RAW_ePhy_LIST__[i], color='black', alpha=0.2)
    temp__voltage_step_average = []
    for j in np.arange(5, VOLTAGE_STEP_NUMBER):
        PROTOCOL_START_POINT_ePhy = int(PROTOCOL_START_TIME_ePhy* sampling__)+j*(STEP_LENGTH_time*sampling__)
        temp__voltage_step_average.append(np.nanmedian(__RAW_ePhy_LIST__[i][int(PROTOCOL_START_POINT_ePhy):int(PROTOCOL_START_POINT_ePhy+(0.1*sampling__))]))
    AVERAGE_VOLTAGE_STEP_VALUES.append(temp__voltage_step_average)
ax.set_xlabel('Time(sec.)')
ax.set_ylabel('Membrane voltage(mV)')

AVERAGE_fluo_STEP_VALUES = []
ax2 = plt.subplot(212, sharex=ax)
AVERAGE_FLUO_STEP_VALUES = []
for i in range(len(__IMAGING_DATA_LIST__)):
    BSL__ = np.nanmedian(__IMAGING_DATA_LIST__[i][100:120])
    SIG__ = sp.signal.medfilt(__IMAGING_DATA_LIST__[i], 1)/BSL__
    ax2.plot(np.linspace(0, len(SIG__)*Camera_SAMPLING_INTERVAL, len(SIG__)), SIG__*100, color='black', alpha=0.3, lw=1)


    temp__fluo_step_average = []
    for j in  np.arange(5, VOLTAGE_STEP_NUMBER):
        if np.nanmedian(SIG__[int(PROTOCOL_START_POINT_imaging +j*STEP_LENGTH_imaging):int(PROTOCOL_START_POINT_imaging+j*STEP_LENGTH_imaging+WINDOW_LENGTH_imaging)]) >= 1:
            temp__fluo_step_average.append(np.nanmax(SIG__[int(PROTOCOL_START_POINT_imaging +j*STEP_LENGTH_imaging):int(PROTOCOL_START_POINT_imaging+j*STEP_LENGTH_imaging+WINDOW_LENGTH_imaging)]))
        else:
            temp__fluo_step_average.append(np.nanmin(SIG__[int(PROTOCOL_START_POINT_imaging +j*STEP_LENGTH_imaging):int(PROTOCOL_START_POINT_imaging+j*STEP_LENGTH_imaging+WINDOW_LENGTH_imaging)]))
    AVERAGE_fluo_STEP_VALUES.append(temp__fluo_step_average)
ax2.set_xlabel('Time(sec.)')
ax2.set_ylabel('F/F')
plt.tight_layout()
#######################

#STEP3 - FINAL PLOT
plt.figure(figsize=(4, 4))
BIN_SIZE = (MAX_RANGE- MIN_RANGE)/ NBINS__
AVERAGE_DYE1 = []
AVERAGE_DYE4 = []
for m in range(len(Dye1_CELL_ID)+len(Dye4_CELL_ID)):
    AVERAGE_BINNED_TRACE = []
    BINNED_TRACE_TEMP__ = []
    for i in range(len(AVERAGE_VOLTAGE_STEP_VALUES)):
        if True:
            if ID__LIST[i] == m:
                counter = 0
                for k in np.linspace(MIN_RANGE, MAX_RANGE, NBINS__):
                    BIN_TEMP__ = []
                    for l in range(len(AVERAGE_VOLTAGE_STEP_VALUES[i])):
                        if k-BIN_SIZE/2 < AVERAGE_VOLTAGE_STEP_VALUES[i][l] < k+BIN_SIZE/2:
                            BIN_TEMP__.append(AVERAGE_fluo_STEP_VALUES[i][l])
                    try:
                        BINNED_TRACE_TEMP__[counter] = np.append(BINNED_TRACE_TEMP__[counter], BIN_TEMP__)
                    except:
                        BINNED_TRACE_TEMP__.append(BIN_TEMP__)
                    counter += 1
    for i in range(len(BINNED_TRACE_TEMP__)):
        AVERAGE_BINNED_TRACE.append(np.nanmean(BINNED_TRACE_TEMP__[i]))
    if m in Dye1_CELL_ID:
        plt.scatter(np.linspace(MIN_RANGE+(MAX_RANGE-MIN_RANGE)/NBINS__-BIN_SIZE/2, MAX_RANGE+(MAX_RANGE-MIN_RANGE)/NBINS__-BIN_SIZE/2, len(AVERAGE_BINNED_TRACE)), AVERAGE_BINNED_TRACE, color='orange', alpha=0.4, s=5)
        if len(AVERAGE_BINNED_TRACE) > 0:
            AVERAGE_DYE1.append(AVERAGE_BINNED_TRACE)
    else:
        plt.scatter(np.linspace(MIN_RANGE+(MAX_RANGE-MIN_RANGE)/NBINS__-BIN_SIZE/2, MAX_RANGE+(MAX_RANGE- MIN_RANGE)/ NBINS__-BIN_SIZE/2, len(AVERAGE_BINNED_TRACE)), AVERAGE_BINNED_TRACE, color='purple', alpha=0.4, s=5)
        if len(AVERAGE_BINNED_TRACE) > 0:
            AVERAGE_DYE4.append(AVERAGE_BINNED_TRACE)
plt.scatter(np.linspace(MIN_RANGE+(MAX_RANGE-MIN_RANGE)/NBINS__-BIN_SIZE/2, MAX_RANGE+(MAX_RANGE-MIN_RANGE)/NBINS__-BIN_SIZE/2, len(AVERAGE_DYE4[0])), np.nanmean(AVERAGE_DYE4, axis=0), facecolors='none', edgecolors='purple', alpha=1)
plt.scatter(np.linspace(MIN_RANGE+(MAX_RANGE-MIN_RANGE)/NBINS__-BIN_SIZE/2, MAX_RANGE+(MAX_RANGE-MIN_RANGE)/NBINS__-BIN_SIZE/2, len(AVERAGE_DYE1[0])), np.nanmean(AVERAGE_DYE1, axis=0), facecolors='none', edgecolors='orange', alpha=1)
plt.xlabel('Membrane voltage (mV)')
plt.ylabel('F/F0')
plt.tight_layout()

plt.figure(figsize=(4, 4), num='di1-ANNINE-6plus ')
ax = plt.subplot(212)
for m in Dye1_CELL_ID:
    AVG__ = []
    for i in range(len(AVERAGE_VOLTAGE_STEP_VALUES)):
        if ID__LIST[i] == m:
            BSL__ = np.nanmedian(__IMAGING_DATA_LIST__[i][0:126])
            SIG__ = __IMAGING_DATA_LIST__[i]/BSL__
            plt.plot(np.linspace(0, len(SIG__)*Camera_SAMPLING_INTERVAL, len(SIG__)), SIG__, alpha=0.1, color='orange', linewidth=1)
            if len(SIG__) == SAMPLE_NUMBER_PER_TRIAL:
                AVG__.append(SIG__)
    MEAN = np.nanmean(AVG__, axis=0)
    SEM = np.std(AVG__, axis=0)
    X = np.linspace(0, len(np.nanmean(AVG__, axis=0))*Camera_SAMPLING_INTERVAL, len(np.nanmean(AVG__, axis=0)))
    plt.plot(np.linspace(0, len(np.nanmean(AVG__, axis=0))*Camera_SAMPLING_INTERVAL, len(np.nanmean(AVG__, axis=0))), np.nanmean(AVG__, axis=0), color='orange', alpha=0.8)
    plt.fill_between(X, MEAN+SEM, MEAN-SEM, color='orange', alpha=0.2)

ax2 = plt.subplot(211, sharex=ax)
for m in Dye1_CELL_ID:
    AVG__ = []
    TEMP_RAW = []
    for i in range(len(__RAW_ePhy_LIST__)):
        if ID__LIST[i] == m:
            ax2.plot(np.linspace(-ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])/ (sampling__) -ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])), __RAW_ePhy_LIST__[i], color='orange', alpha=0.2, linewidth=1)
            TEMP_RAW.append(__RAW_ePhy_LIST__[i])
    MEAN = np.nanmean(TEMP_RAW, axis=0)
    SEM = np.std(TEMP_RAW, axis=0)
    X = np.linspace(-ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])/ (sampling__) -ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i]))
    ax2.fill_between(X, MEAN+SEM, MEAN-SEM, color='orange', alpha=0.2)
    ax2.plot(X, MEAN, color='orange', alpha=0.8)
ax2.set_xlabel('Time(sec.)')
ax2.set_ylabel('Membrane voltage(mV)')
ax.set_xlabel('Membrane voltage (mV)')
ax.set_ylabel('F/F0')
ax.set_xlim(2.5, 8)
ax.set_ylim(0.6, 1.6)
ax2.set_ylim(-400, 400)
plt.tight_layout()


plt.figure(figsize=(4, 4), num='di4-ANNINE-6plus ')
ax = plt.subplot(212)
for m in Dye4_CELL_ID:
    AVG__ = []
    for i in range(len(AVERAGE_VOLTAGE_STEP_VALUES)):
        if ID__LIST[i] == m:
            BSL__ = np.nanmedian(__IMAGING_DATA_LIST__[i][0:126])
            SIG__ = __IMAGING_DATA_LIST__[i]/BSL__
            ax.plot(np.linspace(0, len(SIG__)*Camera_SAMPLING_INTERVAL, len(SIG__)), SIG__, alpha=0.1, color='purple', linewidth=1)
            if len(SIG__) == SAMPLE_NUMBER_PER_TRIAL:
                AVG__.append(SIG__)
    MEAN = np.nanmean(AVG__, axis=0)
    SEM = np.std(AVG__, axis=0)
    X = np.linspace(0, len(np.nanmean(AVG__, axis=0))*Camera_SAMPLING_INTERVAL, len(np.nanmean(AVG__, axis=0)))
    ax.plot(np.linspace(0, len(np.nanmean(AVG__, axis=0))*Camera_SAMPLING_INTERVAL, len(np.nanmean(AVG__, axis=0))), np.nanmean(AVG__, axis=0), color='purple', alpha=0.8)
    ax.fill_between(X, MEAN+SEM, MEAN-SEM, color='purple', alpha=0.2)

ax2 = plt.subplot(211, sharex=ax)
for m in Dye1_CELL_ID:
    AVG__ = []
    TEMP_RAW = []
    for i in range(len(__RAW_ePhy_LIST__)):
        if ID__LIST[i] == m:
            ax2.plot(np.linspace(-ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])/ (sampling__) -ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])), __RAW_ePhy_LIST__[i], color='purple', alpha=0.2, linewidth=1)
            TEMP_RAW.append(__RAW_ePhy_LIST__[i])
    MEAN = np.nanmean(TEMP_RAW, axis=0)
    SEM = np.std(TEMP_RAW, axis=0)
    X = np.linspace(-ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i])/ (sampling__) -ePhy_Imaging_DELAY__, len(__RAW_ePhy_LIST__[i]))
    ax2.fill_between(X, MEAN+SEM, MEAN-SEM, color='purple', alpha=0.2)
    ax2.plot(X, MEAN, color='purple', alpha=0.8)
ax2.set_xlabel('Time(sec.)')
ax2.set_ylabel('Membrane voltage(mV)')
ax.set_xlabel('Time (sec.)')
ax.set_ylabel('F/F0')
ax.set_xlim(2.5, 8)
ax.set_ylim(0.6, 1.6)
ax2.set_ylim(-400, 400)
plt.tight_layout()

#AVERAGE C.I. figure
AVG__X = []
AVG__Y = []
plt.figure(num='All observations CI', figsize=(4, 4))
for i in range(len(AVERAGE_fluo_STEP_VALUES)):
    if (ID__LIST[i] in [1, 2, 3]) == False:
        AVG__Y = np.append(AVG__Y, AVERAGE_fluo_STEP_VALUES[i])
        AVG__X = np.append(AVG__X, AVERAGE_VOLTAGE_STEP_VALUES[i])
x = AVG__X
y = AVG__Y
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
fit = p(x)
c_x = np.linspace(-400, 400, 200)
c_y = c_x*z[0]+z[1]
p_y = z[0] * x + z[1]
y_err = y -p_y
p_x = np.arange(np.min(x), np.max(x)+1, 1)
mean_x = np.mean(x)         # mean of x
n = len(x)              # number of samples in origional fit
t = 9               # appropriate t value (where n=9, two tailed 95%)
s_err = np.sum(np.power(y_err, 2))   # sum of the squares of the residuals

confs = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x), 2)/((np.sum(np.power(x, 2)))-n*(np.power(mean_x, 2))))))

#predict y based on test x-values
p_y = z[0]*p_x+z[1]
lower = p_y - abs(confs)
upper = p_y + abs(confs)
print('SLOPE:'+str(z[0]))
plt.xlabel('Membrane voltage (mV)')
plt.ylabel('F/F0')
plt.scatter(x, y, color='orange', alpha=0.8, s=1)
plt.plot(c_x, c_y, color='orange',label='Regression line')
plt.plot(p_x, lower, color='orange', lw=0.5, ls='--',)
plt.plot(p_x, upper, color='orange', ls='--', lw=0.5)
plt.fill_between(p_x, lower, upper,  color='orange', alpha=0.05)
plt.show()



#AVERAGE C.I. figure
AVG__X = []
AVG__Y = []
plt.figure(num='All observations CI', figsize=(4, 4))
for i in range(len(AVERAGE_fluo_STEP_VALUES)):
    if (ID__LIST[i] in [4, 5]) == False:
        AVG__Y = np.append(AVG__Y, AVERAGE_fluo_STEP_VALUES[i])
        AVG__X = np.append(AVG__X, AVERAGE_VOLTAGE_STEP_VALUES[i])
       
x = AVG__X
y = AVG__Y
z = np.polyfit(x, y, 1)
p = np.poly1d(z)
fit = p(x)
c_x = np.linspace(-400, 400, 200)
c_y = c_x*z[0]+z[1]
p_y = z[0] * x + z[1]
y_err = y -p_y
p_x = np.arange(np.min(x), np.max(x)+1, 1)
mean_x = np.mean(x)         # mean of x
n = len(x)              # number of samples in origional fit
t = 9               # appropriate t value (where n=9, two tailed 95%)
s_err = np.sum(np.power(y_err, 2))   # sum of the squares of the residuals

confs = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x), 2)/
            ((np.sum(np.power(x, 2)))-n*(np.power(mean_x, 2))))))

#predict y based on test x-values
p_y = z[0]*p_x+z[1]
lower = p_y - abs(confs)
upper = p_y + abs(confs)
print('SLOPE:'+str(z[0]))
plt.xlabel('Membrane voltage (mV)')
plt.ylabel('F/F0')
plt.scatter(x, y, color='purple', alpha=0.8, s=1)
plt.plot(c_x, c_y, color='purple',label='Regression line')
plt.plot(p_x, lower, color='purple', lw=0.5, ls='--',)
plt.plot(p_x, upper, color='purple', ls='--', lw=0.5)
plt.fill_between(p_x, lower, upper,  color='purple', alpha=0.05)
plt.show()
