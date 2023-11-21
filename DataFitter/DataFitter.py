# Copyright 2014-2021 Keysight Technologies
import InstrumentDriver
from InstrumentConfig import InstrumentQuantity
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import random
import os
from scipy.fftpack import fft, fftfreq, fftshift

class Error(Exception):
    pass

class Driver(InstrumentDriver.InstrumentWorker):
    """ This class implements a demodulation driver"""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection"""
        #global here
        self.fitPerformed = False


    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation"""
        pass


    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function should
        return the actual value set by the instrument"""
        # set to false if called
        self.fitPerformed = False
        return value

    def performGetValue(self, quant, options={}):
        """Check what kind of fit to perform here then return the fit parameters"""
        xUnit = self.getValue('X-axis Unit')
        if xUnit == 'Frequency [Hz]':
            funcType = self.getValue('Function-Hz')
            if funcType == 'Single Lorentzian':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-Hz')
                        if waveformIn['shape'][0] == 0: #Check if waveform in is empty, if so return
                            value = quant.getValue()
                            return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dF = waveformIn['dt']
                        self.fInitial = waveformIn['t0']
                        f_temp = np.arange(self.fInitial, float(len(self.waveformInYValues))*self.dF, self.dF)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_lorentzian(waveformIn)
                        self.fitSignal =  np.array(self.lorentzian(f_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2], self.fitParams[0][3])).copy()
                        self.setValue('Fit Parameters-F0-Hz', self.fitParams[0][0])
                        self.setValue('Fit Parameters-F0 Std-Hz', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Width-Hz', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Width Std-Hz', np.sqrt(self.fitParams[1][2][2]))
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][3])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][3][3]))
                        self.setValue('Fit Parameters-Ql-Hz', self.fitParams[0][0]/self.fitParams[0][2])
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-Hz'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.fInitial, dt=self.dF)
                        return trace
                    elif quant.name.endswith('Residuals-Hz'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.fInitial,dt=self.dF)
                        return trace
                    elif quant.name.endswith('F0-Hz'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('F0 Std-Hz'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Width-Hz'):
                        value = self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Ql-Hz'):
                        value = self.fitParams[0][0]/self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Width Std-Hz'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][3]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][3][3])
                        return value
                else:
                    value = quant.getValue()
                    return value

            elif funcType == 'Single Gaussian':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-Hz')
                        if waveformIn['shape'][0] == 0: #Check if waveform in is empty, if so return
                            value = quant.getValue()
                            return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dF = waveformIn['dt']
                        self.fInitial = waveformIn['t0']
                        f_temp = np.arange(self.fInitial, float(len(self.waveformInYValues))*self.dF, self.dF)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_gaussian(waveformIn)
                        self.fitSignal =  np.array(self.gaussian(f_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2], self.fitParams[0][3])).copy()
                        self.setValue('Fit Parameters-Mean-Hz', self.fitParams[0][0])
                        self.setValue('Fit Parameters-Mean Std-Hz', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Sigma-Hz', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Sigma Std-Hz', np.sqrt(self.fitParams[1][2][2]))
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][3])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][3][3]))
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-Hz'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.fInitial, dt=self.dF)
                        return trace
                    elif quant.name.endswith('Residuals-Hz'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.fInitial,dt=self.dF)
                        return trace
                    elif quant.name.endswith('Mean-Hz'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('Mean Std-Hz'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Sigma-Hz'):
                        value = self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Sigma Std-Hz'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][3]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][3][3])
                        return value
                else:
                    value = quant.getValue()
                    return value

            elif funcType == 'Exponential':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-Hz')
                        if waveformIn['shape'][0] == 0: #Check if waveform in is empty, if so return
                            value = quant.getValue()
                            return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dF = waveformIn['dt']
                        self.fInitial = waveformIn['t0']
                        f_temp = np.arange(self.fInitial, float(len(self.waveformInYValues))*self.dF, self.dF)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_exp_lin(waveformIn)
                        self.fitSignal =  np.array(self.exponential(f_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2])).copy()
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][0])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Decay Constant-Hz', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Decay Constant Std-Hz', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][2][2]))
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-Hz'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.fInitial, dt=self.dF)
                        return trace
                    elif quant.name.endswith('Residuals-Hz'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.fInitial,dt=self.dF)
                        return trace
                    elif quant.name.endswith('Decay Constant-Hz'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Decay Constant Std-Hz'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                else:
                    value = quant.getValue()
                    return value

            elif funcType == 'Exponential W/ Sinusoid':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-Hz')
                        if waveformIn['shape'][0] == 0: #Check i waveform in is empty, if so return
                            value = quant.getValue()
                            return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dF = waveformIn['dt']
                        self.fInitial = waveformIn['t0']
                        f_temp = np.arange(self.fInitial, float(len(self.waveformInYValues))*self.dF, self.dF)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_exp_lin_with_sinusoid(waveformIn)
                        self.fitSignal =  np.array(self.exp_with_sinusoid(f_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2], self.fitParams[0][3], self.fitParams[0][4])).copy()
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][0])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Decay Constant-Hz', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Decay Constant Std-Hz', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Detuning-Hz', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Detuning Std-Hz', np.sqrt(self.fitParams[1][2][2]))
                        self.setValue('Fit Parameters-Phase Offset-Hz', self.fitParams[0][3]*180.0/np.pi)
                        self.setValue('Fit Parameters-Phase Offset Std-Hz', np.sqrt(self.fitParams[1][3][3])*180.0/np.pi)
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][4])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][4][4]))
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-Hz'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.fInitial, dt=self.dF)
                        return trace
                    elif quant.name.endswith('Residuals-Hz'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.fInitial,dt=self.dF)
                        return trace
                    elif quant.name.endswith('Decay Constant-Hz'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Decay Constant Std-Hz'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Detuning-Hz'):
                        value = self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Detuning Std-Hz'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                    elif quant.name.endswith('Phase Offset-Hz'):
                        value = 180.0/np.pi*self.fitParams[0][3]
                        return value
                    elif quant.name.endswith('Phase Offset Std-Hz'):
                        value = 180.0/np.pi*np.sqrt(self.fitParams[1][3][3])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][4]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][4][4])
                        return value
                else:
                    value = quant.getValue()
                    return value

        if xUnit == 'Time [s]':
            funcType = self.getValue('Function-s')
            if funcType == 'Single Lorentzian':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-s')
                        if waveformIn['shape'][0] == 0: #Check if waveform in is empty, if so return
                                value = quant.getValue()
                                return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dt = waveformIn['dt']
                        self.tInitial = waveformIn['t0']
                        t_temp = np.arange(self.tInitial, float(len(self.waveformInYValues))*self.dt, self.dt)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_lorentzian(waveformIn)
                        self.fitSignal =  self.lorentzian(t_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2], self.fitParams[0][3])
                        self.setValue('Fit Parameters-F0-s', self.fitParams[0][0])
                        self.setValue('Fit Parameters-F0 Std-s', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Width-s', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Width Std-s', np.sqrt(self.fitParams[1][2][2]))
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][3])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][3][3]))
                        self.setValue('Fit Parameters-Ql-s', self.fitParams[0][0]/self.fitParams[0][2])
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-s'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.tInitial, dt=self.dt)
                        return trace
                    elif quant.name.endswith('Residuals-s'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.tInitial,dt=self.dt)
                        return trace
                    elif quant.name.endswith('F0-s'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('F0 Std-s'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Width-s'):
                        value = self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Ql-s'):
                        value = self.fitParams[0][0]/self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Width Std-s'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][3]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][3][3])
                        return value
                else:
                    value = quant.getValue()
                    return value

            elif funcType == 'Single Gaussian':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-s')
                        if waveformIn['shape'][0] == 0: #Check if waveform in is empty, if so return
                            value = quant.getValue()
                            return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dT = waveformIn['dt']
                        self.tInitial = waveformIn['t0']
                        f_temp = np.arange(self.tInitial, float(len(self.waveformInYValues))*self.dT, self.dT)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_gaussian(waveformIn)
                        self.fitSignal =  np.array(self.gaussian(f_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2], self.fitParams[0][3])).copy()
                        self.setValue('Fit Parameters-Mean-s', self.fitParams[0][0])
                        self.setValue('Fit Parameters-Mean Std-s', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Sigma-s', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Sigma Std-s', np.sqrt(self.fitParams[1][2][2]))
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][3])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][3][3]))
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-s'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.tInitial, dt=self.dT)
                        return trace
                    elif quant.name.endswith('Residuals-s'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.tInitial,dt=self.dT)
                        return trace
                    elif quant.name.endswith('Mean-s'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('Mean Std-s'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Sigma-s'):
                        value = self.fitParams[0][2]
                        return value
                        return value
                    elif quant.name.endswith('Sigma Std-s'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][3]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][3][3])
                        return value
                else:
                    value = quant.getValue()
                    return value

            elif funcType == 'Exponential':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-s')
                        if waveformIn['shape'][0] == 0: #Check i waveform in is empty, if so return
                            value = quant.getValue()
                            return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dT = waveformIn['dt']
                        self.tInitial = waveformIn['t0']
                        t_temp = np.arange(self.tInitial, float(len(self.waveformInYValues))*self.dT, self.dT)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_exp_lin(waveformIn)
                        self.fitSignal =  np.array(self.exponential(t_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2])).copy()
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][0])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Decay Constant-s', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Decay Constant Std-s', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][2][2]))
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-s'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.tInitial, dt=self.dT)
                        return trace
                    elif quant.name.endswith('Residuals-s'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.tInitial,dt=self.dT)
                        return trace
                    elif quant.name.endswith('Decay Constant-s'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Decay Constant Std-s'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                else:
                    value = quant.getValue()
                    return value

            elif funcType == 'Exponential W/ Sinusoid':
                if quant.name.startswith('Fit Parameters-'):
                    if self.fitPerformed == False:
                        waveformIn = self.getValue('Waveform In-s')
                        if waveformIn['shape'][0] == 0: #Check i waveform in is empty, if so return
                            value = quant.getValue()
                            return value
                        #use values from incoming waveform
                        self.waveformInYValues = np.array(waveformIn['y']).copy()
                        self.dT = waveformIn['dt']
                        self.tInitial = waveformIn['t0']
                        t_temp = np.arange(self.tInitial, float(len(self.waveformInYValues))*self.dT, self.dT)
                        #Fit to incoming waveform and create an array based on the fit
                        self.fitParams = self.fit_single_exp_lin_with_sinusoid(waveformIn)
                        self.fitSignal =  np.array(self.exp_with_sinusoid(t_temp, self.fitParams[0][0], self.fitParams[0][1], self.fitParams[0][2], self.fitParams[0][3], self.fitParams[0][4])).copy()
                        self.setValue('Fit Parameters-Amplitude', self.fitParams[0][0])
                        self.setValue('Fit Parameters-Amplitude Std', np.sqrt(self.fitParams[1][0][0]))
                        self.setValue('Fit Parameters-Decay Constant-s', self.fitParams[0][1])
                        self.setValue('Fit Parameters-Decay Constant Std-s', np.sqrt(self.fitParams[1][1][1]))
                        self.setValue('Fit Parameters-Detuning-s', self.fitParams[0][2])
                        self.setValue('Fit Parameters-Detuning Std-s', np.sqrt(self.fitParams[1][2][2]))
                        self.setValue('Fit Parameters-Phase Offset-s', self.fitParams[0][3]*180.0/np.pi)
                        self.setValue('Fit Parameters-Phase Offset Std-s', np.sqrt(self.fitParams[1][3][3])*180.0/np.pi)
                        self.setValue('Fit Parameters-Offset', self.fitParams[0][4])
                        self.setValue('Fit Parameters-Offset Std', np.sqrt(self.fitParams[1][4][4]))
                        self.fitPerformed = True
                    else:
                        self.fitPerformed = True
                    if quant.name.endswith('Waveform-s'):#return waveform based on fit parameters
                        trace = quant.getTraceDict(self.fitSignal, t0=self.tInitial, dt=self.dT)
                        return trace
                    elif quant.name.endswith('Residuals-s'):#return residuals
                        trace = quant.getTraceDict(self.fitParams[2], t0=self.tInitial,dt=self.dT)
                        return trace
                    elif quant.name.endswith('Decay Constant-s'):
                        value = self.fitParams[0][1]
                        return value
                    elif quant.name.endswith('Decay Constant Std-s'):
                        value = np.sqrt(self.fitParams[1][1][1])
                        return value
                    elif quant.name.endswith('Detuning-s'):
                        value = self.fitParams[0][2]
                        return value
                    elif quant.name.endswith('Detuning Std-s'):
                        value = np.sqrt(self.fitParams[1][2][2])
                        return value
                    elif quant.name.endswith('Phase Offset-s'):
                        value = 180.0/np.pi*self.fitParams[0][3]
                        return value
                    elif quant.name.endswith('Phase Offset Std-s'):
                        value = 180.0/np.pi*np.sqrt(self.fitParams[1][3][3])
                        return value
                    elif quant.name.endswith('Amplitude'):
                        value = self.fitParams[0][0]
                        return value
                    elif quant.name.endswith('Amplitude Std'):
                        value = np.sqrt(self.fitParams[1][0][0])
                        return value
                    elif quant.name.endswith('Offset'):
                        value = self.fitParams[0][4]
                        return value
                    elif quant.name.endswith('Offset Std'):
                        value = np.sqrt(self.fitParams[1][4][4])
                        return value
                else:
                    value = quant.getValue()
                    return value

        else:
            value = quant.getValue()
            return value

    def lorentzian(self, x, x0, a, gam, offset):
        return a * gam**2 / ( gam**2 + ( x - x0 )**2)+offset

    def fit_single_lorentzian(self, amps):
        amplitudes = np.array(amps['y']).copy()
        blamplitudes = np.array(amps['y']).copy()
        amps_i = amplitudes.mean()
        dF = amps['dt']
        fInitial = amps['t0']
        f_temp = np.arange(fInitial, float(len(amplitudes))*dF, dF)

        amp_sign = np.ones(4)
        filter_window = 51 #must be an odd integer
        guess = np.zeros(4) #guess array format: f0, amplitude, linewidth, offset.
        ####Check user amplitude here
        if self.getValue('Guess Amplitude?'):
            guess[1] = self.getValue('Amplitude Guess')

            if guess[1] < 0.0:
                amplitudes *= -1.0
                amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
        else: #Use Second derivative to check if band stop and if so then invert, else do nothing
            secondDerivative = savgol_filter(amplitudes, filter_window, 4, 2)
            maxPos = np.abs(secondDerivative.max())
            maxNeg = np.abs(secondDerivative.min())
            if maxPos > maxNeg:
                amplitudes *= -1.0
                amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])

            else:
                amplitudes *= 1.0
                amp_sign *= np.ones(4)


        #Check user frequency here
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess F0-Hz?'):
                guess[0] = self.getValue('F0 Guess-Hz')
            else:
                #Use a filter to the raw data to better estimate the resonance frequency
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                f0_guess_index = sgFilter.argmax()

                #From trial and error some fixes to the guess f0
                if (f0_guess_index == len(f_temp)-1):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()
                elif (f0_guess_index == 0):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()

                #Input index into frequency array to have an intial resonance frequency guess
                guess[0] = f_temp[f0_guess_index]

        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess F0-s?'):
                guess[0] = self.getValue('F0 Guess-s')
            else:
                #Use a filter to the raw data to better estimate the resonance frequency
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                f0_guess_index = sgFilter.argmax()

                #From trial and error some fixes to the guess f0
                if (f0_guess_index == len(f_temp)-1):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()
                elif (f0_guess_index == 0):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()

                #Input index into frequency array to have an intial resonance frequency guess
                guess[0] = f_temp[f0_guess_index]

        #Check user offset guess here
        if self.getValue('Guess Offset?'):
            guess[3] = self.getValue('Offset Guess')
        #### Check if any component is negative, if so add an offset to ensure it the waveform is always positive
        offset_correction = np.zeros(4)
        if amplitudes.min() < 0.0:
            offset_correction[3] += -2.0*amplitudes.min()
            amplitudes += offset_correction[3]
        else:
            offset_correction += offset_correction

        #Check user guess width here
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess Width-Hz?'):
                guess[2] = self.getValue('Width Guess-Hz')
            else:
                index_width = np.abs(guess[0]-f_temp).argmin()
                #Come up with an initial guess for the spectral width
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-0.5*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                #Handling use case where it finds a center frequency in the frequency sweep that is wrong and has a wide width
                if (guess[2] > 0.3*f_temp[-1]) and (filter_window != 101):
                    filter_window = 201
                    amplitudes -= offset_correction[3]
                    amplitudes *= amp_sign[3]
                    #amplitudes = (amplitudes-offset_correction[3])*amp_sign[3]
                    secondDerivative = savgol_filter(amplitudes, filter_window, 4, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])

                    else:
                        amplitudes = np.array(amps['y']).copy() #Recheck this
                        amp_sign *= np.ones(4)

                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    index_width = sgFilter.argmax()

                    #### Check if any component is negative, if so add an offset to ensure it the waveform is always positive
                    if amplitudes.min() < 0.0:
                        offset_correction[3] += -2.0*amplitudes.min()
                        amplitudes += offset_correction[3]
                    else:
                        offset_correction += offset_correction

                    #Come up with an initial guess for the spectral width ###################Recheck the absolute value
                    guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-0.5*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                else:
                    guess[2] *= 1.0
        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess Width-s?'):
                guess[2] = self.getValue('Width Guess-s')
            else:
                index_width = np.abs(guess[0]-f_temp).argmin()
                #Come up with an initial guess for the spectral width
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-0.5*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                #Handling use case where it finds a center frequency in the frequency sweep that is wrong and has a wide width
                if (guess[2] > 0.3*f_temp[-1]) and (filter_window != 101):
                    filter_window = 201
                    amplitudes -= offset_correction[3]
                    amplitudes *= amp_sign[3]
                    #amplitudes = (amplitudes-offset_correction[3])*amp_sign[3]
                    secondDerivative = savgol_filter(amplitudes, filter_window, 4, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])

                    else:
                        amplitudes = np.array(amps['y']).copy() #Recheck this
                        amp_sign *= np.ones(4)

                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    index_width = sgFilter.argmax()

                    #### Check if any component is negative, if so add an offset to ensure it the waveform is always positive
                    if amplitudes.min() < 0.0:
                        offset_correction[3] += -2.0*amplitudes.min()
                        amplitudes += offset_correction[3]
                    else:
                        offset_correction += offset_correction

                    #Come up with an initial guess for the spectral width ###################Recheck the absolute value
                    guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-0.5*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                else:
                    guess[2] *= 1.0

        #Come up with initial guess of the offset
        if guess[3] == 0.0:
            guess[3] = amplitudes.mean()

        #Come up with initial guess of the amplitude.  Seems the Raw data gives a better fit
        if guess[1] == 0:
            index = np.abs(f_temp-guess[0]).argmin()
            guess[1] = amplitudes[index]-guess[3]



        #bounds for fitting
        bounds = np.array([ [f_temp[0], -np.inf, 0.0, -np.inf], [f_temp[-1], np.inf, f_temp[-1], np.inf]])

        try:
            popt, pcov = curve_fit(self.lorentzian, f_temp, amplitudes, guess, bounds=bounds)
        except RuntimeError:
            self.log("Error with fit for Guess settings {}".format(guess), level=30)
            popt, pcov = np.array([[1.0, 1.0, 1.0, 1.0],[[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0]]])

        amps_f = ((amplitudes-offset_correction[3])*amp_sign[1]).mean()
        popt -= offset_correction
        popt *= amp_sign
        if (amps_f - amps_i) > 1e-12:
            amp_index = (np.abs(popt[0]-f_temp)).argmin()
            guess = np.array([popt[0], blamplitudes[amp_index], popt[2], blamplitudes.mean()])
            popt, pcov = curve_fit(self.lorentzian, f_temp, blamplitudes, guess, bounds=bounds)

        residuals = (blamplitudes - self.lorentzian(f_temp, popt[0], popt[1], popt[2], popt[3]))
        if np.abs(residuals.mean()) > 1.0e-11:
            if amp_sign[1] < 0: #false flip
                amplitudes = np.array(blamplitudes).copy()
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[0] = f_temp[sgFilter.argmax()]
                guess[1] = amplitudes.max()-amplitudes.mean()
                index_width = np.abs(guess[0]-f_temp).argmin()
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-0.5*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                guess[3] = amplitudes.mean()
                popt, pcov = curve_fit(self.lorentzian, f_temp, amplitudes, guess, bounds=bounds)
                popt[3] += offset_correction

            else: #missed flip
                amplitudes = np.array(blamplitudes).copy()
                amplitudes *= -1.0
                amp_sign[1] *= -1.0
                amp_sign[3] *= -1.0
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[0] = f_temp[sgFilter.argmax()]
                guess[1] = amplitudes.max()-amplitudes.mean()
                index_width = np.abs(guess[0]-f_temp).argmin()
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-0.5*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                guess[3] = amplitudes.mean()
                popt, pcov = curve_fit(self.lorentzian, f_temp, amplitudes, guess, bounds=bounds)
                popt[3] += offset_correction
                popt[3] *= amp_sign[3]
                popt[1] *= amp_sign[1]

        residuals = (blamplitudes - self.lorentzian(f_temp, popt[0], popt[1], popt[2], popt[3]))

        if np.abs(residuals.mean()) > 1.0e-11:
            if amp_sign[1] < 0: #missed flip
                amplitudes = np.array(blamplitudes).copy()
                amplitudes *= -1.0
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[0] = f_temp[sgFilter.argmax()]
                guess[1] = amplitudes.max()-amplitudes.mean()
                index_width = np.abs(guess[0]-f_temp).argmin()
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-0.5*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                guess[3] = amplitudes.mean()
                popt, pcov = curve_fit(self.lorentzian, f_temp, amplitudes, guess, bounds=bounds)
                popt[3] += offset_correction
                popt[3] *= amp_sign[3]
                popt[1] *= amp_sign[1]

        residuals = (blamplitudes - self.lorentzian(f_temp, popt[0], popt[1], popt[2], popt[3]))
        return popt, pcov, residuals

    def gaussian(self, x, mean, amp, sigma, offset):
        return amp*np.exp(-(x-mean)**2/(2*sigma**2))+offset

    def fit_single_gaussian(self, amps):
        amplitudes = np.array(amps['y']).copy()
        blamplitudes = np.array(amps['y']).copy()
        amps_i = amplitudes.mean()
        dF = amps['dt']
        fInitial = amps['t0']
        f_temp = np.arange(fInitial, float(len(amplitudes))*dF, dF)

        amp_sign = np.ones(4)
        filter_window = 51 #must be an odd integer
        sgFilter = []
        guess = np.zeros(4) #guess array format: f0, amplitude, linewidth, offset.
        ####Check user amplitude here
        if self.getValue('Guess Amplitude?'):
            guess[1] = self.getValue('Amplitude Guess')

            if guess[1] < 0.0:
                amplitudes *= -1.0
                amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
        else: #Use Second derivative to check if band stop and if so then invert, else do nothing
            secondDerivative = savgol_filter(amplitudes, filter_window, 4, 2)
            maxPos = np.abs(secondDerivative.max())
            maxNeg = np.abs(secondDerivative.min())
            if maxPos > maxNeg:
                amplitudes *= -1.0
                amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])

            else:
                amplitudes *= 1.0
                amp_sign *= np.ones(4)


        #Check user frequency here
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess Mean-Hz?'):
                guess[0] = self.getValue('Mean Guess-Hz')
            else:
                #Use a filter to the raw data to better estimate the resonance frequency
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                f0_guess_index = sgFilter.argmax()

                #From trial and error some fixes to the guess f0
                if (f0_guess_index == len(f_temp)-1):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()
                elif (f0_guess_index == 0):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()

                #Input index into frequency array to have an intial resonance frequency guess
                guess[0] = f_temp[f0_guess_index]

        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess Mean-s?'):
                guess[0] = self.getValue('Mean Guess-s')
            else:
                #Use a filter to the raw data to better estimate the resonance frequency
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                f0_guess_index = sgFilter.argmax()

                #From trial and error some fixes to the guess f0
                if (f0_guess_index == len(f_temp)-1):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()
                elif (f0_guess_index == 0):
                    filter_window = 101
                    secondDerivative = savgol_filter(amplitudes, filter_window, 2, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])
                    else:
                        amplitudes *= 1.0
                        amp_sign *= np.ones(4)
                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    f0_guess_index = sgFilter.argmax()

                #Input index into frequency array to have an intial resonance frequency guess
                guess[0] = f_temp[f0_guess_index]

        #Check user offset guess here
        if self.getValue('Guess Offset?'):
            guess[3] = self.getValue('Offset Guess')
        #### Check if any component is negative, if so add an offset to ensure it the waveform is always positive
        offset_correction = np.zeros(4)
        if amplitudes.min() < 0.0:
            offset_correction[3] += -2.0*amplitudes.min()
            amplitudes += offset_correction[3]
        else:
            offset_correction += offset_correction

        #Check user guess width here
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess Sigma-Hz?'):
                guess[2] = self.getValue('Sigma Guess-Hz')
            else:
                index_width = np.abs(guess[0]-f_temp).argmin()
                #Come up with an initial guess for the spectral width
                if len(sgFilter) == 0:
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-np.exp(-1)*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                #Handling use case where it finds a center frequency in the frequency sweep that is wrong and has a wide width
                if (guess[2] > 0.3*f_temp[-1]) and (filter_window != 101):
                    filter_window = 201
                    amplitudes -= offset_correction[3]
                    amplitudes *= amp_sign[3]
                    #amplitudes = (amplitudes-offset_correction[3])*amp_sign[3]
                    secondDerivative = savgol_filter(amplitudes, filter_window, 4, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])

                    else:
                        amplitudes = np.array(amps['y']).copy() #Recheck this
                        amp_sign *= np.ones(4)

                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    index_width = sgFilter.argmax()

                    #### Check if any component is negative, if so add an offset to ensure it the waveform is always positive
                    if amplitudes.min() < 0.0:
                        offset_correction[3] += -2.0*amplitudes.min()
                        amplitudes += offset_correction[3]
                    else:
                        offset_correction += offset_correction

                    #Come up with an initial guess for the spectral width ###################Recheck the absolute value
                    guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-np.exp(-1)*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                else:
                    guess[2] *= 1.0
        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess Sigma-s?'):
                guess[2] = self.getValue('Sigma Guess-s')
            else:
                index_width = np.abs(guess[0]-f_temp).argmin()
                #Come up with an initial guess for the spectral width
                if len(sgFilter) == 0:
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-np.exp(-1)*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                #Handling use case where it finds a center frequency in the frequency sweep that is wrong and has a wide width
                if (guess[2] > 0.3*f_temp[-1]) and (filter_window != 101):
                    filter_window = 201
                    amplitudes -= offset_correction[3]
                    amplitudes *= amp_sign[3]
                    #amplitudes = (amplitudes-offset_correction[3])*amp_sign[3]
                    secondDerivative = savgol_filter(amplitudes, filter_window, 4, 2)
                    maxPos = np.abs(secondDerivative.max())
                    maxNeg = np.abs(secondDerivative.min())
                    if maxPos > maxNeg:
                        amplitudes *= -1.0
                        amp_sign *= np.array([1.0, -1.0, 1.0, -1.0])

                    else:
                        amplitudes = np.array(amps['y']).copy() #Recheck this
                        amp_sign *= np.ones(4)

                    #Use a filter to the raw data to better estimate the resonance frequency
                    sgFilter = savgol_filter(amplitudes, filter_window, 2)
                    index_width = sgFilter.argmax()

                    #### Check if any component is negative, if so add an offset to ensure it the waveform is always positive
                    if amplitudes.min() < 0.0:
                        offset_correction[3] += -2.0*amplitudes.min()
                        amplitudes += offset_correction[3]
                    else:
                        offset_correction += offset_correction

                    #Come up with an initial guess for the spectral width ###################Recheck the absolute value
                    guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-np.exp(-1)*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                else:
                    guess[2] *= 1.0

        #Come up with initial guess of the offset
        if guess[3] == 0.0:
            guess[3] = amplitudes.mean()

        #Come up with initial guess of the amplitude.  Seems the Raw data gives a better fit
        if guess[1] == 0:
            index = np.abs(f_temp-guess[0]).argmin()
            guess[1] = amplitudes[index]-guess[3]

        #occasionally filter misses when varying widths from very narrow to very large this catches it
        amp_temp = self.getValue('Guess Amplitude?')
        if amp_temp == False:
            if guess[1] < 0.0:
                guess[1] *= -1.0
                guess[0] = f_temp[amplitudes.argmax()]

        #bounds for fitting
        bounds = np.array([ [f_temp[0], -np.inf, 0.0, -np.inf], [f_temp[-1], np.inf, f_temp[-1], np.inf]])

        try:
            popt, pcov = curve_fit(self.gaussian, f_temp, amplitudes, guess, bounds=bounds)
        except RuntimeError:
            self.log("Error with fit for Guess settings {}".format(guess), level=30)
            popt, pcov = np.array([[1.0, 1.0, 1.0, 1.0],[[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0]]])

        amps_f = ((amplitudes-offset_correction[3])*amp_sign[1]).mean()
        popt -= offset_correction
        popt *= amp_sign
        if (amps_f - amps_i) > 1e-12:
            amp_index = (np.abs(popt[0]-f_temp)).argmin()
            guess = np.array([popt[0], blamplitudes[amp_index], popt[2], blamplitudes.mean()])
            popt, pcov = curve_fit(self.gaussian, f_temp, blamplitudes, guess, bounds=bounds)

        residuals = (blamplitudes - self.gaussian(f_temp, popt[0], popt[1], popt[2], popt[3]))

        if np.abs(residuals.mean()) > 1.0e-11:
            if amp_sign[1] < 0: #false flip
                amplitudes = np.array(blamplitudes).copy()
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[0] = f_temp[sgFilter.argmax()]
                guess[1] = amplitudes.max()-amplitudes.mean()
                index_width = np.abs(guess[0]-f_temp).argmin()
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-np.exp(-1)*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                guess[3] = amplitudes.mean()
                popt, pcov = curve_fit(self.gaussian, f_temp, amplitudes, guess, bounds=bounds)
                popt[3] += offset_correction

            else: #missed flip
                amplitudes = np.array(blamplitudes).copy()
                amplitudes *= -1.0
                amp_sign[1] *= -1.0
                amp_sign[3] *= -1.0
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[0] = f_temp[sgFilter.argmax()]
                guess[1] = amplitudes.max()-amplitudes.mean()
                index_width = np.abs(guess[0]-f_temp).argmin()
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-np.exp(-1)*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                guess[3] = amplitudes.mean()
                popt, pcov = curve_fit(self.gaussian, f_temp, amplitudes, guess, bounds=bounds)
                popt[3] += offset_correction
                popt[3] *= amp_sign[3]
                popt[1] *= amp_sign[1]

        residuals = (blamplitudes - self.gaussian(f_temp, popt[0], popt[1], popt[2], popt[3]))

        if np.abs(residuals.mean()) > 1.0e-11:
            if amp_sign[1] < 0: #missed flip
                amplitudes = np.array(blamplitudes).copy()
                amplitudes *= -1.0
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
                sgFilter = savgol_filter(amplitudes, filter_window, 2)
                guess[0] = f_temp[sgFilter.argmax()]
                guess[1] = amplitudes.max()-amplitudes.mean()
                index_width = np.abs(guess[0]-f_temp).argmin()
                guess[2] = np.abs(index_width -(np.abs(sgFilter-sgFilter.mean()-np.exp(-1.0)*(sgFilter[index_width]-sgFilter.mean())).argmin()))*dF
                guess[3] = amplitudes.mean()
                popt, pcov = curve_fit(self.gaussian, f_temp, amplitudes, guess, bounds=bounds)
                popt[3] += offset_correction
                popt[3] *= amp_sign[3]
                popt[1] *= amp_sign[1]

        residuals = (blamplitudes - self.gaussian(f_temp, popt[0], popt[1], popt[2], popt[3]))
        return popt, pcov, residuals

    def exponential(self, x, a, tau, offset):
        return a*np.exp(-(x/tau))+offset

    def fit_single_exp_lin(self,amps):
        amplitudes = np.array(amps['y']).copy()
        blamplitudes = np.array(amps['y']).copy()
        amps_i = amplitudes.mean()
        dT = amps['dt']
        tInitial = amps['t0']
        t_temp = np.arange(tInitial, float(len(amplitudes))*dT, dT)
        amp_sign = np.ones(3)
        filter_window = 51 #must be an odd integer
        firstDerivative = []

        guess = np.zeros(3) #guess array format: amp, decay constant, offset.

        ####Check user amplitude here
        if self.getValue('Guess Amplitude?'):
            guess[0] = self.getValue('Amplitude Guess')
        else:
            firstDerivative = savgol_filter(amplitudes, filter_window, 3, 1) #first derivative
            if np.abs(firstDerivative.max()) > np.abs(firstDerivative.min()):
                amplitudes *= -1.0
                amp_sign *= np.array([-1.0, 1.0, -1.0])
            guess[0] = amplitudes.max()-amplitudes.min()

        ####Check user offset guess
        if self.getValue('Guess Offset?'):
            guess[2] = self.getValue('Offset Guess')
        else:
            if amplitudes.min() < 0.0:
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
            else:
                offset_correction = 0.0


            guess[2] = amplitudes.mean()

        ######## Decay constant guess ###############3
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess Decay Constant-Hz?'):
                guess[1] = self.getValue('Decay Constant Guess-Hz')
            else:
                response = np.abs(fftshift(fft(amplitudes*np.sin(2.0*np.pi*t_temp/(10.0*dT))) / (len(amplitudes)/2.0)))
                xf = fftshift(fftfreq(len(amplitudes), dT))
                amps = 1.0*response[len(amplitudes)//2:]
                freqs = 1.0*xf[len(amplitudes)//2:]
                dFreqs = np.abs(freqs[1]-freqs[0])
                index_f0 = np.abs(1.0/(10.0*dT)-freqs).argmin()
                guess_width = np.abs(index_f0 -(np.abs(amps-amps.mean()-0.5*(amps[index_f0]-amps.mean())).argmin()))*(freqs[1]-freqs[0])
                guess_decay = np.array([1.0/(10.0*dT), amps.max()-amps.min(), guess_width, amps.min()])
                temp, temp1 = curve_fit(self.lorentzian, freqs, amps, guess_decay)
                if temp[2] < dFreqs:
                    guess[1] = 1/(2.0*np.pi*dFreqs)
                else:
                    guess[1] = 1/(2.0*np.pi*temp[2])
        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess Decay Constant-s?'):
                guess[1] = self.getValue('Decay Constant Guess-s')
            else:
                response = np.abs(fftshift(fft(amplitudes*np.sin(2.0*np.pi*t_temp/(10.0*dT))) / (len(amplitudes)/2.0)))
                xf = fftshift(fftfreq(len(amplitudes), dT))
                amps = 1.0*response[len(amplitudes)//2:]
                freqs = 1.0*xf[len(amplitudes)//2:]
                dFreqs = np.abs(freqs[1]-freqs[0])
                index_f0 = np.abs(1.0/(10.0*dT)-freqs).argmin()
                guess_width = np.abs(index_f0 -(np.abs(amps-amps.mean()-0.5*(amps[index_f0]-amps.mean())).argmin()))*(freqs[1]-freqs[0])
                guess_decay = np.array([1.0/(10.0*dT), amps.max()-amps.min(), guess_width, amps.min()])
                temp, temp1 = curve_fit(self.lorentzian, freqs, amps, guess_decay)
                if temp[2] < dFreqs:
                    guess[1] = 1/(2.0*np.pi*dFreqs)
                else:
                    guess[1] = 1/(2.0*np.pi*temp[2])

        bounds = np.array([ [-np.inf, t_temp[0], -np.inf], [np.inf, np.inf, np.inf]])

        try:
            popt, pcov = curve_fit(self.exponential, t_temp, amplitudes, guess, bounds=bounds)
        except RuntimeError:
            self.log("Error with fit for Guess settings {}".format(guess), level=30)
            popt, pcov = np.array([[1.0, 1.0, 1.0],[[1.0, 1.0, 1.0],[1.0, 1.0, 1.0],[1.0, 1.0, 1.0],[1.0, 1.0, 1.0]]])

        popt[2] += offset_correction
        popt[2] *= amp_sign[2]
        popt[0] *= amp_sign[0]
        residuals = (blamplitudes - self.exponential(t_temp, popt[0], popt[1], popt[2]))

        return popt, pcov, residuals

    def exp_with_sinusoid(self, x, a, tau, det, phi, offset):
        return a*np.exp(-x/tau)*np.cos(2*np.pi*det*x + phi) + offset

    def fit_single_exp_lin_with_sinusoid(self, amps):
        amplitudes = np.array(amps['y']).copy()
        blamplitudes = np.array(amps['y']).copy()
        amps_i = amplitudes.mean()
        dT = amps['dt']
        tInitial = amps['t0']
        t_temp = np.arange(tInitial, float(len(amplitudes))*dT, dT)
        amp_sign = np.ones(5)
        filter_window = 51 #must be an odd integer
        firstDerivative = []
        self.fitLorGuessParams = []
        self.temp = np.array([])

        guess = np.zeros(5) #guess array format: amp, decay constant, offset.
        ####Check user amplitude here
        if self.getValue('Guess Amplitude?'):
            guess[0] = self.getValue('Amplitude Guess')
        else:
            firstDerivative = savgol_filter(amplitudes, filter_window, 3, 1) #first derivative
            if firstDerivative[0] > 0.0:
                 amplitudes *= -1.0
                 amp_sign *= np.array([-1.0, 1.0, 1.0, 1.0, -1.0])
            guess[0] = amplitudes.max()-amplitudes.min()

        ####Check user offset guess
        if self.getValue('Guess Offset?'):
            guess[4] = self.getValue('Offset Guess')
        else:
            if amplitudes.min() < 0.0:
                offset_correction = 2.0*amplitudes.min()
                amplitudes -= offset_correction
            else:
                offset_correction = 0.0

            guess[4] = amplitudes.mean()



        ######## Decay constant guess ###############3
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess Decay Constant-Hz?'):
                guess[1] = self.getValue('Decay Constant Guess-Hz')
            else:
                response = np.abs(fftshift(fft(amplitudes) / (len(amplitudes)/2.0)))
                xf = fftshift(fftfreq(len(amplitudes), dT))
                amps = 1.0*response[len(amplitudes)//2:]
                amps[0] = 0.0
                freqs = 1.0*xf[len(amplitudes)//2:]
                dFreqs = np.abs(freqs[1]-freqs[0])
                index_f0 = amps.argmax()
                guess_width = np.abs(index_f0 -(np.abs(amps-amps.mean()-0.5*(amps[index_f0]-amps.mean())).argmin()))*dFreqs
                guess_decay = np.array([freqs[index_f0], amps.max()-amps.min(), guess_width, amps.min()])
                self.temp, temp1 = curve_fit(self.lorentzian, freqs, amps, guess_decay)
                guess[1] = 1/(2.0*np.pi*self.temp[2])
                if guess[1] < dT:
                    guess[1] = 5.0*dT
                    self.fitLorGuessParams = []

        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess Decay Constant-s?'):
                guess[1] = self.getValue('Decay Constant Guess-s')
            else:
                response = np.abs(fftshift(fft(amplitudes) / (len(amplitudes)/2.0)))
                xf = fftshift(fftfreq(len(amplitudes), dT))
                amps = 1.0*response[len(amplitudes)//2:]
                amps[0] = 0.0
                freqs = 1.0*xf[len(amplitudes)//2:]
                dFreqs = np.abs(freqs[1]-freqs[0])
                index_f0 = amps.argmax()
                guess_width = np.abs(index_f0 -(np.abs(amps-amps.mean()-0.5*(amps[index_f0]-amps.mean())).argmin()))*dFreqs
                guess_decay = np.array([freqs[index_f0], amps.max()-amps.min(), guess_width, amps.min()])
                self.temp, temp1 = curve_fit(self.lorentzian, freqs, amps, guess_decay, maxfev=10000)
                guess[1] = 1/(2.0*np.pi*self.temp[2])
                if guess[1] < dT:
                    guess[1] = 5.0*dT
                    self.fitLorGuessParams = []

        ######## Detuning guess ###############3
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess Detuning-Hz?'):
                guess[2] = self.getValue('Detuning Guess-Hz')
            else:
                if self.temp.any():
                    guess[2] = self.temp[0]
                else:
                    response = np.abs(fftshift(fft(amplitudes) / (len(amplitudes)/2.0)))
                    xf = fftshift(fftfreq(len(amplitudes), dT))
                    amps = 1.0*response[len(amplitudes)//2:]
                    amps[0] = 0.0
                    freqs = 1.0*xf[len(amplitudes)//2:]
                    dFreqs = np.abs(freqs[1]-freqs[0])
                    index_f0 = amps.argmax()
                    guess_width = np.abs(index_f0 -(np.abs(amps-amps.mean()-0.5*(amps[index_f0]-amps.mean())).argmin()))*dFreqs
                    guess_decay = np.array([freqs[index_f0], amps.max()-amps.min(), guess_width, amps.min()])
                    temp, temp1 = curve_fit(self.lorentzian, freqs, amps, guess_decay)
                    guess[2] = temp[0]
        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess Detuning-s?'):
                guess[2] = self.getValue('Detuning Guess-s')
            else:
                if self.temp.any():
                    guess[2] = self.temp[0]
                else:
                    response = np.abs(fftshift(fft(amplitudes) / (len(amplitudes)/2.0)))
                    xf = fftshift(fftfreq(len(amplitudes), dT))
                    amps = 1.0*response[len(amplitudes)//2:]
                    amps[0] = 0.0
                    freqs = 1.0*xf[len(amplitudes)//2:]
                    dFreqs = np.abs(freqs[1]-freqs[0])
                    index_f0 = amps.argmax()
                    guess_width = np.abs(index_f0 -(np.abs(amps-amps.mean()-0.5*(amps[index_f0]-amps.mean())).argmin()))*dFreqs
                    guess_decay = np.array([freqs[index_f0], amps.max()-amps.min(), guess_width, amps.min()])
                    temp, temp1 = curve_fit(self.lorentzian, freqs, amps, guess_decay,maxfev=10000) #increased maxfev from default 1000 to 10000 for small Rabi Frequencies
                    guess[2] = temp[0]

        ######## Phase offset guess ###############3
        if self.getValue('X-axis Unit') == 'Frequency [Hz]':
            if self.getValue('Guess Phase Offset-Hz?'):
                guess[3] = np.pi/180.0*self.getValue('Phase Offset-Hz')
            else:
                guess[3] = 0.0

        elif self.getValue('X-axis Unit') == 'Time [s]':
            if self.getValue('Guess Phase Offset-s?'):
                guess[3] = np.pi/180.0*self.getValue('Phase Offset-s')
            else:
                guess[3] = 0.0

        bounds = np.array([ [-np.inf, dT, 0.0, -2.0*np.pi, -np.inf], [np.inf, np.inf, 1.0/dT, 2.0*np.pi, np.inf]])
        try:
            popt, pcov = curve_fit(self.exp_with_sinusoid, t_temp, amplitudes, guess, bounds=bounds)
        except RuntimeError:
            self.log("Error with fit for Guess settings {}".format(guess), level=30)
            popt, pcov = np.array([[1.0, 1.0, 1.0, 1.0, 1.0],[[1.0, 1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0, 1.0],[1.0, 1.0, 1.0, 1.0, 1.0]]])
        popt[4] += offset_correction
        popt[4] *= amp_sign[4]
        popt[0] *= amp_sign[0]
        residuals = (blamplitudes - self.exp_with_sinusoid(t_temp, *popt))

        return popt, pcov, residuals
