import time
import numpy as np

import InstrumentDriver
import AlazarTech_SignalDemodulator_Wrapper as AlazarDig


class Driver(InstrumentDriver.InstrumentWorker):
    """This class implements the AlazarTech card driver."""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection."""
        # keep track of sampled records
        self.data = {}
        #self.lTrace = [np.array([]), np.array([])]
        #self.lSignalNames = ['Channel A - Demodulated values']#extra line for heralding
        self._dt = 1.0
        # open connection
        boardId = int(self.comCfg.address)
        timeout = self.dComCfg['Timeout']
        self.dig = AlazarDig.AlazarTechDigitizer(systemId=1,
                boardId=boardId, timeout=timeout)   
        self.dig.testLED()

    def performClose(self, bError=False, options={}):
        """Perform the close instrument connection operation."""
        # try to remove buffers
        self.dig.removeBuffersDMA()
        # remove digitizer object
        del self.dig
   
    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation. This function
        should return the actual value set by the instrument.
        """
        # start with setting current quant value
        quant.setValue(value)
        # don't do anything until all options are set, then set
        # complete config
        if self.isFinalCall(options):
            self.setConfiguration(options)
        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        # only implemented for records
        if quant.name in ('Channel A - Records',
                          'Channel B - Records',
                          'Channel A - Average record',
                          'Channel B - Average record',
                          'Channel A - Demodulated values',
                          'Channel B - Demodulated values',
                          'Channel A - SNR',
                          'Channel B - SNR',
                          'Channel A - Average demodulated value',
                          'Channel B - Average demodulated value',
                          'Channel A,B - IQ Average demodulated value',
                          'Channel A - Average piecewise demodulated values',
                          'Channel B - Average piecewise demodulated values',
                          'Channel A,B - IQ Average piecewise demodulated values',
                          'Channel A - Average buffer',
                          'Channel B - Average buffer',
                          'Channel A - Average buffer demodulated values',
                          'Channel B - Average buffer demodulated values',
                          'Channel A,B - IQ Average buffer demodulated values',):

            # special case for hardware looping
            if self.isHardwareLoop(options):
                self.getSignalHardwareLoop(quant, options)
                # cash demodulation parameters
                if self._mode.startswith('Individual') or \
                        self._mode.startswith('Referenced Individual'):
                    self.cashDemodulationParametersForIndividual()
                self.log('called performGetValue and isHardwareLoop')
                return self.getData(quant, options)
            
            # check if first call, if so get new records
            if self.isFirstCall(options):
                # clear trace dictionary
                self.data = {}
                # read records
                if self.getValue('NPT AsyncDMA Enabled'):
                    self.getRecordsDMA(
                        hardware_trig=self.isHardwareTrig(options))
                else:
                    self.getRecordsSinglePort()

                # cash demodulation parameters
                if self._mode.startswith('Individual') or \
                        self._mode.startswith('Referenced Individual'):
                    self.cashDemodulationParametersForIndividual()
                    
            
            # return correct data
            return self.getData(quant, options)
        else:
            # just return the quantity value
            return quant.getValue()

    def performArm(self, quant_names, options={}):
        """Perform the instrument arm operation."""
        # configure and start acquisition
        if self.isHardwareLoop(options):
            (seq_no, n_seq) = self.getHardwareLoopIndex(options)
            self.sendValueToOther('Records per buffer', n_seq)
            nRecords = int(round(self.getValue('Number of records')))
            self._nReadout = self.getValue('Number of readout pulse in one waveform')
            self._nRecords = n_seq * nRecords
            # disable trig timeout (set to 1 minute)
            self.dig.AlazarSetTriggerTimeOut(self.dComCfg['Timeout'] + 5*60.0)
            # need to re-configure the card since record size was not
            # known at config
            self.dig.readRecordsDMA(self._mode, self._nSamples,
                    self._nRecords, n_seq,
                    bConfig=True, bArm=True, bMeasure=False,
                    maxBuffers=self._nMaxBuffers,
                    maxBufferSize=self._maxBufferSize)
        else:
            # if not hardware looping, just trigger the card, buffers
            # are already configured
            self.dig.readRecordsDMA(self._mode, self._nSamples,
                    self._nRecords, self._nRecordsPerBuffer,
                    bConfig=False, bArm=True, bMeasure=False,
                    maxBuffers=self._nMaxBuffers,
                    maxBufferSize=self._maxBufferSize)

    def _callbackProgress(self, progress):
        """Report progress to server, as text string."""
        s = 'Acquiring records (%.0f%%).' % (100 * progress)
        self.reportStatus(s)

    def getSignalHardwareLoop(self, quant, options):
        """Get data from round-robin type averaging."""
        (seq_no, n_seq) = self.getHardwareLoopIndex(options)
        # if first sequence call, get data
        if seq_no == 0 and self.isFirstCall(options):
            # show status before starting acquisition
            self.reportStatus('AlazarTech Digitizer - Starting...')
            # get data
            self.data = self.dig.readRecordsDMA(self._mode,
                    self._nSamples, self._nRecords, self._nRecordsPerBuffer,
                    bConfig=False, bArm=False, bMeasure=True,
                    funcStop=self.isStopped,
                    funcProgress=self._callbackProgress,
                    firstTimeout=self.dComCfg['Timeout'] + 5*60.,
                    maxBuffers=self._nMaxBuffers,
                    maxBufferSize=self._maxBufferSize)

    def setConfiguration(self, options={}):
        """Set digitizer configuration based on driver settings."""
        sampleRateIndex = self.getValueIndex('Sample rate')
        # clock configuration
        SourceId = int(self.getCmdStringFromValue('Clock source'))
        if self.getValue('Clock source') == 'Internal':
            # internal
            SampleRateId = int(self.getCmdStringFromValue('Sample rate'),
                    16)
            lFreq = [1E3, 2E3, 5E3, 10E3, 20E3, 50E3, 100E3, 200E3,
                     500E3, 1E6, 2E6, 5E6, 10E6, 20E6, 50E6, 100E6,
                     200E6, 500E6, 1E9, 1.2E9, 1.5E9, 2E9, 2.4E9, 3E9,
                     3.6E9, 4E9]
            Decimation = 0
        elif self.getValue('Clock source') == '10 MHz Reference' and \
                self.getModel() in ('9373', '9360'):
            # 10 MHz ref, for 9373 - decimation is 1
            # for now don't allow DES mode; talk to Simon about best
            # implementation
            lFreq = [1E3, 2E3, 5E3, 10E3, 20E3, 50E3, 100E3, 200E3, 500E3,
                     1E6, 2E6, 5E6, 10E6, 20E6, 50E6, 100E6, 200E6, 500E6, 1E9,
                     1.2E9, 1.5E9, 2E9, 2.4E9, 3E9, 3.6E9, 4E9]
            SampleRateId = int(lFreq[self.getValueIndex('Sample rate')])
            Decimation = 0

        elif self.getValue('Clock source') == '10 MHz Reference' and \
                self.getModel() in ('9870'):
            # 10 MHz ref, use 1GHz rate + divider; divider must be 1, 2,
            # 4, or multiple of 10
            SampleRateId = int(1E9)
            lFreq = [1E3, 2E3, 5E3, 10E3, 20E3, 50E3, 100E3, 200E3,
                     500E3, 1E6, 2E6, 5E6, 10E6, 20E6, 50E6, 100E6,
                     250E6, 500E6, 1E9]
            Decimation = int(1E9 / lFreq[sampleRateIndex])
        elif self.getValue('Clock source') in ('Slow External',
                                               'Fast External') and \
                self.getModel() in ('9870', '9360'):
            SampleRateId = 0x00000040
            lFreq = [1E3, 2E3, 5E3, 10E3, 20E3, 50E3, 100E3, 200E3,
                     500E3, 1E6, 2E6, 5E6, 10E6, 20E6, 50E6, 100E6,
                     250E6, 500E6, 1E9]
            Decimation = 0
        try:
            self.dig.AlazarSetCaptureClock(SourceId, SampleRateId, 0,
                                           Decimation)
        except:
            self.log("'ApiPllNotLocked' has been thrown.")
            time.sleep(.1)
            self.dig.AlazarSetCaptureClock(SourceId, SampleRateId, 0,
                                           Decimation)
        # define time step from sample rate
        self._dt = 1.0 / lFreq[sampleRateIndex]
        # configure inputs
        chnls = {1: 'Channel A', 2: 'Channel B'}
        for n in (1, 2):
            # coupling and range
            if self.getModel() in ('9373', '9360'):
                # these options are not available for these models,
                # set to default
                Coupling = 2
                InputRange = 7
                Impedance = 2
            else:
                Coupling = int(self.getCmdStringFromValue('%s - '
                        'Coupling' % chnls[n]))
                InputRange = int(self.getCmdStringFromValue('%s - '
                        'Range' % chnls[n]))
                Impedance = int(self.getCmdStringFromValue('%s - '
                        'Impedance' % chnls[n]))
            # set coupling, input range, impedance
            self.dig.AlazarInputControl(n, Coupling, InputRange,
                                        Impedance)
            # bandwidth limit, only for model 9870
            if self.getModel() in ('9870',):
                self.dig.AlazarSetBWLimit(n, 0)

        # configure trigger
        Source = int(self.getCmdStringFromValue('Trigger source'))
        Slope = int(self.getCmdStringFromValue('Trigger slope'))
        Delay = self.getValue('Trigger delay')
        timeout = self.dComCfg['Timeout']
        # trigger level is relative to full range
        trigLevel = self.getValue('Trigger level')
        vAmp = np.array([4, 2, 1, 0.4, 0.2, 0.1, .04], dtype=float)
        if self.getValue('Trigger source') == 'Channel A':
            maxLevel = vAmp[self.getValueIndex('Channel A - Range')]
        elif self.getValue('Trigger source') == 'Channel B':
            maxLevel = vAmp[self.getValueIndex('Channel B - Range')]
        elif self.getValue('Trigger source') == 'External':
            maxLevel = 5.0
        elif self.getValue('Trigger source') == 'Immediate':
            maxLevel = 5.0
            # set timeout to very short with immediate triggering
            timeout = 0.001
        # convert relative level to U8
        if abs(trigLevel) > maxLevel:
            trigLevel = maxLevel * np.sign(trigLevel)
        Level = int(128 + 127 * trigLevel / maxLevel)
        # set configuration
        self.dig.AlazarSetTriggerOperation(0, 0, Source, Slope, Level)

        # configure external input, if in use
        if self.getValue('Trigger source') == 'External':
            Coupling = int(self.getCmdStringFromValue('Trigger coupling'))
            self.dig.AlazarSetExternalTrigger(Coupling)

        # set trigger delay and timeout
        delay = int(round(self.getValue('Trigger delay') / self._dt))
        self.dig.AlazarSetTriggerDelay(delay)
        self.dig.AlazarSetTriggerTimeOut(time=timeout)

        # configure memory buffers, only possible when using DMA read
        self._nSamples = int(round(self.getValue('Number of samples')))
        self._nRecords = int(round(self.getValue('Number of records')))
        self._nRecordsPerBuffer = int(round(self.getValue('Records per buffer')))
        self._maxBufferSize = int(round(self.getValue('Max buffer size')))
        self._nMaxBuffers = int(round(self.getValue('Max number of buffers')))
        self._mode = self.getValue('Acquisition mode')
        
        if self._mode != 'Raw':
            self._length = int(round(self.getValue('Demodulation length') / self._dt))
            self._dFreq = self.getValue('Demodulation frequency')
            self._skip = int(round(self.getValue('Skip start') / self._dt))
            if self._length > self._nSamples:
                self._length = self._nSamples
            if self._length < 1:
                raise Exception('Nothing to demodulate: increase '
                    'the demodulation length or record length.')
                   
        if self._mode.startswith('Referenced'):
            self._bRef = True
        else:
            self._bRef = False

        self.perform_IQdemodulation = self.getValue("Use IQ demodulation")

        if not self.isHardwareLoop(options):
            # configure DMA read
            self.dig.readRecordsDMA(self._mode, self._nSamples,
                    self._nRecords, self._nRecordsPerBuffer,
                    bConfig=True, bArm=False, bMeasure=False,
                    maxBuffers=self._nMaxBuffers,
                    maxBufferSize=self._maxBufferSize)

    def getRecordsDMA(self, hardware_trig=False):
        """Resample the data with DMA."""
        # in hardware triggering mode, there is no need to re-arm
        # the card
        bArm = not hardware_trig
        # get data
        self.data = self.dig.readRecordsDMA(self._mode, self._nSamples,
                self._nRecords, self._nRecordsPerBuffer,
                bConfig=False, bArm=bArm, bMeasure=True,
                funcStop=self.isStopped,
                funcProgress=self._callbackProgress,
                firstTimeout=self.dComCfg['Timeout'] + 300.,
                maxBuffers=self._nMaxBuffers,
                maxBufferSize=self._maxBufferSize)

    def getRecordsSinglePort(self):
        """Resample the data."""
        # get channels in use
        nPreSize = int(round(self.getValue('Pre-trigger samples')))
        nPostSize = self._nSamples

        model = self.getModel()
        if model == '9870' or model == '9360':
            pAlgn = 64
            rAlgn = 256
            nPreAlgn = pAlgn * ((nPreSize + pAlgn - 1) // pAlgn)
            nPostAlgn = pAlgn * ((nPostSize + pAlgn - 1) // pAlgn)
            nTotal = rAlgn * ((nPreAlgn + nPostAlgn + rAlgn - 1) // rAlgn)
            nPostAlgn = nTotal - nPreAlgn
            start = nPreAlgn - nPreSize
            end = start + nPreSize + nPostSize
        else:
            raise NotImplementedError('Model ATS%s is not yet supported.'
                                      % model)

        self.dig.AlazarSetRecordSize(nPreAlgn, nPostAlgn)
        self.dig.AlazarSetRecordCount(self._nRecords)
        # start aquisition
        self.dig.AlazarStartCapture()
        nTry = self.dComCfg['Timeout'] / 0.05
        while nTry and self.dig.AlazarBusy() and not self.isStopped():
            # sleep for a while to save resources, then try again
            self.wait(0.050)
            nTry -= 1
        # check if timeout occurred
        if nTry <= 0:
            self.dig.AlazarAbortCapture()
            raise TimeoutError('Acquisition timed out.')
        # check if user stopped
        if self.isStopped():
            self.dig.AlazarAbortCapture()
            return

        # read data for channels in use
        mode = self._mode

        records = self.dig.readRecordsSinglePort(1)
        self.data['Channel A'] = records[:,start:end]
        if mode.startswith('Average') or \
                mode.startswith('Referenced Average'):
            self.data['Channel A - Average record'] = \
                    np.mean(self.data['Channel A'])
        
        records = self.dig.readRecordsSinglePort(2)
        self.data['Channel B'] = records[:,start:end]
        if mode.startswith('Average') or \
                mode.startswith('Referenced Average'):
            self.data['Channel B - Average record'] = \
                    np.mean(self.data['Channel B'])

    def _raiseError(self, name, mode):
            raise NotImplementedError("Output data '%s' could not be "
                    "acquired in acquisition mode '%s'." % (name, mode))

    def getData(self, quant, options={}):
        """
        Return data that corresponds to the selected data acquisition
        mode.
        """
        name = quant.name
        mode = self._mode
        if mode == 'Raw':
            if name in ('Channel A - Records',
                        'Channel B - Records'):
                flattened = self.data[name[:9]].ravel()
                return quant.getTraceDict(flattened, dt=self._dt)
            else:
                raise self._raiseError(name, mode)
                
        elif mode == 'Individual Record Demodulation':
            if name in ('Channel A - Records',
                        'Channel B - Records'):
                flattened = self.data[name[:9]].ravel()
                return quant.getTraceDict(flattened, dt=self._dt)
            elif name in ('Channel A - Demodulated values',
                          'Channel B - Demodulated values',
                          'Channel A,B - IQ Demodulated values'):
                if name not in self.data:
                    self.getDemodulatedValuesFromIndividual()  
                return quant.getTraceDict(self.data[name], dt=1)
               
            elif name in ('Channel A - Average demodulated value',
                          'Channel B - Average demodulated value',
                          'Channel A,B - IQ Average demodulated value'):
                if name not in self.data:
                    t0 = time.clock()
                    self.getDemodulatedValuesFromIndividual()
                    self.log('Data processing %.6f s.' % (time.clock() - t0))
                    if hasattr(self, '_tcycle'):
                        self.log('Total time %.6f s.' % (time.clock() - self._tcycle))
                    self._tcycle = time.clock()
                return self.data[name]
            elif name in ('Channel A - Average record',
                          'Channel B - Average record'):
                if name not in self.data:
                    self.getAverageRecordFromIndiviual()
                return quant.getTraceDict(self.data[name], dt=self._dt)
            elif name in ('Channel A - Average piecewise demodulated values',
                          'Channel B - Average piecewise demodulated values',
                          'Channel A,B - IQ Average piecewise demodulated values'):
                if name not in self.data:
                    if ('%s - Average record' % name[:9]) not in self.data:
                        self.getAverageRecordFromIndiviual() 
                    self.getPiecewiseDemodulatedValuesFromAverage()
                dt = self.getValue('Demodulation length')
                return quant.getTraceDict(self.data[name], dt=dt)
            elif name in ('Channel A - SNR', 'Channel B - SNR'):
                if name not in self.data:
                        self.getDemodulatedValuesFromIndividual()
                return self.data[name]
            else:
                raise self._raiseError(name, mode)
                
        elif mode == 'Referenced Individual Record Demodulation':
            if name in ('Channel A - Records',
                        'Channel B - Records'):
                flattened = self.data[name[:9]].ravel()
                return quant.getTraceDict(flattened, dt=self._dt)
            elif name == 'Channel A - Demodulated values':
                if self.isHardwareLoop(options):
                    (seq_no, n_seq) = self.getHardwareLoopIndex(options)
                    nReadout = int(self._nReadout)
                    if seq_no == 0:
                        self.getDemodulatedValuesFromIndividual()
                        self.data_reshape = self.data[name].reshape((int(self._nRecords / n_seq / nReadout), n_seq * nReadout))
                    self.log('get data!')
                    self.log(self.data_reshape.shape)
                    data_for_seq_no = self.data_reshape[:, seq_no:seq_no + nReadout].ravel()
                    self.log(data_for_seq_no.shape)
                    return data_for_seq_no
                
                if name not in self.data:
                    self.getDemodulatedValuesFromIndividual()
                return quant.getTraceDict(self.data[name], dt=1)
            elif name == 'Channel A - Average demodulated value':
                if name not in self.data:
                    t0 = time.clock()
                    self.getDemodulatedValuesFromIndividual()
                    self.log('Data processing %.6f s.' % (time.clock() - t0))
                    if hasattr(self, '_tcycle'):
                        self.log('Total time %.6f s.' % (time.clock() - self._tcycle))
                    self._tcycle = time.clock()
                return self.data[name]
            elif name == 'Channel A - Average record':
                if name not in self.data:
                    self.getAverageRecordFromIndiviual()
                return quant.getTraceDict(self.data[name], dt=self._dt)
            elif name == 'Channel A - Average piecewise demodulated values':
                t0 = time.clock()
                if name not in self.data:
                    if 'Channel A - Average record' not in self.data:
                        self.getAverageRecordFromIndiviual() 
                    self.getPiecewiseDemodulatedValuesFromAverage()
                self.log('Data processing %.6f s.' % (time.clock() - t0))
                if hasattr(self, '_tcycle'):
                    self.log('Total time %.6f s.' % (time.clock() - self._tcycle))
                self._tcycle = time.clock()
                dt = self.getValue('Demodulation length')
                return quant.getTraceDict(self.data[name], dt=dt)
            elif name == 'Channel A - SNR':
                if name not in self.data:
                    self.getDemodulatedValuesFromIndividual()
                return self.data[name]
            else:
                raise self._raiseError(name, mode)
                
        elif mode == 'Average Record Demodulation':
            if name in ('Channel A - Average record',
                        'Channel B - Average record'):
                return quant.getTraceDict(self.data[name], dt=self._dt)
            elif name in ('Channel A - Average demodulated value',
                          'Channel B - Average demodulated value',
                          'Channel A,B - IQ Average demodulated value'):
                self.getDemodulatedValueFromAverage()
                return self.data[name]
            elif name in ('Channel A - Average piecewise demodulated values',
                          'Channel B - Average piecewise demodulated values',
                          'Channel A,B - IQ Average piecewise demodulated values'):
                if name not in self.data:
                    self.getPiecewiseDemodulatedValuesFromAverage()
                dt = self.getValue('Demodulation length')
                return quant.getTraceDict(self.data[name], dt=dt)
            else:
                raise self._raiseError(name, mode)
                              
        elif mode == 'Referenced Average Record Demodulation':
            if name in ('Channel A - Average record',
                        'Channel B - Average record'):
                return quant.getTraceDict(self.data[name], dt=self._dt)
            elif name == 'Channel A - Average demodulated value':
                t0 = time.clock()
                self.getDemodulatedValueFromAverage()
                self.log('Data processing %.6f s.' % (time.clock() - t0))
                if hasattr(self, '_tcycle'):
                    self.log('Total time %.6f s.' % (time.clock() - self._tcycle))
                self._tcycle = time.clock()
                return self.data[name]
            elif name == 'Channel A - Average piecewise demodulated values':
                t0 = time.clock()
                if name not in self.data:
                    self.getPiecewiseDemodulatedValuesFromAverage()
                dt = self.getValue('Demodulation length')
                self.log('Data processing %.6f s.' % (time.clock() - t0))
                if hasattr(self, '_tcycle'):
                    self.log('Total time %.6f s.' % (time.clock() - self._tcycle))
                self._tcycle = time.clock()
                return quant.getTraceDict(self.data[name], dt=dt)
            else:
                raise self._raiseError(name, mode)

        elif mode == 'Average Buffer Demodulation':
            if name in ('Channel A - Average buffer',
                        'Channel B - Average buffer'):
                return quant.getTraceDict(self.data[name].flatten(),
                                          dt=self._dt)
            elif name in ('Channel A - Average buffer demodulated values',
                          'Channel B - Average buffer demodulated values',
                          'Channel A,B - IQ Average buffer demodulated values'):
                if name not in self.data:
                    self.getDemodulatedValuesFromBuffer()
                dt = self.getValue('Sequence time step')
                return quant.getTraceDict(self.data[name], dt=dt)            
            elif self.isHardwareLoop(options) and name in \
                    ('Channel A - Average demodulated value',
                     'Channel B - Average demodulated value',
                     'Channel A,B - IQ Average demodulated value'):
                (seq_no, n_seq) = self.getHardwareLoopIndex(options)
                if seq_no == 0:
                    self.getDemodulatedValuesFromBuffer()
                if name.startswith('Channel A -'):
                    return self.data['Channel A - Average buffer demodulated values'][seq_no]
                elif name.startswith('Channel B -'):
                    return self.data['Channel B - Average buffer demodulated values'][seq_no]
                elif name.startswith('Channel A,B -'):
                    return self.data['Channel A,B - IQ Average buffer demodulated values'][seq_no]
            elif self.isHardwareLoop(options) and \
                    name in ('Channel A - Average record',
                             'Channel B - Average record'):
                (seq_no, n_seq) = self.getHardwareLoopIndex(options)
                buff_name = name.replace('record', 'buffer')
                return quant.getTraceDict(self.data[buff_name][seq_no],
                        dt=self._dt)
            else:
                raise self._raiseError(name, mode)
                              
        elif mode == 'Referenced Average Buffer Demodulation':
            if name in ('Channel A - Average buffer',
                        'Channel B - Average buffer'):
                return quant.getTraceDict(self.data[name].flatten(),
                                          dt=self._dt)
            elif name == 'Channel A - Average buffer demodulated values':
                if name not in self.data:
                    self.getDemodulatedValuesFromBuffer()
                dt = self.getValue('Sequence time step')
                return quant.getTraceDict(self.data[name], dt=dt)
            elif self.isHardwareLoop(options) and \
                    name == 'Channel A - Average demodulated value':
                (seq_no, n_seq) = self.getHardwareLoopIndex(options)
                if seq_no == 0:
                    self.getDemodulatedValuesFromBuffer()
                return self.data['Channel A - Average buffer demodulated values'][seq_no]
            elif self.isHardwareLoop(options) and \
                    name in ('Channel A - Average record',
                             'Channel B - Average record'):
                (seq_no, n_seq) = self.getHardwareLoopIndex(options)
                buff_name = name.replace('record', 'buffer')
                return quant.getTraceDict(self.data[buff_name][seq_no],
                        dt=self._dt)
            else:
                raise self._raiseError(name, mode)
        else:
            raise self._raiseError(name, mode)

    def cashDemodulationParametersForIndividual(self):
        """Cash parameters that are used for calculating output values."""
        # get parameters
        nTotLength = self.data['Channel A'].size
        dt = self._dt
        dFreq = self._dFreq
        skip = self._skip
        length = self._length
        nRecords = self._nRecords
        recordLength = nTotLength // nRecords
        length = min(length, recordLength - skip)

        # calculate cos/sin vectors, allow segmenting
        if not hasattr(self, '_prev_dt') or \
                dt != self._prev_dt or \
                not hasattr(self, '_prev_skip') or \
                skip != self._prev_skip or \
                not hasattr(self, '_prev_length') or \
                length != self._prev_length or \
                not hasattr(self, '_prev_dFreq') or \
                dFreq != self._prev_dFreq:
            vTime = dt * (skip + np.arange(length, dtype=np.float32))
            self._vExp = np.exp(2.j * np.pi * vTime * dFreq).view('complex64')
            self._prev_dt = dt
            self._prev_skip = skip
            self._prev_length = length
            self._prev_dFreq = dFreq
        
        if self._bRef:
            vChB = self.data['Channel B'][:,skip:skip+length]
            vDemodRefs = np.dot(vChB, self._vExp)
            self._vExpRef = np.exp(-1.j * np.angle(vDemodRefs))

    def getAverageRecordFromIndiviual(self):
        """Calculate time average from data and reference."""
        if self._bRef:
            timeAverage = np.dot(self.data['Channel A'].T, self._vExpRef)
            timeAverage /= self._nRecords
            self.data['Channel A - Average record'] = timeAverage
        else:
            for ch in ('Channel A', 'Channel B'):
                timeAverage = np.sum(self.data[ch], axis=0)
                timeAverage /= self._nRecords
                self.data['%s - Average record' % ch] = timeAverage

    def getDemodulatedValuesFromIndividual(self):
        """Calculate complex signal vector from data and reference."""
        skip = self._skip
        length = self._length

        if self._bRef:
            vChA = self.data['Channel A'][:,skip:skip+length]
            vDemodVals = np.dot(vChA, self._vExp)
            vDemodVals /= np.float32(.5) * np.float32(length)
            vDemodVals *= self._vExpRef
            self.data['Channel A - Demodulated values'] = vDemodVals
            meanDemodVal = np.mean(vDemodVals)
            self.data['Channel A - Average demodulated value'] = meanDemodVal
            self.data['Channel A - SNR'] = np.abs(meanDemodVal) / np.std(vDemodVals)
        else:
            for ch in ('Channel A', 'Channel B'):
                vCh = self.data[ch][:,skip:skip+length]
                vDemodVals = np.dot(vCh, self._vExp)
                vDemodVals /= np.float32(.5) * np.float32(length)
                self.data['%s - Demodulated values' % ch] = vDemodVals
                meanDemodVal = np.mean(vDemodVals)
                self.data['%s - Average demodulated value' % ch] = meanDemodVal
                self.data['%s - SNR' % ch] = np.abs(meanDemodVal) / np.std(vDemodVals)
        if self.perform_IQdemodulation:
            self.data['Channel A,B - IQ Demodulated value'] = \
                self.data['Channel A - Demodulated value'] - 1.j * self.data['Channel B - Demodulated value']
            self.data['Channel A,B - IQ Average demodulated value'] = \
                self.data['Channel A - Average demodulated value'] - 1.j * self.data['Channel B - Average demodulated value']

    def getDemodulatedValueFromAverage(self):
        """Calculate complex signal vector from data and reference."""
        # get parameters
        dFreq = self._dFreq
        skip = self._skip
        length = self._length
        vTime = self._dt * (skip + np.arange(length, dtype=np.float32))
        vExp = np.exp(2.j * np.pi * vTime * dFreq).view('complex64')

        if self._bRef:
            vChA = self.data['Channel A - Average record'][skip:skip+length]
            vDemodVal = np.dot(vChA, vExp)
            vDemodVal /= np.float32(.5) * np.float32(length)
            vChB = self.data['Channel B - Average record'][skip:skip+length]
            vDemodRef = np.dot(vChB, vExp)
            vDemodVal *= np.exp(-1.j * np.angle(vDemodRef))
            self.data['Channel A - Average demodulated value'] = vDemodVal
        else:
            for ch in ('Channel A', 'Channel B'):
                vCh = self.data['%s - Average record' % ch][skip:skip+length]
                vDemodVal = np.dot(vCh, vExp)
                vDemodVal /= .5 * np.float32(length -1)
                self.data['%s - Average demodulated value' % ch] = vDemodVal
        if self.perform_IQdemodulation:
            self.data['Channel A,B - IQ Average demodulated value'] = \
                self.data['Channel A - Average demodulated value'] - 1.j * self.data['Channel B - Average demodulated value']

    def getPiecewiseDemodulatedValuesFromAverage(self):
        # get parameters
        dFreq = self._dFreq
        skip = self._skip
        length = self._length
        nSamples = self._nSamples
        vTime = self._dt * (skip + np.arange(nSamples - skip, dtype=np.float32))
        vExp = np.exp(2.j * np.pi * vTime * dFreq).view('complex64')
        nSegments = (nSamples - skip) // length
        total_length = nSegments * length

        if total_length != nSamples:
            vExp = vExp[:total_length]
        vExp.shape = (nSegments, length)

        data = {}
        for ch in ('Channel A', 'Channel B'):
            vDemodVals = np.empty(nSegments, dtype=np.complex64)
            vCh = self.data['%s - Average record' % ch]
            if skip > 0 or total_length != nSamples:
                vCh = vCh[skip:skip+total_length]
            vCh.shape = (nSegments, length)
            for idx in range(nSegments):
                vDemodVals[idx] = np.dot(vCh[idx], vExp[idx])
            vDemodVals /= np.float32(.5) * np.float32(length)
            data[ch] = vDemodVals
  
        if self._bRef:
            data['Channel A'] *= np.exp(-1.j * np.angle(data['Channel B']))
        else:         
            self.data['Channel B - Average piecewise '
                      'demodulated values'] = data['Channel B']
        self.data['Channel A - Average piecewise '
                  'demodulated values'] = data['Channel A']

        if self.perform_IQdemodulation:
            self.data['Channel A,B - IQ Average piecewise '
                  'demodulated values'] = self.data['Channel A - Average piecewise '
                  'demodulated values'] - 1.j * self.data['Channel B - Average piecewise '
                  'demodulated values']

    def getDemodulatedValuesFromBuffer(self):
        # get parameters
        dFreq = self._dFreq
        skip = self._skip
        length = self._length
        vTime = self._dt * (skip + np.arange(length, dtype=np.float32))
        vExp = np.exp(2.j * np.pi * vTime * dFreq).view('complex64')

        if self._bRef:
            vChA = self.data['Channel A - Average buffer'][:,skip:skip+length]
            vDemodVal = np.dot(vChA, vExp)
            vDemodVal /= np.float32(.5) * np.float32(length)
            vChB = self.data['Channel B - Average buffer'][:,skip:skip+length]
            vDemodRef = np.dot(vChB, vExp)
            vDemodVal *= np.exp(-1.j * np.angle(vDemodRef))
            self.data['Channel A - Average buffer demodulated values'] = vDemodVal
        else:
            for ch in ('Channel A', 'Channel B'):
                vCh = self.data['%s - Average buffer' % ch][:,skip:skip+length]
                vDemodVal = np.dot(vCh, vExp)
                vDemodVal /= .5 * np.float32(length -1)
                self.data['%s - Average buffer demodulated values' % ch] = vDemodVal
        if self.perform_IQdemodulation:
            self.data['Channel A,B - IQ Average buffer demodulated values'] = \
                self.data['Channel A - Average buffer demodulated values'] - 1.j * self.data['Channel B - Average buffer demodulated values']

if __name__ == '__main__':
    pass
