import time
import numpy as np
import ctypes
from ctypes import (c_int, c_uint8, c_uint16, 
                    c_char_p, c_void_p, c_long, byref)

# add logger, to allow logging to Labber's instrument log
import logging
log = logging.getLogger('LabberDriver')

from atsapi import *


class AlazarTechDigitizer():
    """Represents the AlazarTech digitizer, redefines the DLL functions
    in Python.
    """

    def __init__(self, systemId=1, boardId=1, timeout=10.0):
        """Defines a session ID, used to identify the instrument."""
        # range settings; default value of 400 mV for model 9373;
        # will be overwritten if model is 9870 and AlazarInputControl is
        # called
        self.dRange = {1: 0.4, 2: 0.4}
        self.buffers = []
        self.timeout = timeout
        # create a session ID
        func = getattr(ats, 'AlazarNumOfSystems')
        func.restype = U32
        func = getattr(ats, 'AlazarGetBoardBySystemID')
        func.restype = c_void_p
        handle = func(U32(systemId), U32(boardId))
        if handle is None:
            raise Exception('Device with system ID=%d and board ID=%d '
                            'could not be found.' % (systemId, boardId))
        self.handle = handle
        # get memory and bitsize
        (self.memorySize, bitsPerSample) = self.AlazarGetChannelInfo()
        self.bytesPerSample = (bitsPerSample + 7) // 8
        # range and zero for conversion to voltages
        self.codeZero = np.float32(1 << (bitsPerSample - 1))
        self.codeZero -= np.float32(.5)
        self.codeRange = np.float32(1 << (bitsPerSample - 1))
        self.codeRange -= np.float32(.5)
        self.nChannels = 2

    def testLED(self):
        self.callFunc('AlazarSetLED', self.handle, U32(1))
        time.sleep(0.1)
        self.callFunc('AlazarSetLED', self.handle, U32(0))

    def callFunc(self, sFunc, *args, **kargs):
        """General function caller with restype=status,
        also checks for errors.
        """
        # get function from DLL
        func = getattr(ats, sFunc)
        func.restype = c_int
        # call function, raise error if needed
        status = func(*args)
        if 'bIgnoreError' in kargs:
            bIgnoreError = kargs['bIgnoreError']
        else:
            bIgnoreError = False
        if status is not None and status > 512 and not bIgnoreError:
            raise Exception(self.getError(status))

    def getError(self, status):
        """Convert the error in status to a string."""
        func = getattr(ats, 'AlazarErrorToText')
        func.restype = c_char_p
        # const char* AlazarErrorToText(RETURN_CODE retCode)
        errorText = func(c_int(status))
        return str(errorText)

    def AlazarGetChannelInfo(self):
        """Get the on-board memory in samples per channel and sample
        size in bits per sample.
        """
        memorySize_samples = U32(0)
        bitsPerSample = U8(0)
        self.callFunc('AlazarGetChannelInfo', self.handle,
                      byref(memorySize_samples), byref(bitsPerSample))
        return (int(memorySize_samples.value), int(bitsPerSample.value))

    # RETURN_CODE AlazarSetCaptureClock(HANDLE h, U32 Source, U32 Rate,
    #       U32 Edge, U32 Decimation);
    def AlazarSetCaptureClock(self, SourceId, SampleRateId, EdgeId=0,
                              Decimation=0):
        self.callFunc('AlazarSetCaptureClock', self.handle, U32(SourceId),
                      U32(SampleRateId), U32(EdgeId), U32(Decimation))

    # RETURN_CODE AlazarInputControl(HANDLE h, U8 Channel,
    #       U32 Coupling, U32 InputRange, U32 Impedance);
    def AlazarInputControl(self, Channel, Coupling, InputRange, Impedance):
        # keep track of input range
        dConv = {12: 4.0, 11: 2.0, 10: 1.0, 7: 0.4,
                 6: 0.2, 5: 0.1, 2: 0.04}
        self.dRange[Channel] = dConv[InputRange]
        self.callFunc('AlazarInputControl', self.handle, U8(Channel),
                      U32(Coupling), U32(InputRange), U32(Impedance))

    # RETURN_CODE AlazarSetBWLimit(HANDLE h, U8 Channel, U32 enable);
    def AlazarSetBWLimit(self, Channel, enable):
        self.callFunc('AlazarSetBWLimit', self.handle, U32(Channel),
                      U32(enable))

    # RETURN_CODE AlazarSetTriggerOperation(HANDLE h, U32 TriggerOperation,
    #       U32 TriggerEngine1/*J,K*/, U32 Source1, U32 Slope1, U32 Level1,
    #       U32 TriggerEngine2/*J,K*/, U32 Source2, U32 Slope2, U32 Level2);
    def AlazarSetTriggerOperation(self, TriggerOperation=0,
            TriggerEngine1=0, Source1=0, Slope1=1, Level1=128,
            TriggerEngine2=1, Source2=3, Slope2=1, Level2=128):
        self.callFunc('AlazarSetTriggerOperation', self.handle,
            U32(TriggerOperation),
            U32(TriggerEngine1), U32(Source1), U32(Slope1), U32(Level1),
            U32(TriggerEngine2), U32(Source2), U32(Slope2), U32(Level2))

    # RETURN_CODE AlazarSetExternalTrigger(HANDLE h, U32 Coupling, U32 Range);
    def AlazarSetExternalTrigger(self, Coupling, Range=3):# Range=0): MODIFICATION
        self.callFunc('AlazarSetExternalTrigger', self.handle,
                      U32(Coupling), U32(Range))

    # RETURN_CODE AlazarSetTriggerDelay( HANDLE h, U32 Delay);
    def AlazarSetTriggerDelay(self, Delay=0):
        self.callFunc('AlazarSetTriggerDelay', self.handle, U32(Delay))

    # RETURN_CODE AlazarSetTriggerTimeOut( HANDLE h, U32 to_ns);
    def AlazarSetTriggerTimeOut(self, time=0.0):
        tick = U32(int(time * 1E5))
        self.callFunc('AlazarSetTriggerTimeOut', self.handle, tick)

    # RETURN_CODE AlazarSetRecordSize(HANDLE h, U32 PreSize, U32 PostSize);
    def AlazarSetRecordSize(self, PreSize, PostSize):
        self.nPreSize = int(PreSize)
        self.nPostSize = int(PostSize)
        self.callFunc('AlazarSetRecordSize', self.handle, U32(PreSize),
                      U32(PostSize))

    # RETURN_CODE AlazarSetRecordCount(HANDLE h, U32 Count);
    def AlazarSetRecordCount(self, Count):
        self.nRecord = int(Count)
        self.callFunc('AlazarSetRecordCount', self.handle, U32(Count))

    # RETURN_CODE AlazarStartCapture(HANDLE h);
    def AlazarStartCapture(self):
        self.callFunc('AlazarStartCapture', self.handle)

    # RETURN_CODE AlazarAbortCapture(HANDLE h);
    def AlazarAbortCapture(self):
        self.callFunc('AlazarAbortCapture', self.handle)

    # U32 AlazarBusy(HANDLE h);
    def AlazarBusy(self):
        # get function from DLL
        func = getattr(ats, 'AlazarBusy')
        func.restype = U32
        # call function, return result
        return bool(func(self.handle))

    # U32 AlazarRead(HANDLE h, U32 Channel, void *buffer, int ElementSize,
    #       long Record, long TransferOffset, U32 TransferLength);
    def AlazarRead(self, Channel, buffer, ElementSize,
                   Record, TransferOffset, TransferLength):
        self.callFunc('AlazarRead', self.handle, U32(Channel), buffer,
            c_int(ElementSize), c_long(Record), c_long(TransferOffset),
            U32(TransferLength))

    def AlazarBeforeAsyncRead(self, channels, transferOffset,
                              samplesPerRecord, nRecordsPerBuffer,
                              recordsPerAcquisition, flags):
        """Prepares the board for an asynchronous acquisition."""
        self.callFunc('AlazarBeforeAsyncRead', self.handle,
                      channels, transferOffset, samplesPerRecord,
                      nRecordsPerBuffer, recordsPerAcquisition, flags)

    # RETURN_CODE AlazarAbortAsyncRead( HANDLE h);
    def AlazarAbortAsyncRead(self):
        """Cancels any asynchronous acquisition running on a board."""
        self.callFunc('AlazarAbortAsyncRead', self.handle)

    def AlazarPostAsyncBuffer(self, buffer, bufferLength):
        """Posts a DMA buffer to a board."""
        self.callFunc('AlazarPostAsyncBuffer', self.handle, buffer,
                      bufferLength)

    def AlazarWaitAsyncBufferComplete(self, buffer, timeout_ms):
        """Blocks until the board confirms that buffer is filled with
        data.
        """
        self.callFunc('AlazarWaitAsyncBufferComplete', self.handle,
                      buffer, timeout_ms)

    def readRecordsDMA(self, mode, nSamples, nRecord, nRecordsPerBuffer,
                bConfig=True, bArm=True, bMeasure=True,
                funcStop=None, funcProgress=None,
                timeout=None, firstTimeout=None,
                maxBuffers=8, maxBufferSize=(2**26)):
        """
        Reads records in NPT AutoDMA mode, converts to float,
        demodulates/averages to a single trace.
        """
        data = {}
        # use global timeout if not given
        timeout = self.timeout if timeout is None else timeout
        # first timeout can be different in case of slow initial arming
        firstTimeout = timeout if firstTimeout is None else firstTimeout

        # change alignment to be 64
        samplesPerRecord = 64 * ((nSamples + 63) // 64)

        # compute the number of bytes per record and per buffer
        bytesPerRecord = self.bytesPerSample * samplesPerRecord

        if self.nChannels * bytesPerRecord > maxBufferSize:
            raise MemoryError('Maximum allowed buffer size '
                    'is too small to contain even a single record.')

        expansion = 1
        if nRecord % nRecordsPerBuffer != 0:
            log.info('Recomputing the number of records per buffer.')
            nBuffersPerAcquisition = 1
            nRecordsPerBuffer = nRecord
            nResidual = nRecord
            consistencyFlag = False
            while not consistencyFlag:
                bytesPerBuffer = int(self.nChannels *
                        nRecordsPerBuffer * bytesPerRecord)
                if bytesPerBuffer > maxBufferSize:
                    noFactorFlag = True
                    for factor in range(2, int(np.sqrt(nResidual)) + 2):
                        if nResidual % factor == 0:
                            nRecordsPerBuffer = int(nRecordsPerBuffer / factor)
                            nResidual = int(nResidual / factor)
                            noFactorFlag = False
                            break
                    if noFactorFlag:
                        break
                else:
                    consistencyFlag = True
            if not consistencyFlag:
                raise MemoryError('The number of records and the number'
                    ' of samples in a single record is not consistent '
                    'with the allowed maximum buffer size. Try to use '
                    'parameter values that could be factorized into '
                    'small prime numbers (values that are powers of 2 '
                    'are most preferable). Alternatively, try to '
                    'increase the maximum buffer size.')
        else:
            bytesPerBuffer = int(self.nChannels * nRecordsPerBuffer *
                    bytesPerRecord)
            for factor in range(2, int(maxBufferSize / bytesPerBuffer) + 1):
                if nRecord % (factor * nRecordsPerBuffer) == 0 and \
                        bytesPerBuffer * factor < maxBufferSize:
                    nRecordsPerBuffer *= factor
                    bytesPerBuffer *= factor
                    expansion *= factor
        
        log.info('Number of records: %d' % nRecord)
        log.info('Number of records per buffer: %d' % nRecordsPerBuffer)

        nBuffersPerAcquisition = int(nRecord / nRecordsPerBuffer)
        log.info('Number of buffers per acquisition: %d'
                % nBuffersPerAcquisition)
        # do not allocate more buffers than needed for all data
        bufferCount = int(min(2 * ((nBuffersPerAcquisition + 1) // 2),
                              maxBuffers))                          
        if mode.startswith('Average Record') or \
                mode.startswith('Referenced Average Record'):
            mode = 'a'
        elif mode.startswith('Average Buffer') or \
                mode.startswith('Referenced Average Buffer'):
            mode = 'b'
        else:
            mode = 'i'

        # configure board
        if not hasattr(self, '_samplesPerRecord') or \
                samplesPerRecord != self._samplesPerRecord:
            bConfig = True
            self._samplesPerRecord = samplesPerRecord
        if not hasattr(self, '_nRecord') or nRecord != self._nRecord:
            bConfig = True
            self._nRecord = nRecord
        if not hasattr(self, '_bytesPerBuffer') or \
                bytesPerBuffer != self._bytesPerBuffer:
            bConfig = True
            self._bytesPerBuffer = bytesPerBuffer
        if not hasattr(self, '_bufferCount') or \
                bufferCount != self._bufferCount:
            bConfig = True
            self._bufferCount = bufferCount

        if bConfig:
            log.info('Requested number of buffers: %d' % bufferCount)
            log.info('Buffer size [MB]: %f' % (float(bytesPerBuffer / 2**20)))
            log.info('Records per buffer: %d' % nRecordsPerBuffer)
            self.AlazarSetRecordSize(0, samplesPerRecord)
            self.AlazarSetRecordCount(nRecord)
            # allocate DMA buffers
            if self.bytesPerSample > 1:
                sample_type = ctypes.c_uint16
            else:
                sample_type = ctypes.c_uint8
            # clear old buffers
            self.removeBuffersDMA()
            # create new buffers
            self.buffers = []
            for i in range(bufferCount):
                self.buffers.append(DMABuffer(sample_type,
                                              bytesPerBuffer))
            # initialize data arrays
            if mode == 'i':
                samples = self.nChannels * nRecord * samplesPerRecord
                if self.bytesPerSample == 1:
                    self._records = np.empty(samples, dtype=np.uint8)
                elif self.bytesPerSample == 2:
                    self._records = np.empty(samples, dtype=np.uint16)
                else:
                    self._records = np.empty(samples, dtype=np.uint32)

        # arm and start capture, if wanted
        if bArm:
            # configure the board to make an NPT AutoDMA acquisition
            self.AlazarBeforeAsyncRead(2, 0,
                    samplesPerRecord,
                    nRecordsPerBuffer,
                    nRecord,
                    ADMA_EXTERNAL_STARTCAPTURE | ADMA_NPT)
            # Post DMA buffers to board
            for buf in self.buffers:
                self.AlazarPostAsyncBuffer(buf.addr, buf.size_bytes)
            try:
                self.AlazarStartCapture()
            except:
                # make sure buffers release memory if failed
                self.removeBuffersDMA()
                raise

        # if not waiting for result, return here
        if not bMeasure:
            return

        samplesPerBuffer = (self.nChannels * \
                            nRecordsPerBuffer * \
                            samplesPerRecord)
        if mode == 'i':
            records = self._records
        else:
            if self.bytesPerSample == 1 and nRecord < 2**24:
                avgRecord = np.zeros(samplesPerBuffer, dtype=np.uint32)
            else:
                avgRecord = np.zeros(samplesPerBuffer, dtype=np.uint64)

        # range and zero for each channel, combined with bit shifting
        rangeA = np.float32(self.dRange[1] / self.codeRange)
        rangeB = np.float32(self.dRange[2] / self.codeRange)

        buffersCompleted = 0
        timeout_ms = int(1000. * firstTimeout + 5.0*60)
        t0 = time.clock()
        try:
            while buffersCompleted < nBuffersPerAcquisition:
                # wait for the buffer at the head of the list of
                # available buffers to be filled by the board
                buf = self.buffers[buffersCompleted % bufferCount]
                self.AlazarWaitAsyncBufferComplete(buf.addr,
                                                   timeout_ms=timeout_ms)

                if mode == 'i':
                    bufferPosition = samplesPerBuffer * buffersCompleted
                    records[bufferPosition:bufferPosition+samplesPerBuffer] = \
                        buf.buffer
                else:
                    avgRecord += buf.buffer

                # add the buffer to the end of the list of available
                # buffers
                self.AlazarPostAsyncBuffer(buf.addr, buf.size_bytes)
                
                buffersCompleted += 1
                
                # break if stopped from outside
                if funcStop is not None and funcStop():
                    break
                # reset timeout time, can be different than first call
                timeout_ms = int(1000. * timeout)
                # report progress
                if funcProgress is not None:
                    funcProgress(float(buffersCompleted) /
                                 float(nBuffersPerAcquisition))
        finally:
            # release resources
            self.AlazarAbortAsyncRead()
        
        log.info('Data acquisition: %.6f s.' %(time.clock() - t0))
        t0 = time.clock()

        if mode == 'i':
            recordsView = records.view()
            recordsView.shape = (nBuffersPerAcquisition,
                                 self.nChannels,
                                 nRecordsPerBuffer,
                                 samplesPerRecord)
            recordsView = np.rollaxis(recordsView,
                    1, 0).reshape(self.nChannels,
                                  nRecord,
                                  samplesPerRecord)

            if samplesPerRecord != nSamples:
                recordsView = recordsView[:,:,:nSamples]
            recordsFloat = recordsView.astype(dtype=np.float32, copy=True)
            recordsFloat -= self.codeZero

            recordsFloat[0] *= rangeA
            recordsFloat[1] *= rangeB
            data['Channel A'] = recordsFloat[0]
            data['Channel B'] = recordsFloat[1]

        elif mode == 'a':
            avgRecord.shape = (self.nChannels,
                               nRecordsPerBuffer,
                               samplesPerRecord)
            avgRecord = np.sum(avgRecord, axis=1)
            if samplesPerRecord != nSamples:
                avgRecord = avgRecord[:,:nSamples]

            avgRecordFloat = avgRecord.astype(dtype=np.float32, copy=True)
            avgRecordFloat /= np.float32(nRecord)
            avgRecordFloat -= self.codeZero
            avgRecordFloat[0] *= rangeA
            avgRecordFloat[1] *= rangeB
            
            data['Channel A - Average record'] = avgRecordFloat[0]
            data['Channel B - Average record'] = avgRecordFloat[1]

        elif mode == 'b':
            if expansion > 1:
                avgRecord.shape = (self.nChannels,
                                   expansion,
                                   int(nRecordsPerBuffer / expansion),
                                   samplesPerRecord)
                avgRecord = np.sum(avgRecord, axis=1) / float(expansion)
            else:
                avgRecord.shape = (self.nChannels,
                                   nRecordsPerBuffer,
                                   samplesPerRecord)

            if samplesPerRecord != nSamples:
                avgRecord = avgRecord[:,:,:nSamples]

            avgRecordFloat = avgRecord.astype(dtype=np.float32, copy=True)
            avgRecordFloat /= np.float32(nBuffersPerAcquisition)
            avgRecordFloat -= self.codeZero
            avgRecordFloat[0] *= rangeA
            avgRecordFloat[1] *= rangeB

            data['Channel A - Average buffer'] = avgRecordFloat[0]
            data['Channel B - Average buffer'] = avgRecordFloat[1]

        log.info('Pre-processing: %.6f s.' %(time.clock() - t0))

        return data

    def removeBuffersDMA(self):
        """Clear and remove DMA buffers to release memory."""
        # make sure buffers release memory
        for buf in self.buffers:
            buf.__exit__()
        # remove all
        self.buffers = []

    def readRecordsSinglePort(self, Channel):
        """Reads records, converts to float."""
        # define sizes
        samplesPerRecord = self.nPreSize + self.nPostSize
        dataBuffer = (c_uint8 * samplesPerRecord)()
        # define scale factors
        voltScale = self.dRange[Channel] / self.codeRange
        # initialize a scaled float vector
        records = np.empty((self.nRecord, samplesPerRecord), dtype=np.float32)
        for n in range(self.nRecord):
            self.AlazarRead(Channel, dataBuffer, self.bytesPerSample,
                            n+1, -self.nPreSize, samplesPerRecord)
            # convert and scale to float
            records[n] = np.array(dataBuffer, dtype=np.float32)
        records -= self.codeZero
        records *= voltScale
        return records


if __name__ == '__main__':
    # test driver
    Digitizer = AlazarTechDigitizer()
