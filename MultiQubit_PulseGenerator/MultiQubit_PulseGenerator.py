# Copyright 2014-2021 Keysight Technologies
#!/usr/bin/env python
import importlib
import os
import sys
import copy
import re

import numpy as np

from BaseDriver import LabberDriver
from sequence_builtin import (
    CPMG, PulseTrain, Rabi, SpinLocking, ReadoutTraining)
from sequence_rb import SingleQubit_RB, TwoQubit_RB
from sequence import SequenceToWaveforms
import logging
log = logging.getLogger('LabberDriver')

# dictionary with built-in sequences
SEQUENCES = {'Rabi': Rabi,
             'CP/CPMG': CPMG,
             'Pulse train': PulseTrain,
             '1-QB Randomized Benchmarking': SingleQubit_RB,
             '2-QB Randomized Benchmarking': TwoQubit_RB,
             'Spin-locking': SpinLocking,
             'Readout training': ReadoutTraining,
             'Custom': type(None)}


class Driver(LabberDriver):
    """This class implements a multi-qubit pulse generator."""

    def performOpen(self, options={}):
        """Perform the operation of opening the instrument connection."""
        # init variables
        self.n_qubits = int(self.getValue('Number of qubits'))
        self.sequence = None
        self.sequence_to_waveforms = SequenceToWaveforms(self.n_qubits)
        self.waveforms = {}
        # always create a sequence at startup
        name = self.getValue('Sequence')
        self.sendValueToOther('Sequence', name)

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        """Perform the Set Value instrument operation."""
        # only do something here if changing the sequence type
        if quant.name == 'Sequence':
            # create new sequence if sequence type changed
            new_type = SEQUENCES[value]

            if not isinstance(self.sequence, new_type):
                # create a new sequence object
                if value == 'Custom':
                    # for custom python files
                    path = self.getValue('Custom Python file')
                    (path, modName) = os.path.split(path)
                    sys.path.append(path)
                    modName = modName.split('.py')[0]  # strip suffix
                    self.mod = importlib.import_module(modName)
                    # the custom sequence class has to be named
                    # 'CustomSequence'
                    if not isinstance(self.sequence, self.mod.CustomSequence):
                        self.sequence = self.mod.CustomSequence(self.n_qubits)
                else:
                    # standard built-in sequence
                    self.sequence = new_type(self.n_qubits)

        elif (quant.name == 'Custom Python file' and
              self.getValue('Sequence') == 'Custom'):
            # for custom python files
            path = self.getValue('Custom Python file')
            (path, modName) = os.path.split(path)
            modName = modName.split('.py')[0]  # strip suffix
            sys.path.append(path)
            self.mod = importlib.import_module(modName)
            # the custom sequence class has to be named 'CustomSequence'
            if not isinstance(self.sequence, self.mod.CustomSequence):
                self.sequence = self.mod.CustomSequence(self.n_qubits)
                self.bCfgUpdated = True

        elif quant.name == 'Number of qubits':
            self.existing_n_qubits = self.getValue('Number of qubits')
            n_qubits = int(value)
            self.sequence.__init__(n_qubits)
            self.sequence_to_waveforms = SequenceToWaveforms(n_qubits)
            # force recalculation of the sequence in performGetValue
            self.bCfgUpdated = True

        elif (quant.name == 'Generate waveform primitives' and value and (self.getValue('Sample rate') != 1e9)) \
                or (quant.name == 'Sample rate' and value != 1e9 and self.getValue('Generate waveform primitives')):
            raise ValueError('Only 1GS/s is supported for waveform primitives.')

        return value

    def performGetValue(self, quant, options={}):
        """Perform the Get Value instrument operation."""
        # ignore if no sequence
        if self.sequence is None:
            return quant.getValue()

        if self.getValue('Sequence') == 'Custom':
            # reload module in case the sequence definition or name has changed
            custom_file = self.getValue('Custom Python file')
            if custom_file:
                modName = os.path.split(custom_file)[1]
                modName = modName.split('.py')[0]  # strip suffix
                if modName == self.mod.__name__:
                    importlib.reload(self.mod)
                else:
                    self.mod = importlib.import_module(modName)
                self.sequence = self.mod.CustomSequence(self.n_qubits)
                self.bCfgUpdated = True

        # check type of quantity
        if (quant.name.startswith('Voltage, QB') or
                quant.name.startswith('Single-shot, QB')):
            # perform demodulation, check if config is updated
            if self.isConfigUpdated():
                # update sequence object with current driver configuation
                config = self.instrCfg.getValuesDict()
                self.sequence.set_parameters(config)
                self.sequence_to_waveforms.set_parameters(config)
            # get qubit index and waveforms
            n = int(re.search(r'QB(\d+)', quant.name).group(1)) - 1
            demod_iq = self.getValue('Demodulation - IQ')
            if demod_iq:
                signal_i = self.getValue('Demodulation - Input I')
                signal_q = self.getValue('Demodulation - Input Q')
            else:
                signal = self.getValue('Demodulation - Input')
            ref = self.getValue('Demodulation - Reference')
            # perform demodulation
            if demod_iq:
                value = self.sequence_to_waveforms.readout.demodulate_iq(
                    n, signal_i, signal_q, ref)
            else:
                value = self.sequence_to_waveforms.readout.demodulate(
                    n, signal, ref)
            # average values if not single-shot
            if not quant.name.startswith('Single-shot, QB'):
                value = np.mean(value)

        elif quant.isVector():
            # traces, check if waveform needs to be re-calculated
            if self.isConfigUpdated():
                # update sequence object with current driver configuation
                config = self.instrCfg.getValuesDict()
                self.sequence.set_parameters(config)
                self.sequence_to_waveforms.set_parameters(config)

                # check if calculating multiple sequences, for randomization
                multi_rb = config.get('Output multiple sequences', False)
                multi_training = config.get('Train all states at once', False)

                if (multi_rb or multi_training):
                    # create multiple randomizations, store in memory
                    if multi_rb:
                        multi_param = 'Randomize'
                        align_multi_to_end = config.get(
                            'Align RB waveforms to end', False)
                        n_call = int(
                            config.get('Number of multiple sequences', 1))
                    elif multi_training:
                        # readout training
                        multi_param = 'Training, input state'
                        align_multi_to_end = False
                        # for readout, always start with first state
                        config[multi_param] = -1
                        training_type = config['Training type']
                        if training_type == 'Specific qubit':
                            n_call = 2
                        elif training_type == 'All qubits at once':
                            n_call = 2
                        elif training_type == 'All combinations':
                            n_call = 2**self.sequence.n_qubit

                    calls = []
                    for n in range(n_call):
                        config[multi_param] += 1
                        calls.append(
                            copy.deepcopy(
                                self.sequence_to_waveforms.get_waveforms(
                                    self.sequence.get_sequence(config))))

                    # after all calls are done, convert output to matrix form
                    self.waveforms = dict()
                    n_qubit = self.sequence.n_qubit
                    # start with xy, z and gate waveforms, list of data
                    for key in ['xy', 'z', 'gate', 'readout_iq']:
                        # get size of longest waveform
                        self.waveforms[key] = []
                        for n in range(n_qubit):
                            log.info('Generating {} waveform for qubit {}'.format(key, n))

                            length = max([len(call[key][n]) for call in calls])
                            # build matrix
                            datatype = calls[0][key][n].dtype
                            data = np.zeros((n_call, length), dtype=datatype)
                            for m, call in enumerate(calls):
                                if align_multi_to_end:
                                    data[m][-len(call[key][n]):] = call[key][n]
                                else:
                                    data[m][:len(call[key][n])] = call[key][n]
                            self.waveforms[key].append(data)

                    # same for readout trigger
                    key = 'readout_trig'
                    length = max([len(call[key]) for call in calls])
                    datatype = calls[0][key].dtype
                    data = np.zeros((n_call, length), dtype=datatype)
                    for m, call in enumerate(calls):
                        if align_multi_to_end:
                            data[m][-len(call[key]):] = call[key]
                        else:
                            data[m][:len(call[key])] = call[key]
                    self.waveforms[key] = data

                else:
                    # normal operation, calculate waveforms
                    # log.info('generating case 2')
                    self.waveforms = self.sequence_to_waveforms.get_waveforms(
                        self.sequence.get_sequence(config))
                    # log.info('Z waveform max: {}'.format(np.max(self.waveforms['z'])))
            # get correct data from waveforms stored in memory
            value = self.getWaveformFromMemory(quant)
        else:
            # for all other cases, do nothing
            value = quant.getValue()
        return value

    def getWaveformFromMemory(self, quant):
        """Return data from already calculated waveforms."""
        # check which data to return
        if quant.name.startswith('Trace - I'):
            n = int(re.search(r'I(\d+)', quant.name).group(1)) - 1
            if self.getValue('Swap IQ'):
                value = self.waveforms['xy'][n].imag
            else:
                value = self.waveforms['xy'][n].real
        elif quant.name.startswith('Trace - Q'):
            n = int(re.search(r'Q(\d+)', quant.name).group(1)) - 1
            if self.getValue('Swap IQ'):
                value = self.waveforms['xy'][n].real
            else:
                value = self.waveforms['xy'][n].imag
        elif quant.name.startswith('Trace - Z'):
            n = int(re.search(r'Z(\d+)', quant.name).group(1)) - 1
            value = self.waveforms['z'][n]
        elif quant.name.startswith('Trace - G'):
            n = int(re.search(r'G(\d+)', quant.name).group(1)) - 1
            value = self.waveforms['gate'][n]

        elif quant.name == 'Trace - Readout trig':
            value = self.waveforms['readout_trig']
        elif quant.name.startswith('Trace - Readout I'):
            index = int(re.search(r'Readout I(\d+)', quant.name).group(1))
            value = self.waveforms['readout_iq'][index - 1].real
        elif quant.name.startswith('Trace - Readout Q'):
            index = int(re.search(r'Readout Q(\d+)', quant.name).group(1))
            value = self.waveforms['readout_iq'][index - 1].imag

        # primitives
        elif quant.name == 'Waveform primitives':
            value = self.sequence_to_waveforms.primitives.pack_data()
        elif quant.name == 'Control sequence':
            value = self.sequence_to_waveforms.control_sequence.pack_data()

        # return data as dict with sampling information
        dt = 1 / self.sequence_to_waveforms.sample_rate
        value = quant.getTraceDict(value, dt=dt)
        return value


if __name__ == '__main__':
    pass
