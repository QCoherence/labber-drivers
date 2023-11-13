# Copyright 2014-2021 Keysight Technologies
"""This module allows the user to create a new MultiQubit_PulseGenerator driver for an arbitrary number of qubits
between 1 and 64, as well as optionally specify the 2-qubit gate pairs to be included. Specifying the number of
qubits requires supplying an integer value for n_qubits. Specifying the two-qubit gate pairs is done by providing
 a path to a json file with a list of integer qubit pairs. If no custom two-qubit gate file is supplied the default
 gate pairs are nearest neighbors in a line, where qubit n is connected to qubit n-1 and qubit n+1

Usage:
python generate_mqpg.py --n_qubits 64 --gates_file connected_qubits.json
"""

import configparser
import json
import os
import sys
import argparse
import Labber
import ScriptsAndSettings


class MQPGConfig:
    """A class with methods to create a configparser, use it to generate all of the parameters for a
     MultiQubit_PulseGenerator. and write the file back

     Parameters
     --------------
     n_qubits: int
        Number of qubits to be controlled by the pulse generator. Between 1 and 64, default 16.

    connected_qubits: List[List[int]]]
        A list of integer pairs corresponding to connected qubit pairs for 2-qubit gates
     """

    def __init__(self, n_qubits=16, connected_qubits=None):

        if n_qubits > 64:
            raise ValueError("More than 64 qubits not supported")

        self.config = configparser.ConfigParser()
        self.n_qubits = n_qubits
        self.connected_qubits = connected_qubits if connected_qubits else [[x, x + 1] for x in range(1, self.n_qubits)]
        self.state_values = [(f'state_value_{ind + 1}', str(64 - ind)) for ind in range(64)]


    def generate(self):
        self._general_settings()
        self._sequence()
        self._waveform()
        self._primitives()
        self._one_qb_gates_xy()
        self._one_qb_gates_z()
        self._global_z_offset()
        self._z_pulse_during_readout()
        self._two_qb_gates()
        self._qb_spectra()
        self._tomography()
        self._predistortion()
        self._crosstalk()
        self._readout()
        self._output()
        self._demodulation()

    def write(self, output_file):
        with open(output_file, 'w') as configfile:
            self.config.write(configfile)

    def _general_settings(self):
        self.config["General settings"] = {'name': 'Multi-Qubit Pulse Generator',
                                           'version': '1.4',
                                           'driver_path': 'MultiQubit_PulseGenerator',
                                           'signal_generator': 'True',
                                           'signal_analyzer': 'True'}

    def _sequence(self):
        self.config['Sequence'] = {'datatype': 'COMBO',
                                   'combo_def_1': 'Rabi',
                                   'combo_def_2': 'CP/CPMG',
                                   'combo_def_3': 'Pulse train',
                                   'combo_def_4': '1-QB Randomized Benchmarking',
                                   'combo_def_5': '2-QB Randomized Benchmarking',
                                   'combo_def_6': 'Spin-locking',
                                   'combo_def_7': 'Readout training',
                                   'combo_def_8': 'Custom',
                                   'group': 'Sequence',
                                   'section': 'Sequence',
                                   'show_in_measurement_dlg': 'True'}

        self.config["Number of qubits"] = {'datatype': 'COMBO',
                                           'def_value': '{}'.format(self.n_qubits),
                                           'group': 'Sequence',
                                           'section': 'Sequence'}

        for qubit in range(1, self.n_qubits + 1):
            self.config.set("Number of qubits", "combo_def_{}".format(qubit), "{}".format(qubit))

        self.config["Custom Python file"] = {'datatype': 'PATH',
                                             'tooltip': 'Python file with custom Sequence class. The class must be named "CustomSequence"',
                                             'state_quant': 'Sequence',
                                             'state_value_1': 'Custom',
                                             'group': 'Sequence',
                                             'section': 'Sequence'}

        self.config['Pulse spacing'] = {'datatype': 'DOUBLE',
                                        'unit': 's',
                                        'def_value': '20E-9',
                                        'low_lim': '0',
                                        'group': 'Sequence',
                                        'section': 'Sequence',
                                        'show_in_measurement_dlg': 'True'}

        self.config['Local XY Control'] = {'datatype': 'BOOLEAN',
                                           'def_value': 'True',
                                           'group': 'Sequence',
                                           'section': 'Sequence'}

        self.config['Dynamic duration'] = {'datatype': 'BOOLEAN',
                                           'def_value': 'False',
                                           'group': 'Sequence',
                                           'section': 'Sequence',
                                           'tooltip': 'Requires Pathwave TSE'}

        self.config["Pulse"] = {'datatype': 'COMBO',
                                'combo_def_1': 'Xp',
                                'combo_def_2': 'X2p',
                                'combo_def_3': 'Yp',
                                'combo_def_4': 'Y2p',
                                'combo_def_5': 'Zp',
                                'combo_def_6': 'Z2p',
                                'combo_def_7': 'CPh',
                                'state_quant': 'Sequence',
                                'state_value_1': 'Pulse train',
                                'group': 'Pulse train',
                                'section': 'Sequence'}

        self.config["# of pulses"] = {'datatype': 'DOUBLE',
                                      'def_value': '1',
                                      'low_lim': '0',
                                      'state_quant': 'Sequence',
                                      'state_value_1': 'Pulse train',
                                      'group': 'Pulse train',
                                      'section': 'Sequence'}

        self.config["Alternate pulse direction"] = {'datatype': 'BOOLEAN',
                                                    'def_value': 'False',
                                                    'state_quant': 'Sequence',
                                                    'state_value_1': 'Pulse train',
                                                    'group': 'Pulse train',
                                                    'section': 'Sequence'}

        self.config["# of pi pulses"] = {'datatype': 'DOUBLE',
                                         'tooltip': 'Number of refocusing pulses. Set to -1 for generating pulse sequences for T1 experiments',
                                         'def_value': '1',
                                         'low_lim': '-1',
                                         'state_quant': 'Sequence',
                                         'state_value_1': 'CP/CPMG',
                                         'group': 'CPMG',
                                         'section': 'Sequence'}

        self.config["Sequence duration"] = {'datatype': 'DOUBLE',
                                            'unit': 's',
                                            'def_value': '1E-6',
                                            'low_lim': '0',
                                            'state_quant': 'Sequence',
                                            'state_value_1': 'CP/CPMG',
                                            'group': 'CPMG',
                                            'section': 'Sequence'}

        self.config["Add pi pulses to Q"] = {'datatype': 'BOOLEAN',
                                             'def_value': 'True',
                                             'state_quant': 'Sequence',
                                             'state_value_1': 'CP/CPMG',
                                             'group': 'CPMG',
                                             'section': 'Sequence'}

        self.config["Edge-to-edge pulses"] = {'datatype': 'BOOLEAN',
                                              'def_value': 'False',
                                              'state_quant': 'Sequence',
                                              'state_value_1': 'CP/CPMG',
                                              'group': 'CPMG',
                                              'section': 'Sequence'}

        self.config["Add last pi/2 pulse to Q"] = {'datatype': 'BOOLEAN',
                                                   'def_value': 'False',
                                                   'state_quant': 'Sequence',
                                                   'state_value_1': 'CP/CPMG',
                                                   'group': 'CPMG',
                                                   'section': 'Sequence'}

        self.config["Last pulse phase"] = {'datatype': 'DOUBLE',
                                           'unit': 'deg',
                                           'def_value': '0',
                                           'state_quant': 'Sequence',
                                           'state_value_1': 'CP/CPMG',
                                           'group': 'CPMG',
                                           'section': 'Sequence',
                                           'tooltip': 'Phase of second pi/2 pulse'}

        self.config["Dynamic last pulse phase"] = {'datatype': 'BOOLEAN',
                                                   'def_value': 'False',
                                                   'state_quant': 'Sequence',
                                                   'state_value_1': 'CP/CPMG',
                                                   'group': 'CP/CPMG - Hardware-dynamic parameters',
                                                   'section': 'Sequence',
                                                   'tooltip': 'Requires Pathwave TSE'}

        self.config["Drive pulse duration"] = {'datatype': 'DOUBLE',
                                               'unit': 's',
                                               'def_value': '1E-6',
                                               'low_lim': '0',
                                               'state_quant': 'Sequence',
                                               'state_value_1': 'Spin-locking',
                                               'group': 'Spin-locking',
                                               'section': 'Sequence'}

        self.config["Drive pulse phase"] = {'datatype': 'DOUBLE',
                                            'unit': 'deg',
                                            'def_value': '0',
                                            'state_quant': 'Sequence',
                                            'state_value_1': 'Spin-locking',
                                            'group': 'Spin-locking',
                                            'section': 'Sequence'}

        self.config["Pulse sequence"] = {'datatype': 'COMBO',
                                         'combo_def_1': 'SL-3',
                                         'combo_def_2': 'SL-5a',
                                         'combo_def_3': 'SL-5b',
                                         'def_value': 'SL-3',
                                         'state_quant': 'Sequence',
                                         'state_value_1': 'Spin-locking',
                                         'group': 'Spin-locking',
                                         'section': 'Sequence'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Drive pulse amplitude #{}".format(qubit)] = {'datatype': 'DOUBLE',
                                                                      'unit': 'V',
                                                                      'def_value': '0',
                                                                      'state_quant': 'Sequence',
                                                                      'state_value_1': 'Spin-locking',
                                                                      'group': 'Spin-locking',
                                                                      'section': 'Sequence',
                                                                      'tooltip': 'Drive pulse amplitude for qubit #{}'.format(
                                                                          qubit)}



        self.config["Number of Cliffords"] = {'datatype': 'DOUBLE',
                                              'def_value': '17',
                                              'low_lim': '0',
                                              'group': 'Randomized Benchmarking',
                                              'section': 'Sequence',
                                              'state_quant': 'Sequence',
                                              'state_value_1': '1-QB Randomized Benchmarking',
                                              'state_value_2': '2-QB Randomized Benchmarking',
                                              'show_in_measurement_dlg': 'True'}

        self.config["Randomize"] = {'datatype': 'DOUBLE',
                                    'def_value': '0',
                                    'group': 'Randomized Benchmarking',
                                    'tooltip': 'Randomize gate sequence only if "Randomize" value changes or "Number of Cliffords" changes',
                                    'section': 'Sequence',
                                    'state_quant': 'Sequence',
                                    'state_value_1': '1-QB Randomized Benchmarking',
                                    'state_value_2': '2-QB Randomized Benchmarking',
                                    'show_in_measurement_dlg': 'True'}

        self.config["Write sequence as txt file"] = {'datatype': 'BOOLEAN',
                                                     'def_value': '0',
                                                     'group': 'Randomized Benchmarking',
                                                     'section': 'Sequence',
                                                     'tooltip': 'Save txt file for the sequence generated',
                                                     'state_quant': 'Sequence',
                                                     'state_value_1': '1-QB Randomized Benchmarking',
                                                     'state_value_2': '2-QB Randomized Benchmarking'}

        self.config["Output multiple sequences"] = {'datatype': 'BOOLEAN',
                                                    'def_value': '0',
                                                    'group': 'Randomized Benchmarking',
                                                    'tooltip': 'Output multiple randomized gate sequences.',
                                                    'section': 'Sequence',
                                                    'state_quant': 'Sequence',
                                                    'state_value_1': '1-QB Randomized Benchmarking',
                                                    'state_value_2': '2-QB Randomized Benchmarking',
                                                    'show_in_measurement_dlg': 'True'}

        self.config["Number of multiple sequences"] = {'datatype': 'DOUBLE',
                                                       'def_value': '1',
                                                       'low_lim': '1',
                                                       'group': 'Randomized Benchmarking',
                                                       'tooltip': 'Number of multiple randomize gate sequences to output.',
                                                       'section': 'Sequence',
                                                       'state_quant': 'Output multiple sequences',
                                                       'state_value_1': '1',
                                                       'show_in_measurement_dlg': 'True'}

        self.config["Align RB waveforms to end"] = {'datatype': 'BOOLEAN',
                                                    'def_value': '1',
                                                    'group': 'Randomized Benchmarking',
                                                    'section': 'Sequence',
                                                    'state_quant': 'Output multiple sequences',
                                                    'state_value_1': '1',
                                                    'show_in_measurement_dlg': 'True'}

        self.config["Simultaneous pulses"] = {'datatype': 'BOOLEAN',
                                              'tooltip': 'If False, all pulses are separated in time.',
                                              'def_value': '1',
                                              'group': 'Randomized Benchmarking',
                                              'section': 'Sequence',
                                              'state_quant': 'Sequence',
                                              'state_value_1': '1-QB Randomized Benchmarking',
                                              'state_value_2': '2-QB Randomized Benchmarking'}

        self.config['Qubit 1 to Benchmark'] = {'datatype': 'COMBO',
                                               'def_value': '1',
                                               'state_quant': 'Sequence',
                                               'state_value_1': '2-QB Randomized Benchmarking',
                                               'group': 'Randomized Benchmarking',
                                               'section': 'Sequence',
                                               'show_in_measurement_dlg': 'True'}

        for qubit in range(1, self.n_qubits + 1):
            self.config.set("Qubit 1 to Benchmark", "combo_def_{}".format(qubit), "{}".format(qubit))

        if self.n_qubits > 1:
            self.config['Qubit 2 to Benchmark'] = {'datatype': 'COMBO',
                                                   'def_value': '2',
                                                   'state_quant': 'Sequence',
                                                   'state_value_1': '2-QB Randomized Benchmarking',
                                                   'group': 'Randomized Benchmarking',
                                                   'section': 'Sequence',
                                                   'show_in_measurement_dlg': 'True'}

            self.config["Native 2-QB gate"] = {'datatype': 'COMBO',
                                               'def_value': 'CZ',
                                               'group': 'Randomized Benchmarking',
                                               'combo_def_1': 'CZ',
                                               'combo_def_2': 'iSWAP',
                                               'section': 'Sequence',
                                               'state_quant': 'Sequence',
                                               'state_value_1': '2-QB Randomized Benchmarking',
                                               'show_in_measurement_dlg': 'True'}

            for qubit in range(1, self.n_qubits + 1):
                self.config.set("Qubit 2 to Benchmark", "combo_def_{}".format(qubit), "{}".format(qubit))

        self.config["Interleave 1-QB Gate"] = {'datatype': 'BOOLEAN',
                                               'def_value': '0',
                                               'state_quant': 'Sequence',
                                               'state_value_1': '1-QB Randomized Benchmarking',
                                               'group': 'Randomized Benchmarking',
                                               'section': 'Sequence',
                                               'show_in_measurement_dlg': 'True'}

        self.config["Interleaved 1-QB Gate"] = {'datatype': 'COMBO',
                                                'combo_def_1': 'Ref',
                                                'combo_def_2': 'I',
                                                'combo_def_3': 'Xp',
                                                'combo_def_4': 'Xm',
                                                'combo_def_5': 'X2p',
                                                'combo_def_6': 'X2m',
                                                'combo_def_7': 'Yp',
                                                'combo_def_8': 'Ym',
                                                'combo_def_9': 'Y2p',
                                                'combo_def_10': 'Y2m',
                                                'combo_def_11': 'VZp',
                                                'combo_def_12': 'Zp',
                                                'combo_def_13': 'Ilong',
                                                'state_quant': 'Interleave 1-QB Gate',
                                                'state_value_1': '1',
                                                'group': 'Randomized Benchmarking',
                                                'section': 'Sequence',
                                                'show_in_measurement_dlg': 'True'}
        if self.n_qubits > 1:
            self.config["Interleave 2-QB Gate"] = {'datatype': 'BOOLEAN',
                                               'def_value': '0',
                                               'state_quant': 'Sequence',
                                               'state_value_1': '2-QB Randomized Benchmarking',
                                               'group': 'Randomized Benchmarking',
                                               'section': 'Sequence',
                                               'show_in_measurement_dlg': 'True'}

            self.config["Interleaved 2-QB Gate"] = {'datatype': 'COMBO',
                                                'combo_def_1': 'CZ',
                                                'combo_def_2': 'CZEcho',
                                                'combo_def_3': 'iSWAP',
                                                'combo_def_4': 'I',
                                                'state_quant': 'Interleave 2-QB Gate',
                                                'state_value_1': '1',
                                                'group': 'Randomized Benchmarking',
                                                'section': 'Sequence',
                                                'show_in_measurement_dlg': 'True'}

            self.config["Find the cheapest recovery Clifford"] = {'datatype': 'BOOLEAN',
                                                              'def_value': '1',
                                                              'state_quant': 'Sequence',
                                                              'state_value_1': '2-QB Randomized Benchmarking',
                                                              'group': 'Randomized Benchmarking',
                                                              'section': 'Sequence',
                                                              'show_in_measurement_dlg': 'True'}

            self.config["Use a look-up table"] = {'datatype': 'BOOLEAN', 'def_value': '1',
                                              'state_quant': 'Find the cheapest recovery Clifford',
                                              'state_value_1': '1', 'group': 'Randomized Benchmarking',
                                              'section': 'Sequence', 'show_in_measurement_dlg': 'True'}

            self.config["File path of the look-up table"] = {'datatype': 'PATH', 'state_quant': 'Use a look-up table',
                                                         'state_value_1': '1', 'group': 'Randomized Benchmarking',
                                                         'section': 'Sequence', 'show_in_measurement_dlg': 'True'}

        self.config["Training type"] = {'datatype': 'COMBO', 'def_value': 'All combinations',
                                        'combo_def_1': 'Specific qubit', 'combo_def_2': 'All qubits at once',
                                        'combo_def_3': 'All combinations', 'state_quant': 'Sequence',
                                        'state_value_1': 'Readout training',
                                        'tooltip': 'For "All combinations", (Number of states)^(Number of qubits) training data sets are required.  For other options, (Number of states) training sets are required',
                                        'group': 'Readout training', 'section': 'Sequence',
                                        'show_in_measurement_dlg': 'True'}

        self.config["Training, qubit"] = {'datatype': 'DOUBLE', 'state_quant': 'Training type',
                                          'state_value_1': 'Specific qubit', 'def_value': '1', 'low_lim': '-1',
                                          'group': 'Readout training', 'section': 'Sequence',
                                          'show_in_measurement_dlg': 'True'}

        self.config["Training, input state"] = {'datatype': 'DOUBLE', 'def_value': '0', 'low_lim': '-1',
                                                'state_quant': 'Sequence', 'state_value_1': 'Readout training',
                                                'group': 'Readout training', 'section': 'Sequence',
                                                'show_in_measurement_dlg': 'True'}

        self.config["Train all states at once"] = {'datatype': 'BOOLEAN', 'def_value': '0', 'group': 'Readout training',
                                                   'section': 'Sequence',
                                                   'tooltip': 'If checked, the driver will create training waveforms for all possible states',
                                                   'state_quant': 'Sequence', 'state_value_1': 'Readout training',
                                                   'show_in_measurement_dlg': 'True'}

        for param_index in range(1, 11):
            self.config["Parameter #{}".format(param_index)] = {'datatype': 'DOUBLE', 'def_value': '0.0',
                                                                'state_quant': 'Sequence', 'state_value_1': 'Custom',
                                                                'group': 'Custom', 'section': 'Sequence'}

    def _waveform(self):

        self.config['Sample rate'] = {'datatype': 'DOUBLE',
                                      'def_value': '1.0E9',
                                      'unit': 'Hz',
                                      'group': 'Waveform',
                                      'section': 'Waveform',
                                      'show_in_measurement_dlg': 'True'}

        self.config['Trim waveform to sequence'] = {'label': 'Adjust waveform to sequence length',
                                                    'tooltip': 'Automatically adjusts the number of points in the waveform',
                                                    'datatype': 'BOOLEAN',
                                                    'def_value': '1',
                                                    'group': 'Waveform',
                                                    'section': 'Waveform'}

        self.config['Number of points'] = {'datatype': 'DOUBLE',
                                           'def_value': '240E3',
                                           'state_quant': 'Trim waveform to sequence',
                                           'state_value': 'False',
                                           'group': 'Waveform',
                                           'section': 'Waveform',
                                           'show_in_measurement_dlg': 'True'}

        self.config['Align pulses to end of waveform'] = {'datatype': 'BOOLEAN',
                                                          'state_quant': 'Trim waveform to sequence',
                                                          'state_value': 'False',
                                                          'def_value': '0',
                                                          'group': 'Waveform',
                                                          'section': 'Waveform'}

        self.config["First pulse delay"] = {'datatype': 'DOUBLE',
                                            'unit': 's',
                                            'def_value': '100E-9',
                                            'group': 'Delays',
                                            'section': 'Waveform'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Qubit {} XY Delay".format(qubit)] = {'datatype': 'DOUBLE',
                                                              'def_value': '0',
                                                              'unit': 's',
                                                              'group': 'Delays',
                                                              'section': 'Waveform'}

            self.config["Qubit {} Z Delay".format(qubit)] = {'datatype': 'DOUBLE',
                                                             'def_value': '0',
                                                             'unit': 's',
                                                             'group': 'Delays',
                                                             'section': 'Waveform'}

    def _primitives(self):
        self.config['Generate waveform primitives'] = {'datatype': 'BOOLEAN',
                                                       'def_value': 'False',
                                                       'group': 'Primitives',
                                                       'section': 'Primitives',
                                                       'tooltip': 'Requires Pathwave TSE'}

        self.config["Apply phase modulation to qubit primitives"] = {'datatype': 'BOOLEAN',
                                                                     'def_value': 'False',
                                                                     'group': 'Qubit primitives',
                                                                     'section': 'Primitives',
                                                                     'state_quant': 'Generate waveform primitives',
                                                                     'state_value': 'True'}

        self.config["Apply amplitude modulation to qubit primitives"] = {'datatype': 'BOOLEAN',
                                                                         'def_value': 'False',
                                                                         'group': 'Qubit primitives',
                                                                         'section': 'Primitives',
                                                                         'state_quant': 'Generate waveform primitives',
                                                                         'state_value': 'True'}

        self.config["Apply frequency modulation to readout primitives"] = {'datatype': 'BOOLEAN',
                                                                           'def_value': 'True',
                                                                           'group': 'Readout primitives',
                                                                           'section': 'Primitives',
                                                                           'state_quant': 'Generate waveform primitives',
                                                                           'state_value': 'True'}

        self.config["Apply phase modulation to readout primitives"] = {'datatype': 'BOOLEAN',
                                                                       'def_value': 'True',
                                                                       'group': 'Readout primitives',
                                                                       'section': 'Primitives',
                                                                       'state_quant': 'Generate waveform primitives',
                                                                       'state_value': 'True'}

        self.config["Apply amplitude modulation to readout primitives"] = {'datatype': 'BOOLEAN',
                                                                           'def_value': 'True',
                                                                           'group': 'Readout primitives',
                                                                           'section': 'Primitives',
                                                                           'state_quant': 'Generate waveform primitives',
                                                                           'state_value': 'True'}

        self.config["Enable hardware-dynamic amplitudes"] = {'datatype': 'BOOLEAN',
                                                             'def_value': 'False',
                                                             'group': 'Hardware-dynamic parameters',
                                                             'section': 'Primitives',
                                                             'state_quant': 'Generate waveform primitives',
                                                             'state_value': 'True',
                                                             'tooltip': 'Requires Pathwave TSE'}

        self.config["Enable hardware-dynamic frequencies"] = {'datatype': 'BOOLEAN',
                                                              'def_value': 'False',
                                                              'group': 'Hardware-dynamic parameters',
                                                              'section': 'Primitives',
                                                              'state_quant': 'Generate waveform primitives',
                                                              'state_value': 'True',
                                                              'tooltip': 'Requires Pathwave TSE'}

    def _one_qb_gates_xy(self):
        self.config["Pulse type"] = {'datatype': 'COMBO',
                                     'combo_def_1': 'Gaussian',
                                     'combo_def_2': 'Square',
                                     'combo_def_3': 'Ramp',
                                     'combo_def_4': 'Cosine',
                                     'group': 'Pulse settings',
                                     'section': '1-QB gates XY'}

        self.config["Truncation range"] = {'datatype': 'DOUBLE', 'def_value': '3', 'state_quant': 'Pulse type',
                                           'state_value_1': 'Gaussian', 'group': 'Pulse settings',
                                           'section': '1-QB gates XY'}

        self.config["Start at zero"] = {'datatype': 'BOOLEAN', 'def_value': '0', 'state_quant': 'Pulse type',
                                        'state_value_1': 'Gaussian', 'group': 'Pulse settings',
                                        'section': '1-QB gates XY'}

        self.config["Use DRAG"] = {'datatype': 'BOOLEAN', 'def_value': 'False', 'group': 'Pulse settings',
                                   'section': '1-QB gates XY'}

        self.config["Uniform amplitude"] = {'label': 'Uniform amplitude', 'datatype': 'BOOLEAN', 'def_value': 'False',
                                            'group': 'Pulse settings', 'section': '1-QB gates XY'}

        self.config["Amplitude"] = {'datatype': 'DOUBLE', 'unit': 'V', 'def_value': '1.0',
                                    'state_quant': 'Uniform amplitude', 'state_value': 'True',
                                    'group': 'Pulse settings', 'section': '1-QB gates XY',
                                    'show_in_measurement_dlg': 'True'}

        self.config["Uniform pulse shape"] = {'datatype': 'BOOLEAN', 'def_value': 'True', 'group': 'Pulse settings',
                                              'section': '1-QB gates XY'}

        self.config["Width"] = {'label': 'Width', 'datatype': 'DOUBLE', 'unit': 's', 'def_value': '10E-9',
                                'group': 'Pulse settings', 'section': '1-QB gates XY',
                                'state_quant': 'Uniform pulse shape', 'state_value': '1',
                                'show_in_measurement_dlg': 'True'}

        self.config["Plateau"] = {'label': 'Plateau', 'datatype': 'DOUBLE', 'unit': 's', 'group': 'Pulse settings',
                                  'section': '1-QB gates XY', 'state_quant': 'Uniform pulse shape', 'state_value': '1',
                                  'show_in_measurement_dlg': 'True'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Amplitude #{}".format(qubit)] = {'label': 'Amplitude', 'datatype': 'DOUBLE', 'unit': 'V',
                                                          'def_value': '1.0', 'state_quant': 'Uniform amplitude',
                                                          'state_value': 'False', 'group': 'Pulse #{}'.format(qubit),
                                                          'section': '1-QB gates XY', 'show_in_measurement_dlg': 'True'}
            self.config["Width #{}".format(qubit)] = {'label': 'Width', 'datatype': 'DOUBLE', 'unit': 's',
                                                      'def_value': '10E-9', 'group': 'Pulse #{}'.format(qubit),
                                                      'section': '1-QB gates XY', 'state_quant': 'Uniform pulse shape',
                                                      'state_value': '0'}
            self.config["Plateau #{}".format(qubit)] = {'label': 'Plateau', 'datatype': 'DOUBLE', 'unit': 's',
                                                        'group': 'Pulse #{}'.format(qubit), 'section': '1-QB gates XY',
                                                        'state_quant': 'Uniform pulse shape', 'state_value': '0'}
            self.config["Frequency #{}".format(qubit)] = {'label': 'Frequency', 'datatype': 'DOUBLE', 'unit': 'Hz',
                                                          'group': 'Pulse #{}'.format(qubit),
                                                          'section': '1-QB gates XY'}
            self.config["DRAG scaling #{}".format(qubit)] = {'label': 'DRAG scaling', 'datatype': 'DOUBLE', 'unit': 's',
                                                             'def_value': '.25E-9', 'state_quant': 'Use DRAG',
                                                             'state_value': '1', 'group': 'Pulse #{}'.format(qubit),
                                                             'section': '1-QB gates XY'}
            self.config["DRAG frequency detuning #{}".format(qubit)] = {'label': 'DRAG frequency detuning',
                                                                        'datatype': 'DOUBLE', 'unit': 'Hz',
                                                                        'def_value': '0', 'state_quant': 'Use DRAG',
                                                                        'state_value': '1',
                                                                        'group': 'Pulse #{}'.format(qubit),
                                                                        'section': '1-QB gates XY'}
            self.config["Hardware-dynamic XY amplitude #{}".format(qubit)] = {'label': 'Hardware-dynamic amplitude',
                                                                           'datatype': 'BOOLEAN', 'def_value': 'False',
                                                                           'state_quant': 'Enable hardware-dynamic amplitudes',
                                                                           'state_value': 'True',
                                                                           'group': 'Pulse #{}'.format(qubit),
                                                                           'section': '1-QB gates XY',
                                                                           'show_in_measurement_dlg': 'False'}
            self.config["Hardware-dynamic frequency #{}".format(qubit)] = {'label': 'Hardware-dynamic frequency',
                                                                           'datatype': 'BOOLEAN', 'def_value': 'False',
                                                                           'state_quant': 'Enable hardware-dynamic frequencies',
                                                                           'state_value': 'True',
                                                                           'group': 'Pulse #{}'.format(qubit),
                                                                           'section': '1-QB gates XY',
                                                                           'show_in_measurement_dlg': 'False'}

    def _one_qb_gates_z(self):

        self.config["Pulse type, Z"] = {'label': 'Pulse type', 'datatype': 'COMBO', 'combo_def_1': 'Gaussian',
                                        'combo_def_2': 'Square', 'combo_def_3': 'Ramp', 'combo_def_4': 'Cosine',
                                        'group': 'Pulse settings', 'section': '1-QB gates Z'}

        self.config["Truncation range, Z"] = {'label': 'Truncation range', 'datatype': 'DOUBLE', 'def_value': '3',
                                              'state_quant': 'Pulse type, Z', 'state_value_1': 'Gaussian',
                                              'group': 'Pulse settings', 'section': '1-QB gates Z'}

        self.config["Start at zero, Z"] = {'label': 'Start at zero', 'datatype': 'BOOLEAN', 'def_value': '0',
                                           'state_quant': 'Pulse type, Z', 'state_value_1': 'Gaussian',
                                           'group': 'Pulse settings', 'section': '1-QB gates Z'}

        self.config["Uniform amplitude, Z"] = {'label': 'Uniform amplitude', 'datatype': 'BOOLEAN',
                                               'def_value': 'False', 'group': 'Pulse settings',
                                               'section': '1-QB gates Z'}

        self.config["Amplitude, Z"] = {'label': 'Amplitude', 'datatype': 'DOUBLE', 'unit': 'V', 'def_value': '1.0',
                                       'state_quant': 'Uniform amplitude, Z', 'state_value': 'True',
                                       'group': 'Pulse settings', 'section': '1-QB gates Z',
                                       'show_in_measurement_dlg': 'True'}

        self.config["Uniform pulse shape, Z"] = {'datatype': 'BOOLEAN', 'def_value': 'True', 'group': 'Pulse settings',
                                                 'section': '1-QB gates Z'}

        self.config["Width, Z"] = {'label': 'Width', 'datatype': 'DOUBLE', 'unit': 's', 'def_value': '10E-9',
                                   'group': 'Pulse settings', 'section': '1-QB gates Z',
                                   'state_quant': 'Uniform pulse shape, Z', 'state_value': '1',
                                   'show_in_measurement_dlg': 'True'}

        self.config["Plateau, Z"] = {'label': 'Plateau', 'datatype': 'DOUBLE', 'unit': 's', 'group': 'Pulse settings',
                                     'section': '1-QB gates Z', 'state_quant': 'Uniform pulse shape, Z',
                                     'state_value': '1', 'show_in_measurement_dlg': 'True'}


        for qubit in range(1, self.n_qubits + 1):
            self.config["Amplitude #{}, Z".format(qubit)] = {'label': 'Amplitude', 'datatype': 'DOUBLE', 'unit': 'V',
                                                             'def_value': '1.0', 'state_quant': 'Uniform amplitude, Z',
                                                             'state_value': 'False', 'group': 'Pulse #{}'.format(qubit),
                                                             'section': '1-QB gates Z',
                                                             'show_in_measurement_dlg': 'True'}
            self.config["Width #{}, Z".format(qubit)] = {'label': 'Width', 'datatype': 'DOUBLE', 'unit': 's',
                                                         'def_value': '10E-9', 'group': 'Pulse #{}'.format(qubit),
                                                         'section': '1-QB gates Z',
                                                         'state_quant': 'Uniform pulse shape, Z', 'state_value': '0'}
            self.config["Plateau #{}, Z".format(qubit)] = {'label': 'Plateau', 'datatype': 'DOUBLE', 'unit': 's',
                                                           'group': 'Pulse #{}'.format(qubit),
                                                           'section': '1-QB gates Z',
                                                           'state_quant': 'Uniform pulse shape, Z', 'state_value': '0'}

            self.config["Hardware-dynamic Z amplitude #{}".format(qubit)] = {'label': 'Hardware-dynamic amplitude',
                                                                           'datatype': 'BOOLEAN', 'def_value': 'False',
                                                                           'state_quant': 'Enable hardware-dynamic amplitudes',
                                                                           'state_value': 'True',
                                                                           'group': 'Pulse #{}'.format(qubit),
                                                                           'section': '1-QB gates Z',
                                                                           'show_in_measurement_dlg': 'False'}
    def _global_z_offset(self):
        self.config["Use global Z offset"] = {'label': 'Activate', 'datatype': 'BOOLEAN', 'def_value': 'False',
                                              'group': 'Global Z offset', 'section': 'Global Z offset'}

        self.config["Extend Z offset to readout"] = {'label': 'Extend Z offset to readout', 'datatype': 'BOOLEAN',
                                                     'def_value': 'True', 'state_quant': 'Use global Z offset',
                                                     'state_value': 'True', 'group': 'Global Z offset',
                                                     'section': 'Global Z offset'}

        self.config["Time after readout, Z global"] = {'label': 'Time after readout', 'datatype': 'DOUBLE', 'unit': 's',
                                                       'def_value': '10E-9',
                                                       'state_quant': 'Extend Z offset to readout',
                                                       'state_value': 'True', 'group': 'Global Z offset',
                                                       'section': 'Global Z offset'}

        self.config["Ringup, Z global"] = {'label': 'Ring-up time', 'datatype': 'DOUBLE', 'unit': 's',
                                           'def_value': '20E-9', 'state_quant': 'Use global Z offset',
                                           'state_value': 'True', 'group': 'Cosine pulse settings',
                                           'section': 'Global Z offset'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Amplitude #{}, Z global".format(qubit)] = {'label': 'Amplitude #{}'.format(qubit),
                                                                    'datatype': 'DOUBLE', 'unit': 'V',
                                                                    'def_value': '0.',
                                                                    'state_quant': 'Use global Z offset',
                                                                    'state_value': 'True', 'low_lim': '-1.5',
                                                                    'high_lim': '1.5', 'group': 'Z offset amplitudes',
                                                                    'section': 'Global Z offset'}

    def _z_pulse_during_readout(self):

        self.config["Use Z pulse during readout"] = {'label': 'Activate', 'datatype': 'BOOLEAN', 'def_value': 'False',
                                                     'group': 'Z pulse during readout',
                                                     'section': 'Z pulse during readout'}

        self.config["Net zero, Z during readout"] = {'label': 'Make Net Zero', 'datatype': 'BOOLEAN',
                                                     'def_value': 'False', 'group': 'Z pulse during readout',
                                                     'section': 'Z pulse during readout'}

        self.config["Ringup time, Z during readout"] = {'label': 'Ring-up time', 'datatype': 'DOUBLE', 'unit': 's',
                                                        'def_value': '20E-9',
                                                        'state_quant': 'Use Z pulse during readout',
                                                        'state_value': 'True', 'group': 'Cosine pulse settings',
                                                        'section': 'Z pulse during readout'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Amplitude #{}, Z during readout".format(qubit)] = {'label': 'Amplitude #{}'.format(qubit),
                                                                            'datatype': 'DOUBLE', 'unit': 'V',
                                                                            'def_value': '0.',
                                                                            'state_quant': 'Use Z pulse during readout',
                                                                            'state_value': 'True', 'low_lim': '-1.5',
                                                                            'high_lim': '1.5',
                                                                            'group': 'Z amplitudes during readout',
                                                                            'section': 'Z pulse during readout'}

    def _qb_spectra(self):
        for qubit in range(1, self.n_qubits + 1):
            self.config["Ec #{}".format(qubit)] = {'datatype': 'DOUBLE', 'label': 'Ec', 'unit': 'Hz',
                                                   'def_value': '200E6', 'group': 'Qubit #{}'.format(qubit),
                                                   'section': 'QB spectra'}
            self.config["f01 max #{}".format(qubit)] = {'datatype': 'DOUBLE', 'label': 'f01 max', 'unit': 'Hz',
                                                        'def_value': '5E9', 'group': 'Qubit #{}'.format(qubit),
                                                        'section': 'QB spectra'}
            self.config["f01 min #{}".format(qubit)] = {'datatype': 'DOUBLE', 'label': 'f01 min', 'unit': 'Hz',
                                                        'def_value': '4E9', 'group': 'Qubit #{}'.format(qubit),
                                                        'section': 'QB spectra'}
            self.config["Vperiod #{}".format(qubit)] = {'datatype': 'DOUBLE', 'label': 'Voltage period', 'unit': 'V',
                                                        'def_value': '1', 'group': 'Qubit #{}'.format(qubit),
                                                        'section': 'QB spectra'}
            self.config["Voffset #{}".format(qubit)] = {'datatype': 'DOUBLE', 'label': 'Voltage offset', 'unit': 'V',
                                                        'def_value': '0', 'group': 'Qubit #{}'.format(qubit),
                                                        'section': 'QB spectra'}
            self.config["V0 #{}".format(qubit)] = {'datatype': 'DOUBLE', 'label': 'Voltage operating point',
                                                   'unit': 'V', 'def_value': '0', 'group': 'Qubit #{}'.format(qubit),
                                                   'section': 'QB spectra'}

    def _two_qb_gates(self):
        self.config["Pulse type, 2QB"] = {'label': 'Pulse type',
                                          'datatype': 'COMBO',
                                          'combo_def_1': 'Gaussian',
                                          'combo_def_2': 'Square',
                                          'combo_def_3': 'Ramp',
                                          'combo_def_4': 'CZ',
                                          'combo_def_5': 'Cosine',
                                          'combo_def_6': 'NetZero',
                                          'def_value': 'Cosine',
                                          'group': '2-QB pulses',
                                          'section': '2-QB gates'}

        self.config["Truncation range, 2QB"] = {'label': 'Truncation range',
                                                'datatype': 'DOUBLE',
                                                'def_value': '3',
                                                'state_quant': 'Pulse type, 2QB',
                                                'state_value': 'Gaussian',
                                                'group': '2-QB pulses',
                                                'section': '2-QB gates'}

        self.config["Start at zero, 2QB"] = {'label': 'Start at zero',
                                             'datatype': 'BOOLEAN',
                                             'def_value': '0',
                                             'state_quant': 'Pulse type, 2QB',
                                             'state_value': 'Gaussian',
                                             'group': '2-QB pulses',
                                             'section': '2-QB gates'}

        self.config["Uniform 2QB pulses"] = {'datatype': 'BOOLEAN', 'def_value': 'True', 'group': '2-QB pulses',
                                             'section': '2-QB gates'}

        self.config["Width, 2QB"] = {'label': 'Width', 'datatype': 'DOUBLE', 'unit': 's', 'def_value': '50E-9',
                                     'group': '2-QB pulses', 'section': '2-QB gates',
                                     'state_quant': 'Uniform 2QB pulses', 'state_value': '1',
                                     'show_in_measurement_dlg': 'True'}

        self.config["Plateau, 2QB"] = {'label': 'Plateau',
                                       'datatype': 'DOUBLE',
                                       'unit': 's',
                                       'group': '2-QB pulses',
                                       'section': '2-QB gates',
                                       'state_quant': 'Uniform 2QB pulses',
                                       'state_value': '1',
                                       'show_in_measurement_dlg': 'True'}

        self.config["Fourier terms, 2QB"] = {'label': 'Fourier Terms',
                                             'datatype': 'COMBO',
                                             'combo_def_1': 'One',
                                             'combo_def_2': 'Two',
                                             'combo_def_3': 'Three',
                                             'combo_def_4': 'Four',
                                             'def_value': 'One',
                                             'tooltip': 'Number of fourier terms used to define the pulseshape of the C-phase gate',
                                             'group': '2-QB pulses',
                                             'section': '2-QB gates',
                                             'state_quant': 'Pulse type, 2QB',
                                             'state_value_1': 'CZ',
                                             'state_value_2': 'NetZero'}

        for qubit_1, qubit_2 in self.connected_qubits:
            self.config["Amplitude, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'Amplitude',
                                                                              'datatype': 'DOUBLE', 'unit': 'V',
                                                                              'def_value': '1.0',
                                                                              'group': '2-QB pulse #{0}-{1}'.format(
                                                                                  qubit_1, qubit_2),
                                                                              'section': '2-QB gates',
                                                                              'state_quant': 'Pulse type, 2QB',
                                                                              'state_value_1': 'Gaussian',
                                                                              'state_value_2': 'Square',
                                                                              'state_value_3': 'Ramp',
                                                                              'state_value_4': 'Cosine',
                                                                              'show_in_measurement_dlg': 'True'}
            self.config["Width, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'Width', 'datatype': 'DOUBLE',
                                                                          'unit': 's', 'def_value': '50E-9',
                                                                          'group': '2-QB pulse #{0}-{1}'.format(qubit_1,
                                                                                                               qubit_2),
                                                                          'section': '2-QB gates',
                                                                          'state_quant': 'Uniform 2QB pulses',
                                                                          'state_value': '0'}
            self.config["Plateau, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'Plateau', 'datatype': 'DOUBLE',
                                                                            'unit': 's', 'def_value': '0',
                                                                            'group': '2-QB pulse #{0}-{1}'.format(
                                                                                qubit_1, qubit_2),
                                                                            'section': '2-QB gates',
                                                                            'state_quant': 'Uniform 2QB pulses',
                                                                            'state_value': '0'}
            self.config["Assume linear dependence, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {
                'label': 'Assume linear dependence', 'datatype': 'BOOLEAN', 'def_value': 'True',
                'tooltip': 'Assumes a linear dependence between frequency and voltage. If false, use the qubit '
                           'spectrum.',
                'group': '2-QB pulse #{0}-{1}'.format(qubit_1, qubit_2), 'section': '2-QB gates',
                'state_quant': 'Pulse type, 2QB', 'state_value_1': 'CZ', 'state_value_2': 'NetZero'}
            self.config["df/dV, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'df/dV', 'datatype': 'DOUBLE',
                                                                          'unit': 'Hz/V', 'def_value': '0.5E9',
                                                                          'tooltip': 'Translating pulseshape from frequency space to voltage assuming a linear dependence',
                                                                          'group': '2-QB pulse #{0}-{1}'.format(qubit_1,
                                                                                                               qubit_2),
                                                                          'section': '2-QB gates',
                                                                          'state_quant': 'Pulse type, 2QB',
                                                                          'state_value_1': 'CZ',
                                                                          'state_value_2': 'NetZero'}
            self.config["f11-f20 initial, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'f11-f20 initial',
                                                                                    'datatype': 'DOUBLE', 'unit': 'Hz',
                                                                                    'def_value': '0.3E9',
                                                                                    'tooltip': 'Initial frequency splitting between the |11> state and the |02> state',
                                                                                    'group': '2-QB pulse #{0}-{1}'.format(
                                                                                        qubit_1, qubit_2),
                                                                                    'section': '2-QB gates',
                                                                                    'state_quant': 'Pulse type, 2QB',
                                                                                    'state_value_1': 'CZ',
                                                                                    'state_value_2': 'NetZero'}
            self.config["f11-f20 final, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'f11-f20 final',
                                                                              'datatype': 'DOUBLE', 'unit': 'Hz',
                                                                              'def_value': '0.05E9',
                                                                              'tooltip': 'Smallest frequency splitting between the |11> state and the |02> state during the C-phase gate',
                                                                              'group': '2-QB pulse #{0}-{1}'.format(
                                                                                  qubit_1, qubit_2),
                                                                              'section': '2-QB gates',
                                                                              'state_quant': 'Pulse type, 2QB',
                                                                              'state_value_1': 'CZ',
                                                                              'state_value_2': 'NetZero'}
            self.config["Coupling, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'Coupling', 'datatype': 'DOUBLE',
                                                                             'unit': 'Hz', 'def_value': '0.025E9',
                                                                             'tooltip': 'Coupling strength between |11> state and |02> state',
                                                                             'group': '2-QB pulse #{0}-{1}'.format(
                                                                                 qubit_1, qubit_2),
                                                                             'section': '2-QB gates',
                                                                             'state_quant': 'Pulse type, 2QB',
                                                                             'state_value_1': 'CZ',
                                                                             'state_value_2': 'NetZero'}
            self.config["L1, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'L1', 'datatype': 'DOUBLE',
                                                                       'def_value': '1',
                                                                       'tooltip': 'First fourier coefficient used to define the pulse shape of the C-phase gate.',
                                                                       'group': '2-QB pulse #{0}-{1}'.format(qubit_1,
                                                                                                            qubit_2),
                                                                       'section': '2-QB gates',
                                                                       'state_quant': 'Fourier terms, 2QB',
                                                                       'state_value_1': 'One', 'state_value_2': 'Two',
                                                                       'state_value_3': 'Three',
                                                                       'state_value_4': 'Four'}
            self.config["L2, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'L2', 'datatype': 'DOUBLE',
                                                                       'def_value': '0',
                                                                       'group': '2-QB pulse #{0}-{1}'.format(qubit_1,
                                                                                                            qubit_2),
                                                                       'section': '2-QB gates',
                                                                       'state_quant': 'Fourier terms, 2QB',
                                                                       'state_value_1': 'Two', 'state_value_2': 'Three',
                                                                       'state_value_3': 'Four'}
            self.config["L3, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'L3', 'datatype': 'DOUBLE',
                                                                       'def_value': '0',
                                                                       'group': '2-QB pulse #{0}-{1}'.format(qubit_1,
                                                                                                            qubit_2),
                                                                       'section': '2-QB gates',
                                                                       'state_quant': 'Fourier terms, 2QB',
                                                                       'state_value_1': 'Three',
                                                                       'state_value_2': 'Four'}
            self.config["L4, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'L4', 'datatype': 'DOUBLE',
                                                                       'def_value': '0',
                                                                       'group': '2-QB pulse #{0}-{1}'.format(qubit_1,
                                                                                                            qubit_2),
                                                                       'section': '2-QB gates',
                                                                       'state_quant': 'Fourier terms, 2QB',
                                                                       'state_value_1': 'Four'}
            self.config["QB1 Phi, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'QB1 Phase shift (rad)',
                                                                            'datatype': 'DOUBLE', 'def_value': '0',
                                                                            'group': '2-QB pulse #{0}-{1}'.format(
                                                                                qubit_1, qubit_2),
                                                                            'section': '2-QB gates'}
            self.config["QB2 Phi, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'QB2 Phase shift (rad)',
                                                                            'datatype': 'DOUBLE', 'def_value': '0',
                                                                            'group': '2-QB pulse #{0}-{1}'.format(
                                                                                qubit_1, qubit_2),
                                                                            'section': '2-QB gates'}
            self.config["Negative amplitude, 2QB #{0}-{1}".format(qubit_1, qubit_2)] = {'label': 'Negative amplitude',
                                                                                       'datatype': 'BOOLEAN',
                                                                                       'def_value': 'False',
                                                                                       'tooltip': 'Flip the sign of the amplitude of the CZ pulse',
                                                                                       'group': '2-QB pulse #{0}-{1}'.format(
                                                                                           qubit_1, qubit_2),
                                                                                       'section': '2-QB gates',
                                                                                       'state_quant': 'Pulse type, 2QB',
                                                                                       'state_value_1': 'CZ',
                                                                                       'state_value_2': 'NetZero'}

    def _tomography(self):
        self.config["Generate process tomography prepulse"] = {'datatype': 'BOOLEAN', 'def_value': '0',
                                                               'group': 'Tomography', 'section': 'Tomography'}

        self.config["Generate state tomography postpulse"] = {'datatype': 'BOOLEAN', 'def_value': '0',
                                                              'group': 'Tomography', 'section': 'Tomography'}

        if self.n_qubits == 1:
            self.config["Tomography scheme"] = {'label': 'Tomography scheme', 'datatype': 'COMBO',
                                                'def_value': 'Single qubit', 'combo_def_1': 'Single qubit',
                                                'group': 'State tomography setup', 'section': 'Tomography'}

        else:
            self.config["Tomography scheme"] = {'label': 'Tomography scheme', 'datatype': 'COMBO',
                                            'def_value': 'Single qubit', 'combo_def_1': 'Single qubit',
                                            'combo_def_2': 'Two qubit (9 pulse set)',
                                            'combo_def_3': 'Two qubit (30 pulse set)',
                                            'combo_def_4': 'Two qubit (36 pulse set)',
                                            'group': 'State tomography setup', 'section': 'Tomography'}

            self.config["Process tomography prepulse index 2-QB"] = {'label': 'Prepulse index 2-QB',
                                                                     'datatype': 'COMBO',
                                                                     'def_value': '00: I-I', 'combo_def_1': '00: I-I',
                                                                     'combo_def_2': '01: I-Xp',
                                                                     'combo_def_3': '0X: I-Y2p',
                                                                     'combo_def_4': '0Y: I-X2m',
                                                                     'combo_def_5': '10: Xp-I',
                                                                     'combo_def_6': '11: Xp-Xp',
                                                                     'combo_def_7': '1X: Xp-Y2p',
                                                                     'combo_def_8': '1Y: Xp-X2m',
                                                                     'combo_def_9': 'X0: Y2p-I',
                                                                     'combo_def_10': 'X1: Y2p-Xp',
                                                                     'combo_def_11': 'XX: Y2p-Y2p',
                                                                     'combo_def_12': 'XY: Y2p-X2m',
                                                                     'combo_def_13': 'Y0: X2m-I',
                                                                     'combo_def_14': 'Y1: X2m-Xp',
                                                                     'combo_def_15': 'YX: X2m-Y2p',
                                                                     'combo_def_16': 'YY: X2m-X2m',
                                                                     'group': 'State tomography setup',
                                                                     'section': 'Tomography',
                                                                     'show_in_measurement_dlg': 'True',
                                                                     'state_quant': 'Tomography scheme',
                                                                     'state_value_1': 'Two qubit (9 pulse set)',
                                                                     'state_value_2': 'Two qubit (30 pulse set)',
                                                                     'state_value_3': 'Two qubit (36 pulse set)'}

        self.config["Process tomography prepulse index 1-QB"] = {'label': 'Prepulse index 1-QB', 'datatype': 'COMBO',
                                                                 'def_value': '0: I', 'combo_def_1': '0: I',
                                                                 'combo_def_2': '1: Xp', 'combo_def_3': 'X: Y2p',
                                                                 'combo_def_4': 'Y: X2m',
                                                                 'group': 'State tomography setup',
                                                                 'section': 'Tomography',
                                                                 'show_in_measurement_dlg': 'True',
                                                                 'state_quant': 'Tomography scheme',
                                                                 'state_value': 'Single qubit'}



        self.config["Tomography pulse index 1-QB"] = {'label': 'Postpulse index 1-QB', 'datatype': 'COMBO',
                                                      'def_value': 'Z: I', 'combo_def_1': 'Z: I',
                                                      'combo_def_2': 'Y: X2p', 'combo_def_3': 'X: Y2m',
                                                      'group': 'State tomography setup', 'section': 'Tomography',
                                                      'show_in_measurement_dlg': 'True',
                                                      'state_quant': 'Tomography scheme',
                                                      'state_value_1': 'Single qubit'}
        if self.n_qubits > 1:
            self.config["Tomography pulse index 2-QB (9 pulse set)"] = {'label': 'Postpulse index 2-QB',
                                                                    'datatype': 'COMBO', 'def_value': 'XX: Y2m-Y2m',
                                                                    'combo_def_1': 'XX: Y2m-Y2m',
                                                                    'combo_def_2': 'YX: X2p-Y2m',
                                                                    'combo_def_3': 'ZX: I-Y2m',
                                                                    'combo_def_4': 'XY: Y2m-X2p',
                                                                    'combo_def_5': 'YY: X2p-X2p',
                                                                    'combo_def_6': 'ZY: I-X2p',
                                                                    'combo_def_7': 'XZ: Y2m-I',
                                                                    'combo_def_8': 'YZ: X2p-I',
                                                                    'combo_def_9': 'ZZ: I-I',
                                                                    'group': 'State tomography setup',
                                                                    'section': 'Tomography',
                                                                    'show_in_measurement_dlg': 'True',
                                                                    'state_quant': 'Tomography scheme',
                                                                    'state_value': 'Two qubit (9 pulse set)'}

            self.config["Tomography pulse index 2-QB (30 pulse set)"] = {'label': 'Postpulse index 2-QB',
                                                                     'datatype': 'COMBO', 'def_value': 'I-I',
                                                                     'combo_def_1': 'I-I', 'combo_def_2': 'Xp-I',
                                                                     'combo_def_3': 'I-Xp', 'combo_def_4': 'X2p-I',
                                                                     'combo_def_5': 'X2p-X2p', 'combo_def_6': 'X2p-Y2p',
                                                                     'combo_def_7': 'X2p-Xp', 'combo_def_8': 'Y2p-I',
                                                                     'combo_def_9': 'Y2p-X2p',
                                                                     'combo_def_10': 'Y2p-Y2p',
                                                                     'combo_def_11': 'Y2p-Xp', 'combo_def_12': 'I-X2p',
                                                                     'combo_def_13': 'Xp-X2p', 'combo_def_14': 'I-Y2p',
                                                                     'combo_def_15': 'Xp-Y2p', 'combo_def_16': 'I-I',
                                                                     'combo_def_17': 'Xm-I', 'combo_def_18': 'I-Xm',
                                                                     'combo_def_19': 'X2m-I', 'combo_def_20': 'X2m-X2m',
                                                                     'combo_def_21': 'X2m-Y2m',
                                                                     'combo_def_22': 'X2m-Xm', 'combo_def_23': 'Y2m-I',
                                                                     'combo_def_24': 'Y2m-X2m',
                                                                     'combo_def_25': 'Y2m-Y2m',
                                                                     'combo_def_26': 'Y2m-Xm', 'combo_def_27': 'I-X2m',
                                                                     'combo_def_28': 'Xm-X2m', 'combo_def_29': 'I-Y2m',
                                                                     'combo_def_30': 'Xm-Y2m',
                                                                     'group': 'State tomography setup',
                                                                     'section': 'Tomography',
                                                                     'show_in_measurement_dlg': 'True',
                                                                     'state_quant': 'Tomography scheme',
                                                                     'state_value': 'Two qubit (30 pulse set)'}

            self.config["Tomography pulse index 2-QB (36 pulse set)"] = {'label': 'Postpulse index 2-QB',
                                                                     'datatype': 'COMBO', 'def_value': 'I-I',
                                                                     'combo_def_1': 'I-I', 'combo_def_2': 'Xp-I',
                                                                     'combo_def_3': 'X2p-I', 'combo_def_4': 'X2m-I',
                                                                     'combo_def_5': 'Y2p-I', 'combo_def_6': 'Y2m-I',
                                                                     'combo_def_7': 'Id-Xp', 'combo_def_8': 'Xp-Xp',
                                                                     'combo_def_9': 'X2p-Xp', 'combo_def_10': 'X2m-Xp',
                                                                     'combo_def_11': 'Y2p-Xp', 'combo_def_12': 'Y2m-Xp',
                                                                     'combo_def_13': 'I-X2p', 'combo_def_14': 'Xp-X2p',
                                                                     'combo_def_15': 'X2p-X2p',
                                                                     'combo_def_16': 'X2m-X2p',
                                                                     'combo_def_17': 'Y2p-Y2p',
                                                                     'combo_def_18': 'Y2m-Y2p', 'combo_def_19': 'I-X2m',
                                                                     'combo_def_20': 'Xp-X2m',
                                                                     'combo_def_21': 'X2p-X2m',
                                                                     'combo_def_22': 'X2m-X2m',
                                                                     'combo_def_23': 'Y2p-X2m',
                                                                     'combo_def_24': 'Y2m-X2m', 'combo_def_25': 'I-Y2p',
                                                                     'combo_def_26': 'Xp-Y2p',
                                                                     'combo_def_27': 'X2p-Y2p',
                                                                     'combo_def_28': 'X2m-Y2p',
                                                                     'combo_def_29': 'Y2p-Y2p',
                                                                     'combo_def_30': 'Y2m-Y2p', 'combo_def_31': 'I-Y2m',
                                                                     'combo_def_32': 'Xp-Y2m',
                                                                     'combo_def_33': 'X2p-Y2m',
                                                                     'combo_def_34': 'X2m-Y2m',
                                                                     'combo_def_35': 'Y2p-Y2m',
                                                                     'combo_def_36': 'Y2m-Y2m',
                                                                     'group': 'State tomography setup',
                                                                     'section': 'Tomography',
                                                                     'show_in_measurement_dlg': 'True',
                                                                     'state_quant': 'Tomography scheme',
                                                                     'state_value': 'Two qubit (36 pulse set)'}

        self.config["Qubit for tomography"] = {'label': 'Qubit #', 'datatype': 'COMBO', 'def_value': '1',
                                               'state_quant': 'Tomography scheme', 'state_value': 'Single qubit',
                                               'group': 'Qubits for tomography', 'section': 'Tomography'}

        for qubit in range(1, self.n_qubits + 1):
            self.config.set("Qubit for tomography", "combo_def_{}".format(qubit), "{}".format(qubit))



        if self.n_qubits > 1:
            self.config["Qubit 1 for tomography"] = {'label': 'Qubit 1 waveform #', 'datatype': 'COMBO',
                                                     'def_value': '1',
                                                     'state_quant': 'Tomography scheme',
                                                     'state_value_1': 'Two qubit (9 pulse set)',
                                                     'state_value_2': 'Two qubit (30 pulse set)',
                                                     'state_value_3': 'Two qubit (36 pulse set)',
                                                     'group': 'Qubits for tomography', 'section': 'Tomography'}

            self.config["Qubit 2 for tomography"] = {'label': 'Qubit 2 waveform #', 'datatype': 'COMBO', 'def_value': '2',
                                                 'state_quant': 'Tomography scheme',
                                                 'state_value_1': 'Two qubit (9 pulse set)',
                                                 'state_value_2': 'Two qubit (30 pulse set)',
                                                 'state_value_3': 'Two qubit (36 pulse set)',
                                                 'group': 'Qubits for tomography', 'section': 'Tomography'}

            for qubit in range(1, self.n_qubits + 1):
                self.config.set("Qubit 1 for tomography", "combo_def_{}".format(qubit), "{}".format(qubit))
                self.config.set("Qubit 2 for tomography", "combo_def_{}".format(qubit), "{}".format(qubit))

    def _predistortion(self):
        self.config["Predistort waveforms"] = {'datatype': 'BOOLEAN', 'def_value': '0', 'group': 'XY Predistortion',
                                               'section': 'Predistortion'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Transfer function #{}".format(qubit)] = {'datatype': 'PATH', 'group': 'XY Predistortion',
                                                                  'section': 'Predistortion'}

        self.config["Predistort Z"] = {'datatype': 'BOOLEAN', 'def_value': '0', 'group': 'Z predistortion',
                                       'section': 'Predistortion'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Predistort Z{} - A1".format(qubit)] = {'label': 'A1', 'datatype': 'DOUBLE',
                                                                'group': 'Z{}'.format(qubit),
                                                                'section': 'Predistortion'}
            self.config["Predistort Z{} - tau1".format(qubit)] = {'label': 'tau1', 'datatype': 'DOUBLE',
                                                                  'group': 'Z{}'.format(qubit),
                                                                  'section': 'Predistortion'}
            self.config["Predistort Z{} - A2".format(qubit)] = {'label': 'A2', 'datatype': 'DOUBLE',
                                                                'group': 'Z{}'.format(qubit),
                                                                'section': 'Predistortion'}
            self.config["Predistort Z{} - tau2".format(qubit)] = {'label': 'tau2', 'datatype': 'DOUBLE',
                                                                  'group': 'Z{}'.format(qubit),
                                                                  'section': 'Predistortion'}
            self.config["Predistort Z{} - A3".format(qubit)] = {'label': 'A3', 'datatype': 'DOUBLE',
                                                                'group': 'Z{}'.format(qubit),
                                                                'section': 'Predistortion'}
            self.config["Predistort Z{} - tau3".format(qubit)] = {'label': 'tau3', 'datatype': 'DOUBLE',
                                                                  'group': 'Z{}'.format(qubit),
                                                                  'section': 'Predistortion'}
            self.config["Predistort Z{} - A4".format(qubit)] = {'label': 'A4', 'datatype': 'DOUBLE',
                                                                'group': 'Z{}'.format(qubit),
                                                                'section': 'Predistortion'}
            self.config["Predistort Z{} - tau4".format(qubit)] = {'label': 'tau4', 'datatype': 'DOUBLE',
                                                                  'group': 'Z{}'.format(qubit),
                                                                  'section': 'Predistortion'}

    def _crosstalk(self):
        self.config["Compensate cross-talk"] = {'datatype': 'BOOLEAN', 'def_value': '0', 'group': 'Cross-talk',
                                                'section': 'Cross-talk'}

        self.config["Cross-talk (CT) matrix"] = {'datatype': 'PATH', 'group': 'Cross-talk', 'section': 'Cross-talk'}

        self.config["1-1 QB <--> Crosstalk matrix"] = {
            'tooltip': 'One-to-one QB to Cross-talk matrix element correspondence', 'datatype': 'BOOLEAN',
            'def_value': 'True', 'group': 'Cross-talk', 'section': 'Cross-talk'}

        for qubit_index in range(1, self.n_qubits + 1):
            self.config["CT-matrix element #{}".format(qubit_index)] = {
                'tooltip': 'Which QB/output corresponds to which element of the cross-talk matrix', 'datatype': 'COMBO',
                'def_value': 'None', 'state_quant': '1-1 QB <--> Crosstalk matrix', 'state_value': 'False',
                'group': 'Cross-talk', 'section': 'Cross-talk'}
            for target_qubit in range(1, self.n_qubits + 1):
                self.config.set("CT-matrix element #{}".format(qubit_index), "combo_def_{}".format(target_qubit),
                                "{}".format(target_qubit))

            self.config.set("CT-matrix element #{}".format(qubit_index), "combo_def_{}".format(self.n_qubits + 1),
                            "None")

    def _readout(self):
        self.config["Number of readout waveforms"] = {'label': 'Number of waveforms', 'datatype': 'COMBO',
                                                      'def_value': '1', 'group': 'Readout', 'section': 'Readout'}
        for qubit_index in range(1, self.n_qubits + 1):
            self.config.set("Number of readout waveforms", "combo_def_{}".format(qubit_index),
                            "{}".format(qubit_index))

        self.config["Readout pulse type"] = {'datatype': 'COMBO', 'combo_def_1': 'Gaussian', 'combo_def_2': 'Square',
                                             'combo_def_3': 'Ramp', 'combo_def_4': 'Cosine', 'def_value': 'Square',
                                             'group': 'Readout', 'section': 'Readout'}

        self.config["Readout truncation range"] = {'datatype': 'DOUBLE', 'def_value': '3',
                                                   'state_quant': 'Readout pulse type', 'state_value_1': 'Gaussian',
                                                   'group': 'Readout', 'section': 'Readout'}

        self.config["Readout start at zero"] = {'label': 'Start at zero', 'datatype': 'BOOLEAN', 'def_value': '0',
                                                'state_quant': 'Readout pulse type', 'state_value_1': 'Gaussian',
                                                'group': 'Readout', 'section': 'Readout'}

        self.config["Uniform readout amplitude"] = {'label': 'Uniform amplitude', 'datatype': 'BOOLEAN',
                                                    'def_value': 'True', 'group': 'Readout', 'section': 'Readout'}

        self.config["Readout amplitude"] = {'label': 'Amplitude', 'datatype': 'DOUBLE', 'def_value': '0.1', 'unit': 'V',
                                            'state_quant': 'Uniform readout amplitude', 'state_value': 'True',
                                            'group': 'Readout', 'section': 'Readout'}

        self.config["Uniform readout pulse shape"] = {'label': 'Uniform pulse shape', 'datatype': 'BOOLEAN',
                                                      'def_value': 'True', 'group': 'Readout', 'section': 'Readout'}

        self.config["Readout width"] = {'label': 'Width', 'datatype': 'DOUBLE', 'unit': 's', 'def_value': '10E-9',
                                        'state_quant': 'Uniform readout pulse shape', 'state_value_1': 'True',
                                        'group': 'Readout', 'section': 'Readout', 'show_in_measurement_dlg': 'True'}

        self.config["Readout duration"] = {'label': 'Duration', 'datatype': 'DOUBLE', 'def_value': '2E-6', 'unit': 's',
                                           'state_quant': 'Uniform readout pulse shape', 'state_value_1': 'True',
                                           'group': 'Readout', 'section': 'Readout'}

        self.config["Match main sequence waveform size"] = {'datatype': 'BOOLEAN',
                                                            'tooltip': 'If checked, the readout waveform will have '
                                                                       'the same number of point as the main pulse '
                                                                       'waveforms',
                                                            'def_value': 'True', 'group': 'Readout', 'section': 'Readout'}

        self.config["Distribute readout phases"] = {'datatype': 'BOOLEAN',
                                                    'tooltip': 'If checked, the readout tone phases will be '
                                                               'distributed to avoid large peak-to-peak voltages',
                                                    'def_value': '0', 'group': 'Readout', 'section': 'Readout'}

        self.config["Readout delay"] = {'datatype': 'DOUBLE', 'unit': 's', 'group': 'Readout', 'section': 'Readout',
                                        'show_in_measurement_dlg': 'True'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Readout I/Q ratio {}".format(qubit)] = {'datatype': 'DOUBLE', 'def_value': '1.0',
                                                                 'tooltip': 'Ratio of I/Q voltages to compensate for I/Q mixer arm imbalance.',
                                                                 'state_quant': 'Number of readout waveforms',
                                                                 'group': 'Readout waveform {}'.format(qubit),
                                                                 'section': 'Readout'}
            for value in range(65 - qubit):
                self.config.set("Readout I/Q ratio {}".format(qubit), self.state_values[value][0],
                                self.state_values[value][1])

            self.config["Readout IQ skew {}".format(qubit)] = {'label': 'IQ skew {}'.format(qubit), 'unit': 'deg',
                                                               'tooltip': 'Readout IQ mixer phase skew',
                                                               'datatype': 'DOUBLE', 'def_value': '0',
                                                               'state_quant': 'Number of readout waveforms',
                                                               'group': 'Readout waveform {}'.format(qubit),
                                                               'section': 'Readout'}

            for value in range(65 - qubit):
                self.config.set("Readout IQ skew {}".format(qubit), self.state_values[value][0],
                                self.state_values[value][1])

            self.config["Readout offset - I{}".format(qubit)] = {'label': 'Offset I{}'.format(qubit),
                                                                 'tooltip': 'Readout mixer I offset',
                                                                 'datatype': 'DOUBLE', 'def_value': '0', 'unit': 'V',
                                                                 'state_quant': 'Number of readout waveforms',
                                                                 'group': 'Readout waveform {}'.format(qubit),
                                                                 'section': 'Readout'}

            for value in range(65 - qubit):
                self.config.set("Readout offset - I{}".format(qubit), self.state_values[value][0],
                                self.state_values[value][1])

            self.config["Readout offset - Q{}".format(qubit)] = {'label': 'Offset Q{}'.format(qubit),
                                                                 'tooltip': 'Readout mixer Q offset',
                                                                 'datatype': 'DOUBLE', 'def_value': '0', 'unit': 'V',
                                                                 'state_quant': 'Number of readout waveforms',
                                                                 'group': 'Readout waveform {}'.format(qubit),
                                                                 'section': 'Readout'}

            for value in range(65 - qubit):
                self.config.set("Readout offset - Q{}".format(qubit), self.state_values[value][0],
                                self.state_values[value][1])

        for qubit in range(1, self.n_qubits + 1):
            self.config["Readout amplitude #{}".format(qubit)] = {'label': 'Amplitude', 'datatype': 'DOUBLE',
                                                                  'def_value': '0.1', 'unit': 'V',
                                                                  'state_quant': 'Uniform readout amplitude',
                                                                  'state_value_1': 'False',
                                                                  'group': 'Readout #{}'.format(qubit),
                                                                  'section': 'Readout'}
            self.config["Readout width #{}".format(qubit)] = {'label': 'Width', 'datatype': 'DOUBLE', 'unit': 's',
                                                              'def_value': '10E-9',
                                                              'state_quant': 'Uniform readout pulse shape',
                                                              'state_value_1': 'False',
                                                              'group': 'Readout #{}'.format(qubit),
                                                              'section': 'Readout', 'show_in_measurement_dlg': 'True'}
            self.config["Readout duration #{}".format(qubit)] = {'label': 'Duration', 'datatype': 'DOUBLE',
                                                                 'def_value': '2E-6', 'unit': 's',
                                                                 'state_quant': 'Uniform readout pulse shape',
                                                                 'state_value_1': 'False',
                                                                 'group': 'Readout #{}'.format(qubit),
                                                                 'section': 'Readout'}
            self.config["Readout frequency #{}".format(qubit)] = {'label': 'Frequency', 'datatype': 'DOUBLE',
                                                                  'def_value': '0.0', 'unit': 'Hz',
                                                                  'group': 'Readout #{}'.format(qubit),
                                                                  'section': 'Readout'}
            self.config["Readout target #{}".format(qubit)] = {'label': 'Waveform target', 'datatype': 'COMBO',
                                                               'def_value': '1',
                                                               'group': 'Readout #{}'.format(qubit),
                                                               'section': 'Readout'}

            for target in range(1, self.n_qubits + 1):
                self.config.set("Readout target #{}".format(qubit), 'combo_def_{}'.format(target), '{}'.format(target))

        self.config["Predistort readout waveform"] = {'datatype': 'BOOLEAN',
                                                      'tooltip': 'If checked, the readout waveform will be predistorted to increase the rise time',
                                                      'def_value': '0', 'group': 'Readout predistortion',
                                                      'section': 'Readout'}

        self.config["Resonator linewidth"] = {'datatype': 'DOUBLE', 'def_value': '1E6', 'unit': 'Hz',
                                              'tooltip': 'Measured resonator linewidth',
                                              'state_quant': 'Predistort readout waveform', 'state_value': '1',
                                              'group': 'Readout predistortion', 'section': 'Readout'}

        self.config["Target rise time"] = {'datatype': 'DOUBLE', 'def_value': '50E-9', 'unit': 's',
                                           'tooltip': 'Intended resonator ring-up time',
                                           'state_quant': 'Predistort readout waveform', 'state_value': '1',
                                           'group': 'Readout predistortion', 'section': 'Readout'}

        self.config["Generate readout trig"] = {'datatype': 'BOOLEAN', 'group': 'Readout trig', 'section': 'Readout'}

        self.config["Readout trig amplitude"] = {'datatype': 'DOUBLE', 'def_value': '1.0', 'unit': 'V',
                                                 'group': 'Readout trig', 'section': 'Readout'}

        self.config["Readout trig duration"] = {'datatype': 'DOUBLE', 'def_value': '20E-9', 'unit': 's',
                                                'group': 'Readout trig', 'section': 'Readout'}

    def _output(self):
        self.config["Swap IQ"] = {'datatype': 'BOOLEAN', 'group': 'Output', 'section': 'Output'}

        self.config["Generate gate"] = {'datatype': 'BOOLEAN', 'def_value': '1', 'group': 'Microwave gate switch',
                                        'section': 'Output'}

        self.config["Uniform gate"] = {'datatype': 'BOOLEAN', 'def_value': '0', 'state_quant': 'Generate gate',
                                       'state_value': '1', 'group': 'Microwave gate switch', 'section': 'Output'}

        self.config["External gate"] = {'datatype': 'BOOLEAN', 'def_value': '0',
                                        'state_quant': 'Generate gate', 'state_value': '1',
                                        'second_state_quant': 'Generate waveform primitives',
                                        'second_state_value': '1',
                                        'group': 'Microwave gate switch', 'section': 'Output',
                                        'tooltip': 'Gate and control are generated by different hardware. Only applicable when using primitives.'}

        self.config["Gate delay"] = {'datatype': 'DOUBLE', 'unit': 's', 'def_value': '-60E-9', 'low_lim': '-1E-6',
                                     'high_lim': '1E-6', 'group': 'Microwave gate switch', 'section': 'Output'}

        self.config["Gate overlap"] = {'datatype': 'DOUBLE', 'unit': 's', 'def_value': '20E-9',
                                       'group': 'Microwave gate switch', 'section': 'Output'}

        self.config["Minimal gate time"] = {'datatype': 'DOUBLE', 'unit': 's', 'def_value': '20E-9',
                                            'group': 'Microwave gate switch', 'section': 'Output'}

        self.config["Filter gate waveforms"] = {'datatype': 'BOOLEAN', 'def_value': 'False', 'group': 'Output filters',
                                                'section': 'Output'}

        self.config["Gate filter"] = {'datatype': 'COMBO', 'def_value': 'Hanning', 'combo_def_1': 'Rectangular',
                                      'combo_def_2': 'Bartlett', 'combo_def_3': 'Blackman', 'combo_def_4': 'Hamming',
                                      'combo_def_5': 'Hanning', 'combo_def_6': 'Kaiser',
                                      'state_quant': 'Filter gate waveforms', 'state_value_1': '1',
                                      'group': 'Output filters', 'section': 'Output'}

        self.config["Gate - Filter size"] = {'label': 'Filter size', 'datatype': 'DOUBLE', 'def_value': '5',
                                             'low_lim': '1', 'state_quant': 'Filter gate waveforms',
                                             'state_value_1': '1', 'group': 'Output filters', 'section': 'Output'}

        self.config["Gate - Kaiser beta"] = {'label': 'Kaiser beta parameter', 'datatype': 'DOUBLE',
                                             'def_value': '14.0', 'state_quant': 'Gate filter',
                                             'state_value_1': 'Kaiser', 'group': 'Output filters', 'section': 'Output'}

        self.config["Filter Z waveforms"] = {'datatype': 'BOOLEAN', 'def_value': 'False', 'group': 'Output filters',
                                             'section': 'Output'}

        self.config["Z filter"] = {'datatype': 'COMBO', 'def_value': 'Hanning', 'combo_def_1': 'Rectangular',
                                   'combo_def_2': 'Bartlett', 'combo_def_3': 'Blackman', 'combo_def_4': 'Hamming',
                                   'combo_def_5': 'Hanning', 'combo_def_6': 'Kaiser',
                                   'state_quant': 'Filter Z waveforms', 'state_value_1': '1', 'group': 'Output filters',
                                   'section': 'Output'}

        self.config["Z - Filter size"] = {'label': 'Filter size', 'datatype': 'DOUBLE', 'def_value': '5',
                                          'low_lim': '1', 'state_quant': 'Filter Z waveforms', 'state_value_1': '1',
                                          'group': 'Output filters', 'section': 'Output'}

        self.config["Z - Kaiser beta"] = {'label': 'Kaiser beta parameter', 'datatype': 'DOUBLE', 'def_value': '14.0',
                                          'state_quant': 'Z filter', 'state_value_1': 'Kaiser',
                                          'group': 'Output filters', 'section': 'Output'}

        for qubit in range(1, self.n_qubits + 1):

            self.config["Trace - I{}".format(qubit)] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's',
                                                        'datatype': 'VECTOR', 'state_quant': 'Number of qubits',
                                                        'permission': 'READ', 'group': 'Traces', 'section': 'Output',
                                                        'show_in_measurement_dlg': 'True'}

            for value in range(65 - qubit):
                self.config.set("Trace - I{}".format(qubit),
                                'state_value_{}'.format(value + 1),
                                '{}'.format(64 - value))

            self.config["Trace - Q{}".format(qubit)] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's',
                                                        'datatype': 'VECTOR', 'state_quant': 'Number of qubits',
                                                        'permission': 'READ', 'group': 'Traces', 'section': 'Output',
                                                        'show_in_measurement_dlg': 'True'}

            for value in range(65 - qubit):
                self.config.set("Trace - Q{}".format(qubit),
                                'state_value_{}'.format(value + 1),
                                '{}'.format(64 - value))

            self.config["Trace - G{}".format(qubit)] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's',
                                                        'datatype': 'VECTOR', 'state_quant': 'Number of qubits',
                                                        'permission': 'READ', 'group': 'Traces', 'section': 'Output',
                                                        'show_in_measurement_dlg': 'True'}

            for value in range(65 - qubit):
                self.config.set("Trace - G{}".format(qubit),
                                'state_value_{}'.format(value + 1),
                                '{}'.format(64 - value))

            self.config["Trace - Z{}".format(qubit)] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's',
                                                        'datatype': 'VECTOR', 'state_quant': 'Number of qubits',
                                                        'permission': 'READ', 'group': 'Traces', 'section': 'Output',
                                                        'show_in_measurement_dlg': 'True'}

            for value in range(65 - qubit):
                self.config.set("Trace - Z{}".format(qubit),
                                'state_value_{}'.format(value + 1),
                                '{}'.format(64 - value))

        self.config["Trace - Readout trig"] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's', 'datatype': 'VECTOR',
                                               'permission': 'READ', 'group': 'Traces', 'section': 'Output',
                                               'show_in_measurement_dlg': 'True'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Trace - Readout I{}".format(qubit)] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's',
                                                                'datatype': 'VECTOR',
                                                                'permission': 'READ',
                                                                'state_quant': 'Number of readout waveforms',
                                                                'group': 'Traces', 'section': 'Output',
                                                                'show_in_measurement_dlg': 'True'}

            for value in range(65 - qubit):
                self.config.set("Trace - Readout I{}".format(qubit), self.state_values[value][0],
                                self.state_values[value][1])

            self.config["Trace - Readout Q{}".format(qubit)] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's',
                                                                'datatype': 'VECTOR',
                                                                'permission': 'READ',
                                                                'state_quant': 'Number of readout waveforms',
                                                                'group': 'Traces',
                                                                'section': 'Output',
                                                                'show_in_measurement_dlg': 'True'}

            for value in range(65 - qubit):
                self.config.set("Trace - Readout Q{}".format(qubit), self.state_values[value][0],
                                self.state_values[value][1])

        self.config["Waveform primitives"] = {'datatype': 'VECTOR', 'permission': 'READ', 'group': 'Traces',
                                              'section': 'Output', 'show_in_measurement_dlg': 'True',
                                              'state_quant': 'Generate waveform primitives', 'state_value': 'True'}

        self.config["Control sequence"] = {'datatype': 'VECTOR', 'permission': 'READ', 'group': 'Traces',
                                           'section': 'Output', 'show_in_measurement_dlg': 'True',
                                           'state_quant': 'Generate waveform primitives', 'state_value': 'True'}

    def _demodulation(self):
        self.config["Demodulation - Skip"] = {'label': 'Skip start', 'datatype': 'DOUBLE', 'def_value': '0.0',
                                              'unit': 's', 'section': '(Legacy) Demodulation', 'group': '(Legacy) Demodulation'}

        self.config["Demodulation - Length"] = {'label': 'Length', 'datatype': 'DOUBLE', 'def_value': '1E-6',
                                                'unit': 's', 'section': '(Legacy) Demodulation', 'group': '(Legacy) Demodulation'}

        self.config["Demodulation - Frequency offset"] = {'label': 'Frequency offset', 'datatype': 'DOUBLE',
                                                          'def_value': '0', 'unit': 'Hz',
                                                          'tooltip': 'Frequency difference between LO used for up- and downconversion, f_down - f_up',
                                                          'section': '(Legacy) Demodulation', 'group': '(Legacy) Demodulation'}

        self.config["Demodulation - Number of records"] = {'label': 'Number of records', 'datatype': 'DOUBLE',
                                                           'def_value': '1', 'low_lim': '1',
                                                           'tooltip': 'NB!  This will be removed in future versions, records should be set by input waveform',
                                                           'section': '(Legacy) Demodulation', 'group': '(Legacy) Demodulation'}

        self.config["Use phase reference signal"] = {'datatype': 'BOOLEAN', 'def_value': 'True',
                                                     'section': '(Legacy) Demodulation', 'group': '(Legacy) Demodulation'}

        self.config["Demodulation - IQ"] = {'datatype': 'BOOLEAN', 'def_value': 'False', 'section': '(Legacy) Demodulation',
                                            'group': '(Legacy) Demodulation'}

        self.config["Demodulation - Input"] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's', 'datatype': 'VECTOR',
                                               'permission': 'WRITE', 'state_quant': 'Demodulation - IQ',
                                               'state_value': 'False', 'section': '(Legacy) Demodulation',
                                               'group': '(Legacy) Demodulation', 'show_in_measurement_dlg': 'True'}

        self.config["Demodulation - Input I"] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's', 'datatype': 'VECTOR',
                                                 'permission': 'WRITE', 'state_quant': 'Demodulation - IQ',
                                                 'state_value': 'True', 'section': '(Legacy) Demodulation',
                                                 'group': '(Legacy) Demodulation', 'show_in_measurement_dlg': 'True'}

        self.config["Demodulation - Input Q"] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's', 'datatype': 'VECTOR',
                                                 'permission': 'WRITE', 'state_quant': 'Demodulation - IQ',
                                                 'state_value': 'True', 'section': '(Legacy) Demodulation',
                                                 'group': '(Legacy) Demodulation', 'show_in_measurement_dlg': 'True'}

        self.config["Demodulation - Reference"] = {'unit': 'V', 'x_name': 'Time', 'x_unit': 's', 'datatype': 'VECTOR',
                                                   'permission': 'WRITE', 'section': '(Legacy) Demodulation',
                                                   'group': '(Legacy) Demodulation', 'show_in_measurement_dlg': 'True'}

        for qubit in range(1, self.n_qubits + 1):
            self.config["Voltage, QB{}".format(qubit)] = {'unit': 'V', 'datatype': 'COMPLEX',
                                                          'state_quant': 'Number of qubits', 'permission': 'READ',
                                                          'group': '(Legacy) Demodulation', 'section': '(Legacy) Demodulation'}

            for value in range(65 - qubit):
                self.config.set("Voltage, QB{}".format(qubit),
                                'state_value_{}'.format(value + 1),
                                '{}'.format(64 - value))

            self.config["Single-shot, QB{}".format(qubit)] = {'unit': 'V', 'datatype': 'VECTOR_COMPLEX',
                                                              'state_quant': 'Number of qubits', 'permission': 'READ',
                                                              'group': '(Legacy) Demodulation', 'section': '(Legacy) Demodulation',
                                                              'show_in_measurement_dlg': 'True'}

            for value in range(65 - qubit):
                self.config.set("Single-shot, QB{}".format(qubit),
                                'state_value_{}'.format(value + 1),
                                '{}'.format(64 - value))


if __name__ == '__main__':
    """Command line usage"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_qubits',
                        type=int,
                        default=2,
                        help='Number of qubits in the MultiQubit_PulseGenerator')

    parser.add_argument('--gates_file',
                        type=str,

                        help='Filepath to json file containing 2Q gate pairs')

    parser.add_argument('--sequencer',
                        type=bool,
                        default=False,
                        help="Generate matching Keysight PXI Sequencer Driver")

    parser.add_argument('--n_chassis',
                        type=int,
                        default=1,
                        help='Specify number of chassis for sequencer. Only matters if sequencer is set to True')

    args = parser.parse_args()


    if args.n_qubits > 64:
        raise ValueError("More than 64 qubits not supported")

    if args.gates_file:
        with open(args.gates_file, 'r') as f:
            gates_file = json.load(f)

        mqpg = MQPGConfig(n_qubits=args.n_qubits, connected_qubits=gates_file['gates'])

    else:
        mqpg = MQPGConfig(n_qubits=args.n_qubits)
    mqpg.generate()
    mqpg.write('MultiQubit_PulseGenerator.ini')

    if args.sequencer:
        driver_name = 'Keysight_PXI_Sequencer'
        path_driver = ""
        driverPaths = ScriptsAndSettings.get_driver_paths_api()
        for driverPath in driverPaths:
            path_driver = os.path.join(driverPath, driver_name)
            if os.path.isdir(path_driver):
                log.info('Path to driver: {}'.format(path_driver))
                break

        if not os.path.exists(path_driver):
            raise Exception('Sequencer driver location cannot be found')

        sys.path.append(path_driver)
        from generate_sequencer import SequencerConfig
        seq = SequencerConfig(n_qubits=args.n_qubits, n_chassis=args.n_chassis)
        seq.generate_sequencer()
        seq.write_sequencer(path_driver + '/Keysight_PXI_Sequencer.ini')

