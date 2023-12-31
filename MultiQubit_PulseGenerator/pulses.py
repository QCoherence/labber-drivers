# Copyright 2014-2021 Keysight Technologies
#!/usr/bin/env python3
import numpy as np
import logging
import copy
from collections import OrderedDict

log = logging.getLogger('LabberDriver')

# TODO Private methods and variables


class Pulse:
    """Represents physical pulses played by an AWG.

    Parameters
    ----------
    complex : bool
        If True, pulse has both I and Q, otherwise it's real valued.
        Phase, frequency and drag only applies for complex waveforms.

    Attributes
    ----------
    amplitude : float
        Pulse amplitude.
    width : float
        Pulse width.
    plateau : float
        Pulse plateau.
    frequency : float
        SSB frequency.
    phase : float
        Pulse phase.
    use_drag : bool
        If True, applies DRAG correction.
    drag_coefficient : float
        Drag coefficient.
    drag_detuning : float
        Applies a frequnecy detuning for DRAG pulses.
    start_at_zero : bool
        If True, forces the pulse to start in 0.

    """

    def __init__(self, complex=True):

        # set variables
        self.amplitude = 0.5
        self.width = 10E-9
        self.plateau = 0.0
        self.frequency = 0.0
        self.phase = 0.0
        self.use_drag = False
        self.drag_coefficient = 0.0
        self.drag_detuning = 0.0
        self.start_at_zero = False
        self.complex = complex
        self.iq_skew = 0.0
        self.iq_ratio = 1.0

    def get_primitive_parameters(self, with_freq=False, with_phase=False, with_amp=False):
        """Get primitive parameters for pulse.
        """
        d = OrderedDict()
        d['width'] = self.width
        d['plateau'] = self.plateau
        d['use_drag'] = self.use_drag
        d['drag_coefficient'] = self.drag_coefficient
        d['drag_detuning'] = self.drag_detuning
        d['start_at_zero'] = self.start_at_zero
        d['complex'] = self.complex
        d['iq_skew'] = self.iq_skew
        d['iq_ratio'] = self.iq_ratio
        # do not include the following, should be set dynamically
        if with_amp:
            d['amplitude'] = self.amplitude
        if with_freq:
            d['frequency'] = self.frequency
        if with_phase:
            d['phase'] = self.phase
        return d

    def update_parameters_for_primitive(self, d):
        """Update parameters for generating primitives.
        """
        self.width = d.get('width')
        self.plateau = d.get('plateau')
        self.use_drag = d.get('use_drag')
        self.drag_coefficient = d.get('drag_coefficient')
        self.drag_detuning = d.get('drag_detuning')
        self.start_at_zero = d.get('start_at_zero')
        self.complex = d.get('complex')
        self.iq_skew = d.get('iq_skew')
        self.iq_ratio = d.get('iq_ratio')
        # dynamic quantities should be reset to default
        self.amplitude = d.get('amplitude', 1.0)
        self.frequency = d.get('frequency', 0.0)
        self.phase = d.get('phase', 0.0)

    def total_duration(self):
        """Get the total duration for the pulse.

        Returns
        -------
        float
            Total duration in seconds.

        """
        raise NotImplementedError()

    def calculate_envelope(self, t0, t):
        """Calculate pulse envelope.

        Parameters
        ----------
        t0 : float
            Pulse position, referenced to center of pulse.

        t : numpy array
            Array with time values for which to calculate the pulse envelope.

        Returns
        -------
        waveform : numpy array
            Array containing pulse envelope.

        """
        raise NotImplementedError()

    def calculate_primitive(self, sample_rate):
        """Calculate waveform for pulse primitive at given sample rate"""
        total_t = self.total_duration()
        n = int(round(total_t * sample_rate))
        # naturally make waveform end in zero
        # add buffer elements for drags/gradient
        t = np.arange(0, n + 2, dtype=float) / sample_rate
        y = self.calculate_waveform(total_t/2, t)
        # strip buffer elements
        y = y[1:-1]
        # also force last point to zero
        y[-1] = 0
        return y

    def calculate_waveform(self, t0, t, ignore_drag_modulation=False):
        """Calculate pulse waveform including phase shifts and SSB-mixing.

        Parameters
        ----------
        t0 : float
            Pulse position, referenced to center of pulse.

        t : numpy array
            Array with time values for which to calculate the pulse waveform.

        ignore_drag_modulation : bool
            If True, drag and modulation is disabled.

        Returns
        -------
        waveform : numpy array
            Array containing pulse waveform.

        """
        y = self.calculate_envelope(t0, t)
        # Make sure the waveform is zero outside the pulse
        y[t < (t0 - self.total_duration() / 2)] = 0
        y[t > (t0 + self.total_duration() / 2)] = 0

        # speed up generation by not applying transforms if not necessary
        if not self.complex:
            return y

        # if ignoring modulation, just apply global phase and return
        if ignore_drag_modulation:
            y = y * np.exp(- 1j * self.phase)
            return y

        if self.use_drag and self.complex:
            beta = self.drag_coefficient / (t[1] - t[0])
            y = y + 1j * beta * np.gradient(y)
            y = y * np.exp(1j * 2 * np.pi * self.drag_detuning *
                           (t - t0 + self.total_duration() / 2))

        if self.complex:
            # Apply phase and SSB
            phase = self.phase
            # single-sideband mixing, get frequency
            omega = 2 * np.pi * self.frequency
            # apply SSBM transform
            data_i = (y.real * np.cos(omega * t - phase) +
                      -y.imag * np.cos(omega * t - phase + +np.pi / 2))
            data_q = (y.real * np.sin(omega * t - phase + self.iq_skew) +
                      -y.imag * np.sin(omega * t - phase + np.pi / 2 +
                      self.iq_skew))
            y = self.iq_ratio * data_i + 1j * data_q
        return y


class Gaussian(Pulse):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.truncation_range = 5

    def get_primitive_parameters(self, *args, **kwargs):
        d = super().get_primitive_parameters(*args, **kwargs)
        d['truncation_range'] = self.truncation_range
        return d

    def update_parameters_for_primitive(self, d):
        """Update parameters for generating primitives.
        """
        super().update_parameters_for_primitive(d)
        self.truncation_range = d.get('truncation_range')

    def total_duration(self):
        return self.truncation_range * self.width + self.plateau

    def calculate_envelope(self, t0, t):
        # width is two t std
        # std = self.width/2;
        # alternate; std is set to give total pulse area same as a square
        std = self.width / np.sqrt(2 * np.pi)
        values = np.zeros_like(t)
        if self.plateau == 0:
            # pure gaussian, no plateau
            if std > 0:
                values = np.exp(-(t - t0)**2 / (2 * std**2))
        else:
            # add plateau
            values = np.array(
                ((t >= (t0 - self.plateau / 2)) & (t <
                                                   (t0 + self.plateau / 2))),
                dtype=float)
            if std > 0:
                # before plateau
                values += ((t < (t0 - self.plateau / 2)) * np.exp(
                    -(t - (t0 - self.plateau / 2))**2 / (2 * std**2)))
                # after plateau
                values += ((t >= (t0 + self.plateau / 2)) * np.exp(
                    -(t - (t0 + self.plateau / 2))**2 / (2 * std**2)))

        # TODO  Fix this
        if self.start_at_zero:
            values = values - values.min()
            values = values / values.max()
        values = values * self.amplitude

        return values


class Ramp(Pulse):
    def total_duration(self):
        return 2 * self.width + self.plateau

    def calculate_envelope(self, t0, t):
        # rising and falling slopes
        vRise = ((t - (t0 - self.plateau / 2 - self.width)) / self.width)
        vRise[vRise < 0.0] = 0.0
        vRise[vRise > 1.0] = 1.0
        vFall = (((t0 + self.plateau / 2 + self.width) - t) / self.width)
        vFall[vFall < 0.0] = 0.0
        vFall[vFall > 1.0] = 1.0
        values = vRise * vFall

        values = values * self.amplitude

        return values


class Square(Pulse):
    def total_duration(self):
        return self.width + self.plateau

    def calculate_envelope(self, t0, t):
        # reduce risk of rounding errors by putting checks between samples
        if len(t) > 1:
            t0 += (t[1] - t[0]) / 2.0

        values = ((t >= (t0 - (self.width + self.plateau) / 2)) &
                  (t < (t0 + (self.width + self.plateau) / 2)))

        values = values * self.amplitude

        return values


class Cosine(Pulse):
    def total_duration(self):
        return self.width + self.plateau

    def calculate_envelope(self, t0, t):
        tau = self.width
        if self.plateau == 0:
            values = (self.amplitude / 2 *
                      (1 - np.cos(2 * np.pi * (t - t0 + tau / 2) / tau)))
        else:
            values = np.ones_like(t) * self.amplitude
            values[t < t0 - self.plateau / 2] = self.amplitude / 2 * \
                (1 - np.cos(2 * np.pi *
                            (t[t < t0 - self.plateau / 2] - t0 +
                             self.plateau / 2 + tau / 2) / tau))
            values[t > t0 + self.plateau / 2] = self.amplitude / 2 * \
                (1 - np.cos(2 * np.pi *
                            (t[t > t0 + self.plateau / 2] - t0 -
                             self.plateau / 2 + tau / 2) / tau))

        return values


class CZ(Pulse):
    def __init__(self, *args, **kwargs):
        super().__init__(False)
        # For CZ pulses
        self.F_Terms = 1
        self.Coupling = 20E6
        self.Offset = 300E6
        self.Lcoeff = np.array([0.3])
        self.dfdV = 500E6
        self.qubit = None
        self.negative_amplitude = False

        self.t_tau = None

    def get_primitive_parameters(self, *args, **kwargs):
        d = super().get_primitive_parameters(*args, **kwargs)
        d['F_Terms'] = self.F_Terms
        d['Coupling'] = self.Coupling
        d['Offset'] = self.Offset
        d['Lcoeff'] = self.Lcoeff.tolist()
        d['dfdV'] = self.dfdV
        d['qubit'] = self.qubit
        d['negative_amplitude'] = self.negative_amplitude
        d['t_tau'] = self.t_tau
        return d

    def update_parameters_for_primitive(self, d):
        super().update_parameters_for_primitive(d)
        self.F_Terms = d.get('F_Terms', 1)
        self.Coupling = d.get('Coupling', 20E6)
        self.Offset = d.get('Offset', 300E6)
        self.Lcoeff = np.array(d.get('Lcoeff', [0.3]))
        self.dfdV = d.get('dfdV', 500E6)
        self.qubit = d.get('qubit', None)
        self.negative_amplitude = d.get('negative_amplitude', False)
        self.t_tau = d.get('t_tau', None)
        self.calculate_cz_waveform()

    def total_duration(self):
        return self.width+self.plateau

    def calculate_envelope(self, t0, t):
        if self.t_tau is None:
            self.calculate_cz_waveform()

        # Plateau is added as an extra extension of theta_f.
        theta_t = np.ones(len(t)) * self.theta_i
        for i in range(len(t)):
            if 0 < (t[i] - t0 + self.plateau / 2) < self.plateau:
                theta_t[i] = self.theta_f
            elif (0 < (t[i] - t0 + self.width / 2 + self.plateau / 2) <
                  (self.width + self.plateau) / 2):
                theta_t[i] = np.interp(
                    t[i] - t0 + self.width / 2 + self.plateau / 2, self.t_tau,
                    self.theta_tau)

            elif (0 < (t[i] - t0 + self.width / 2 + self.plateau / 2) <
                  (self.width + self.plateau)):
                theta_t[i] = np.interp(
                    t[i] - t0 + self.width / 2 - self.plateau / 2, self.t_tau,
                    self.theta_tau)
        # Clip theta_t to remove numerical outliers:
        theta_t = np.clip(theta_t, self.theta_i, None)

        # clip theta_f to remove numerical outliers
        theta_t = np.clip(theta_t, self.theta_i, None)
        df = 2*self.Coupling * (1 / np.tan(theta_t) - 1 / np.tan(self.theta_i))

        if self.qubit is None:
            # Use linear dependence if no qubit was given
            # log.info('---> df (linear): ' +str(df))
            values = df / self.dfdV
            # values = theta_t
        else:
            values = self.qubit.df_to_dV(df)
        if self.negative_amplitude is True:
            values = -values

        return values

    def calculate_cz_waveform(self):
        """Calculate waveform for c-phase and store in object"""
        # notation and calculations are based on
        # "Fast adiabatic qubit gates using only sigma_z control"
        # PRA 90, 022307 (2014)
        # Initial and final angles on the |11>-|02> bloch sphere
        self.theta_i = np.arctan(2*self.Coupling / self.Offset)
        self.theta_f = np.arctan(2*self.Coupling / self.amplitude)
        # log.log(msg="calc", level=30)

        # Renormalize fourier coefficients to initial and final angles
        # Consistent with both Martinis & Geller and DiCarlo 1903.02492
        Lcoeff = self.Lcoeff
        Lcoeff[0] = (((self.theta_f - self.theta_i) / 2)
                     - np.sum(self.Lcoeff[range(2, self.F_Terms, 2)]))

        # defining helper variabels
        n = np.arange(1, self.F_Terms + 1, 1)
        n_points = 1000  # Number of points in the numerical integration

        # Calculate pulse width in tau variable - See paper for details
        tau = np.linspace(0, 1, n_points)
        self.theta_tau = np.zeros(n_points)
        # This corresponds to the sum in Eq. (15) in Martinis & Geller
        for i in range(n_points):
            self.theta_tau[i] = (
                np.sum(Lcoeff * (1 - np.cos(2 * np.pi * n * tau[i]))) +
                self.theta_i)
        # Now calculate t_tau according to Eq. (20)
        t_tau = np.trapz(np.sin(self.theta_tau), x=tau)
        # log.info('t tau: ' + str(t_tau))
        # t_tau = np.sum(np.sin(self.theta_tau))*(tau[1] - tau[0])
        # Find the width in units of tau:
        Width_tau = self.width / t_tau

        # Calculating time as functions of tau
        # we normalize to width_tau (calculated above)
        tau = np.linspace(0, Width_tau, n_points)
        self.t_tau = np.zeros(n_points)
        self.t_tau2 = np.zeros(n_points)
        for i in range(n_points):
            if i > 0:
                self.t_tau[i] = np.trapz(
                    np.sin(self.theta_tau[0:i+1]), x=tau[0:i+1])
                # self.t_tau[i] = np.sum(np.sin(self.theta_tau[0:i+1]))*(tau[1]-tau[0])

class NetZero(CZ):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.slepian = None

    def total_duration(self):
        return 2*self.slepian.total_duration()

    def calculate_cz_waveform(self):
        self.slepian = CZ()
        self.slepian.__dict__ = copy.copy(self.__dict__)
        self.slepian.width /= 2
        self.slepian.plateau /= 2
        self.slepian.calculate_cz_waveform()

    def calculate_envelope(self, t0, t):
        return (self.slepian.calculate_envelope(t0-self.total_duration()/4, t) -
                self.slepian.calculate_envelope(t0+self.total_duration()/4, t))


if __name__ == '__main__':
    pass
