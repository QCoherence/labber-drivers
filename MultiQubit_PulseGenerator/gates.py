# Copyright 2014-2021 Keysight Technologies
#!/usr/bin/env python3
from copy import copy, deepcopy
import numpy as np
import logging
from sequence import Step
log = logging.getLogger('LabberDriver')

# TODO remove Step dep from CompositeGate


class BaseGate:
    """Base class for a qubit gate.

    """
    def __init__(self, **dynamic_kwargs):
        """
        dynamic_kwargs: dict, optional
            Gate parameters which can be hardware-dynamic
            Keywords should be of the form 'dynamic_{}'
            Values should be int or None
        """
        self.dynamic_kwargs = {}
        self.update_dynamic_kwargs(**dynamic_kwargs)

    def update_dynamic_kwargs(self, **dynamic_kwargs):
        for k, v in dynamic_kwargs.items():
            if not k.startswith('dynamic'):
                raise NameError(f'Unexpected dynamic gate parameter {k}')
        # don't overwrite existing params
        # because there were set in the sequence definition
        for name, value in dynamic_kwargs.items():
            if name in self.dynamic_kwargs:
                value_old = self.dynamic_kwargs[name]
                log.info(f'WARNING: Dynamic gate parameter {name}'
                         + f'not overwritten from {value_old} to {value}')
            else:
                self.dynamic_kwargs[name] = value

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        return pulse

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__str__()

    def copy(self):
        return deepcopy(self)


class OneQubitGate(BaseGate):
    def number_of_qubits(self):
        return 1


class TwoQubitGate(BaseGate):
    def number_of_qubits(self):
        return 2


class SingleQubitXYRotation(OneQubitGate):
    """Single qubit rotations around the XY axes.

    Angles defined as in https://en.wikipedia.org/wiki/Bloch_sphere.

    Parameters
    ----------
    phi : float
        Rotation axis.
    theta : float
        Rotation angle.

    """

    def __init__(self, phi, theta, name=None, **dynamic_kwargs):
        super().__init__(**dynamic_kwargs)
        self.phi = phi
        self.theta = theta
        self.name = name

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        pulse.phase = self.phi
        # pi pulse correspond to the full amplitude
        pulse.amplitude *= self.theta / np.pi
        return pulse

    def __str__(self):
        if self.name is None:
            return "XYPhi={:+.6f}theta={:+.6f}".format(self.phi, self.theta)
        else:
            return self.name

    def __eq__(self, other):
        threshold = 1e-10
        if not isinstance(other, SingleQubitXYRotation):
            return False
        if np.abs(self.phi - other.phi) > threshold:
            return False
        if np.abs(self.theta - other.theta) > threshold:
            return False
        return True


class SingleQubitZRotation(OneQubitGate):
    """Single qubit rotation around the Z axis.

    Parameters
    ----------
    theta : float
        Roation angle.

    """

    def __init__(self, theta, name=None, **dynamic_kwargs):
        super().__init__(**dynamic_kwargs)
        self.theta = theta
        self.name = name

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        # pi pulse correspond to the full amplitude
        pulse.amplitude *= self.theta / np.pi
        return pulse

    def __str__(self):
        if self.name is None:
            return "Ztheta={:+.2f}".format(self.theta)
        else:
            return self.name

    def __eq__(self, other):
        threshold = 1e-10
        if not isinstance(other, SingleQubitZRotation):
            return False
        if np.abs(self.theta - other.theta) > threshold:
            return False
        return True


class IdentityGate(OneQubitGate):
    """Identity gate.

    Does nothing to the qubit. The width can be specififed to
    implement a delay in the sequence. If no width is given, the identity gate
    inherits the width of the given pulse.

    Parameters
    ----------
    width : float
        Width of the I gate in seconds,
        None uses the XY width (the default is None).

    """

    def __init__(self, width=None):
        super().__init__()
        self.width = width

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        pulse.amplitude = 0
        pulse.use_drag = False  # Avoids bug
        if self.width is not None:
            pulse.width = 0
            pulse.plateau = self.width
        return pulse

    def __str__(self):
        return "I"


class VirtualZGate(OneQubitGate):
    """Virtual Z Gate."""

    def __init__(self, theta, name=None):
        super().__init__()
        self.theta = theta
        self.name = name

    def __eq__(self, other):
        threshold = 1e-10
        if not isinstance(other, VirtualZGate):
            return False
        if np.abs(self.theta - other.theta) > threshold:
            return False
        return True

    def __str__(self):
        if self.name is None:
            return "VZtheta={:+.2f}".format(self.theta)
        else:
            return self.name


class CPHASE(TwoQubitGate):
    """ CPHASE gate. """


class iSWAP_no_1qb_phases(TwoQubitGate):
    """ ISWAP gate. """


class ReadoutGate(OneQubitGate):
    """Readouts the qubit state."""
    def __init__(self, amplitude=None,
                 **dynamic_kwargs):
        OneQubitGate.__init__(self, **dynamic_kwargs)
        self.amplitude = amplitude

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        if self.amplitude is not None:
            pulse.amplitude *= self.amplitude
        return pulse

class HeraldingGate(OneQubitGate):
    """Heralding of the qubit state."""
    def __init__(self, amplitude=None,
                 **dynamic_kwargs):
        OneQubitGate.__init__(self, **dynamic_kwargs)
        self.amplitude = amplitude

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        if self.amplitude is not None:
            pulse.amplitude *= self.amplitude
        return pulse


class CustomGate(BaseGate):
    """A gate using a given :obj:`Pulse`.

    Parameters
    ----------
    pulse : :obj:`Pulse`
        The corresponding pulse.

    """

    def __init__(self, pulse, **dynamic_kwargs):
        super().__init__(**dynamic_kwargs)
        self.pulse = pulse


class RabiGate(SingleQubitXYRotation):
    """Creates the Rabi gate used in the spin-locking sequence.

    Parameters
    ----------
    amplitude : Amplitude of the pulse
    plateau : The duration of the pulse.
    phase : Phase of the Rabi gate. 0 corresponds to rotation around X axis.
    frequency: Drive frequency
    width:  Pulse rise/fall time
    use_drag:  Turn on/off drag
    drag_coefficient: DRAG scaling
    drag_detuning: DRAG detuning
    iq_skew: Phase delay between I/Q arms
    iq_ratio: Imbalance between I/Q amplitudes
    """

    def __init__(self, amplitude=None, plateau=None, phase=None,
                 frequency=None, width=None,
                 use_drag=None, drag_coefficient=None,
                 drag_detuning=None, iq_skew=None, iq_ratio=None,
                 **dynamic_kwargs):
        # inherit from BaseGate instead of super()
        # b/c SingleQubitXYRotation.__init__ has irrelevant parameters
        BaseGate.__init__(self, **dynamic_kwargs)
        self.amplitude = amplitude
        self.plateau = plateau
        self.phase = phase
        self.frequency = frequency
        self.width = width
        self.use_drag = use_drag
        self.drag_coefficient = drag_coefficient
        self.drag_detuning = drag_detuning
        self.iq_skew = iq_skew
        self.iq_ratio = iq_ratio

    def get_adjusted_pulse(self, pulse):
        pulse = copy(pulse)
        if self.amplitude is not None:
            pulse.amplitude = self.amplitude
        if self.plateau is not None:
            pulse.plateau = self.plateau
        if self.phase is not None:
            pulse.phase = self.phase
        if self.frequency is not None:
            pulse.frequency = self.frequency
        if self.width is not None:
            pulse.width = self.width
        if self.use_drag is not None:
            pulse.use_drag = self.use_drag
        if self.drag_coefficient is not None:
            pulse.drag_coefficient = self.drag_coefficient
        if self.drag_detuning is not None:
            pulse.drag_detuning = self.drag_detuning
        if self.iq_skew is not None:
            pulse.iq_skew = self.iq_skew
        if self.iq_ratio is not None:
            pulse.iq_ratio = self.iq_ratio
        return pulse


class CompositeGate:
    """Multiple gates in one object.
    Hardware-dynamic quantities should be attached to the component gates.

    Parameters
    ----------
    n_qubit : int
        Number of qubits involved in the composite gate.

    Attributes
    ----------
    sequence : list of :Step:
        Holds the gates involved.

    """

    def __init__(self, n_qubit, name=None):
        self.n_qubit = n_qubit
        self.sequence = []
        self.name = name

    def add_gate(self, gate, qubit=None):
        """Add a set of gates to the given qubit.

        For the qubits with no specified gate, an IdentityGate will be given.
        The length of the step is given by the longest pulse.

        Parameters
        ----------
        qubit : int or list of int
            The qubit(s) to add the gate(s) to.
        gate : :obj:`BaseGate` or list of :obj:`BaseGate`
            The gate(s) to add.
        """
        if qubit is None:
            if self.n_qubit == 1:
                qubit = 0
            else:
                qubit = [n for n in range(self.n_qubit)]

        step = Step()
        if isinstance(gate, list):
            if len(gate) == 1:
                raise ValueError(
                    "For single gates, don't provide gate as a list.")
            if not isinstance(qubit, list):
                raise ValueError(
                    """Please provide qubit indices as a list when adding more
                    than one gate.""")
            if len(gate) != len(qubit):
                raise ValueError(
                    "Length of gate list must equal length of qubit list.")

            for q, g in zip(qubit, gate):
                step.add_gate(q, g)

        else:
            if gate.number_of_qubits() > 1:
                if not isinstance(qubit, list):
                    raise ValueError(
                        """Please provide qubit list for gates with more than
                        one qubit.""")
            else:
                if not isinstance(qubit, int):
                    raise ValueError(
                        "For single gates, give qubit as int (not list).")
            step.add_gate(qubit, gate)

        self.sequence.append(step)

    def number_of_qubits(self):
        return self.n_qubit

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        if self.name is not None:
            return self.name
        else:
            super().__str__()

    def __repr__(self):
        return self.__str__()

    def copy(self):
        return deepcopy(self)


class iSWAP_with_1qb_phases(CompositeGate):
    """iSWAP gate followed by single qubit Z rotations.

    Parameters
    ----------
    phi1 : float
        Z rotation angle for qubit 1.
    phi2 : float
        Z rotation angle for qubit 2.
    """

    def __init__(self, phi1, phi2, **dynamic_kwargs):
        super().__init__(n_qubit=2)
        self.add_gate(iSWAP_no_1qb_phases(**dynamic_kwargs))
        self.add_gate([VirtualZGate(phi1), VirtualZGate(phi2)])

    def new_angles(self, phi1, phi2):
        """Update the angles of the single qubit rotations.

        Parameters
        ----------
        phi1 : float
            Z rotation angle for qubit 1.
        phi2 : float
            Z rotation angle for qubit 2.
        """
        self.__init__(phi1, phi2)

    def __str__(self):
        return "iSWAP"


class CPHASE_with_1qb_phases(CompositeGate):
    """CPHASE gate followed by single qubit Z rotations.

    Parameters
    ----------
    phi1 : float
        Z rotation angle for qubit 1.
    phi2 : float
        Z rotation angle for qubit 2.
    """

    def __init__(self, phi1, phi2, **dynamic_kwargs):
        super().__init__(n_qubit=2)
        self.add_gate(CPHASE(**dynamic_kwargs))
        self.add_gate([VirtualZGate(phi1), VirtualZGate(phi2)])

    def new_angles(self, phi1, phi2):
        """Update the angles of the single qubit rotations.

        Parameters
        ----------
        phi1 : float
            Z rotation angle for qubit 1.
        phi2 : float
            Z rotation angle for qubit 2.
        """
        self.__init__(phi1, phi2)

    def __str__(self):
        return "CZ"


I = IdentityGate(width=None)
I0 = IdentityGate(width=0)
Ilong = IdentityGate(width=75e-9)

# X gates
Xp = SingleQubitXYRotation(phi=0, theta=np.pi, name='Xp')
Xm = SingleQubitXYRotation(phi=0, theta=-np.pi, name='Xm')
X2p = SingleQubitXYRotation(phi=0, theta=np.pi / 2, name='X2p')
X2m = SingleQubitXYRotation(phi=0, theta=-np.pi / 2, name='X2m')

# Y gates
Yp = SingleQubitXYRotation(phi=np.pi / 2, theta=np.pi, name='Yp')
Ym = SingleQubitXYRotation(phi=np.pi / 2, theta=-np.pi, name='Ym')
Y2m = SingleQubitXYRotation(phi=np.pi / 2, theta=-np.pi / 2, name='Y2m')
Y2p = SingleQubitXYRotation(phi=np.pi / 2, theta=np.pi / 2, name='Y2p')

# Z gates
Zp = SingleQubitZRotation(np.pi, name='Zp')
Z2p = SingleQubitZRotation(np.pi / 2, name='Z2p')
Zm = SingleQubitZRotation(-np.pi, name='Zm')
Z2m = SingleQubitZRotation(-np.pi / 2, name='Z2m')

# Virtual Z gates
VZp = VirtualZGate(np.pi, name='VZp')
VZ2p = VirtualZGate(np.pi / 2, name='VZ2p')
VZm = VirtualZGate(-np.pi, name='VZm')
VZ2m = VirtualZGate(-np.pi / 2, name='VZ2m')

# two-qubit gates
CPh = CPHASE()
iSWAP_without_Z = iSWAP_no_1qb_phases()

# Composite gates
CZEcho = CompositeGate(n_qubit=2)
CZEcho.add_gate([X2p, I])
CZEcho.add_gate(CPh)
CZEcho.add_gate([Xp, Xp])
CZEcho.add_gate(CPh)
CZEcho.add_gate([X2p, Xp])

H = CompositeGate(n_qubit=1, name='H')
H.add_gate(VZp)
H.add_gate(Y2p)

# populated by sequence.py with keys = #m-n', with m, n qubit numbers
CZ_dict = dict()
iSWAP_dict = dict()

# references for creating lookup tables
CZ_dict['default'] = CPHASE_with_1qb_phases(0, 0)
iSWAP_dict['default'] = iSWAP_with_1qb_phases(0, 0)

# discontinued: a different CNOT gate need to be defined per qubit pair
#CNOT = CompositeGate(n_qubit=2, name='CNOT')
#CNOT.add_gate(H, 1)
#CNOT.add_gate(CZ, [0, 1])
#CNOT.add_gate(H, 1)

if __name__ == '__main__':
    pass
