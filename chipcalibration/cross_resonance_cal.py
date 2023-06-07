from collections import OrderedDict
from qubic.job_manager import JobManager
import numpy as np

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class CRTomographyExperiment:
    """
    Experiment Class that takes data and stores it for analysis by other classes

    calibrations:
        control_amp
        target_amp
        control_phase
        target_phase
    time_interval contains the twidths of the cr pulse
    """

    def __init__(self, time_interval, control_qid, target_qid, calibration):
        self.time_interval = time_interval
        self.control_amp = calibration['camp']
        self.target_amp = calibration['tamp']
        self.control_phase = calibration['cphase']
        self.target_phase = calibration['tphase']
        self.control_qid = control_qid
        self.target_qid = target_qid

    def run(self, jobmanager, num_shots_per_circuit):
        circ_dict = self.make_tomo_instructions()
        shots = self.run_tomo_circuits(circ_dict, jobmanager, num_shots_per_circuit)
        curves = np.average(shots, axis=-1)
        self.shots = shots
        self.curves = curves
        self.circuits = circ_dict
        return shots, curves, circ_dict

    def make_cr_pulse(self, twidth):
        return [{
            'name': 'pulse',
            "freq": f"{self.target_qid}.freq",
            "phase": self.control_phase,
            "dest": f"{self.control_qid}.qdrv",
            "twidth": twidth,
            "amp": self.control_amp,
            "env": [
                {
                    "env_func": "cos_edge_square",
                    "paradict": {
                        "ramp_fraction": 0.25,
                    }
                }
            ]
        }, {
            'name': 'pulse',
            "freq": f"{self.target_qid}.freq",
            "phase": self.target_phase,
            "dest": f"{self.target_qid}.qdrv",
            "twidth": 4e-07,
            "amp": self.target_amp,
            "env": [
                {
                    "env_func": "cos_edge_square",
                    "paradict": {
                        "ramp_fraction": 0.25,
                    }
                }
            ]
        }

        ]

    def make_tomo_instructions(self):
        """
        makes tomographic circuits instructions for the CR gate and return them

        circuits are stored as an ordered dict,
        indexed by the twidth index, the control preparation, and the Pauli observable
        """

        circuits = OrderedDict()
        control_qubit = self.control_qid
        target_qubit = self.target_qid
        for id_twidth, twidth in enumerate(self.time_interval):
            for ctrl_state in [0, 1]:
                for id_axis, axis in enumerate(['X', 'Y', 'Z']):
                    circ = [{'name': 'delay', 't': 400.e-6}]
                    if ctrl_state == 1:
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append(self.make_cr_pulse(twidth))
                    if axis == 'X':
                        circ.append({'name': 'Y-90', 'qubit': [target_qubit]})
                    elif axis == 'Y':
                        circ.append({'name': 'X90', 'qubit': [target_qubit]})
                    circ.append({'name': 'read', 'qubit': [control_qubit]})
                    circ.append({'name': 'read', 'qubit': [target_qubit]})
                    circuits[(id_twidth, ctrl_state, id_axis)] = circ
        return circuits

    def run_tomo_circuits(self, circuits, jobmanager: JobManager, num_shots_per_circuit):
        """
        runs the tomographic circuits and returns an ndarray of format
        [twidth_idx, control_state, pauli_type, shot_idx]
        """
        shots = jobmanager.collect_classified_shots(list(circuits.values()), num_shots_per_circuit)
        formatted_shots = np.zeros((len(self.time_interval), 2, 3, num_shots_per_circuit))
        for id_key, key in enumerate(circuits.keys()):
            formatted_shots[key[0], key[1], key[2], :] = shots[id_key]
        return formatted_shots

#=== Analyzer ====================================================================================



def get_omega(eDelta, eOmega_x, eOmega_y):
    """Return \Omega from parameter arguments."""
    eOmega = np.sqrt(eDelta ** 2 + eOmega_x ** 2 + eOmega_y ** 2)
    return eOmega


def avg_X(t, eDelta, eOmega_x, eOmega_y):
    """Return average X Pauli measurement vs time t"""
    eOmega = get_omega(eDelta, eOmega_x, eOmega_y)
    eXt = (-eDelta * eOmega_x + eDelta * eOmega_x * np.cos(eOmega * t) + \
           eOmega * eOmega_y * np.sin(eOmega * t)) / eOmega ** 2
    return eXt


def avg_Y(t, eDelta, eOmega_x, eOmega_y):
    """Return average Y Pauli measurement vs time t"""
    eOmega = get_omega(eDelta, eOmega_x, eOmega_y)
    eYt = (eDelta * eOmega_y - eDelta * eOmega_y * np.cos(eOmega * t) - \
           eOmega * eOmega_x * np.sin(eOmega * t)) / eOmega ** 2
    return eYt


def avg_Z(t, eDelta, eOmega_x, eOmega_y):
    """Return average Z Pauli measurement vs time t"""
    eOmega = get_omega(eDelta, eOmega_x, eOmega_y)
    eZt = (eDelta ** 2 + (eOmega_x ** 2 + eOmega_y ** 2) * np.cos(eOmega * t)) / eOmega ** 2
    return eZt


def rt_evol(ts, eDelta, eOmega_x, eOmega_y):
    """Stack average X,Y,Z Pauli measurements vertically."""
    return np.vstack([avg_X(ts, eDelta, eOmega_x, eOmega_y), \
                      avg_Y(ts, eDelta, eOmega_x, eOmega_y), \
                      avg_Z(ts, eDelta, eOmega_x, eOmega_y)])


def rt_flat(ts, eDelta, eOmega_x, eOmega_y):
    """Flatten X,Y,Z Pauli measurement data into 1D array."""
    return rt_evol(ts[0:len(ts) // 3], eDelta, eOmega_x, eOmega_y).flatten()


def fit_rt_evol(ts, eXt, eYt, eZt, p0):
    """Use curve_fit to determine fit parameters of X,Y,Z Pauli measurements together."""
    rt_vec = np.asarray([eXt, eYt, eZt])

    return curve_fit(rt_flat, np.tile(ts, 3), rt_vec.flatten(), p0=p0, method='trf')


def r_diff(tomo_arr):
    """
    returns the r-value curve given the two tomographic curves
    """
    return 0.5 * np.sqrt((tomo_arr[3, :, :] - tomo_arr[0, :, :]) ** 2 +
                         (tomo_arr[4, :, :] - tomo_arr[1, :, :]) ** 2 +
                         (tomo_arr[5, :, :] - tomo_arr[2, :, :]) ** 2
                         )


def trace_difference(tomo_arr):
    """
    Calculate the trace difference between the two target qubit states
    """
    return 0.5 * (abs(tomo_arr[3, :, :] - tomo_arr[0, :, :]) +
                  abs(tomo_arr[4, :, :] - tomo_arr[1, :, :]) +
                  abs(tomo_arr[5, :, :] - tomo_arr[2, :, :])
                  )


def get_interation_rates(ground_fit, excited_fit):
    """Determine interaction rates from fits to ground and excited control qubit data."""
    Delta0 = ground_fit[0]
    Omega0_x = ground_fit[1]
    Omega0_y = ground_fit[2]
    Delta1 = excited_fit[0]
    Omega1_x = excited_fit[1]
    Omega1_y = excited_fit[2]

    rates = dict()
    rates['IX'] = 0.5 * (Omega0_x + Omega1_x)
    rates['IY'] = 0.5 * (Omega0_y + Omega1_y)
    rates['IZ'] = 0.5 * (Delta0 + Delta1)
    rates['ZX'] = 0.5 * (Omega0_x - Omega1_x)
    rates['ZY'] = 0.5 * (Omega0_y - Omega1_y)
    rates['ZZ'] = 0.5 * (Delta0 - Delta1)

    return rates

class CRHamiltonianAnalyzer:
    def __init__(self, tomo_curves, twidth_interval, prior_fits):
        self.tomo_curves = tomo_curves
        self.twidth_interval = twidth_interval
        self.prior_fits = prior_fits
        try:
            self.fits = self._fit_data()
        except:
            raise RuntimeError("Could not fit the curves with the provided priors")

    def _fit_data(self):
        fit_ctl0 = fit_rt_evol(self.twidth_interval, self.tomo_curves[0],
                               self.tomo_curves[1], self.tomo_curves[2], self.prior_fits[0])

        fit_ctl1 = fit_rt_evol(self.twidth_interval, self.tomo_curves[3],
                               self.tomo_curves[4], self.tomo_curves[5], self.prior_fits[1])
        return [fit_ctl0, fit_ctl1]
    
    def plot(self, show_fits=True):
        fig, axs = plt.subplots(4)

        # scatter plots of collected data
        axs[0].scatter(self.twidth_interval, self.tomo_curves[0, :], label='X0')
        axs[0].scatter(self.twidth_interval, self.tomo_curves[3, :], label='X1')
        axs[1].scatter(self.twidth_interval, self.tomo_curves[1, :], label='Y0')
        axs[1].scatter(self.twidth_interval, self.tomo_curves[4, :], label='Y1')
        axs[2].scatter(self.twidth_interval, self.tomo_curves[2, :], label='Z0')
        axs[2].scatter(self.twidth_interval, self.tomo_curves[5, :], label='Z1')
        
        if show_fits:
            # plots of predicted curves with extracted Hamiltonian rates
            axs[0].scatter(self.twidth_interval, avg_X(self.twidth_interval, *self.fits[0][0]))
            axs[0].scatter(self.twidth_interval, avg_X(self.twidth_interval, *self.fits[1][0]))
            axs[1].scatter(self.twidth_interval, avg_Y(self.twidth_interval, *self.fits[0][0]))
            axs[1].scatter(self.twidth_interval, avg_Y(self.twidth_interval, *self.fits[1][0]))
            axs[2].scatter(self.twidth_interval, avg_Z(self.twidth_interval, *self.fits[0][0]))
            axs[2].scatter(self.twidth_interval, avg_Z(self.twidth_interval, *self.fits[1][0]))

        axs[0].set_title('X0 and X1 response')
        axs[1].set_title('Y0 and Y1 response')
        axs[2].set_title('Z0 and Z1 response')
        td = axs[3].plot(self.twidth_interval, self.trace_difference(self.tomo_curves), label='td')
        rdiff = axs[3].plot(self.twidth_interval, self.r_diff(self.tomo_curves), label='rdiff')
        axs[3].set_xlabel('pulse width')

        for i in range(4):
            axs[i].legend()
        plt.tight_layout()