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

    def __init__(self, control_qid, target_qid, time_interval, calibration, cr_type='std'):
        self.time_interval = time_interval
        self.control_amp = calibration['camp']
        self.target_amp = calibration['tamp']
        self.control_phase = calibration['cphase']
        self.target_phase = calibration['tphase']
        self.control_qid = control_qid
        self.target_qid = target_qid
        self.cr_type = cr_type

    def run(self, jobmanager, num_shots_per_circuit):
        # make instructions
        if self.cr_type == 'std':
            circ_dict = self.make_tomo_instructions()
        elif self.cr_type == 'echoed':
            circ_dict = self.make_echoed_tomo_instructions()
        self.circuits = circ_dict
        # collect raw shot data
        shot_dict = self.run_tomo_circuits(circ_dict, jobmanager, num_shots_per_circuit)
        self.shot_dict = shot_dict
        # format and store the target qubit's shots
        qtarget_shots = np.zeros((len(self.time_interval), 2, 3, num_shots_per_circuit))
        for id_key, key in enumerate(circ_dict.keys()):
            qtarget_shots[key[0], key[1], key[2], :] = shot_dict[self.target_qid][key].flatten()
        self.qtarget_shots = qtarget_shots
        # calculate tomographic response on qtarget
        curves = 1-2*np.average(qtarget_shots, axis=-1)
        self.curves = curves
        
        return shot_dict, shot_dict, curves, circ_dict

    def make_cr_pulse(self, twidth, phase_offset=0):
        return [
            {
            'name': 'pulse',
            "freq": f"{self.target_qid}.freq",
            "dest": f"{self.control_qid}.qdrv",
            "twidth": twidth,
            "amp": self.control_amp,
            "phase": self.control_phase + phase_offset,
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
            "dest": f"{self.target_qid}.qdrv",
            "twidth": twidth,
            "amp": self.target_amp,
            "phase": self.target_phase + phase_offset, 
            "env": [
                {
                    "env_func": "cos_edge_square",
                    "paradict": {
                        "ramp_fraction": 0.25,
                    }
                }
            ]
            }, {'name': 'barrier'}
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
                    circ = circ + self.make_cr_pulse(twidth)
                    if axis == 'X':
                        circ.append({'name': 'Y-90', 'qubit': [target_qubit]})
                    elif axis == 'Y':
                        circ.append({'name': 'X90', 'qubit': [target_qubit]})
                    circ.append({'name': 'read', 'qubit': [control_qubit]})
                    circ.append({'name': 'read', 'qubit': [target_qubit]})
                    circuits[(id_twidth, ctrl_state, id_axis)] = circ
        return circuits
    
    def make_echoed_tomo_instructions(self):
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
                    circ = circ + self.make_cr_pulse(twidth/2)
                    circ.append({'name': 'delay', 't': 10.e-9})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'barrier'})
                    circ = circ + self.make_cr_pulse(twidth/2, np.pi)
                    circ.append({'name': 'delay', 't': 10.e-9})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'barrier'})
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
        shot_dict = {qid : { key: shots[qid][idx, :] for idx, key in enumerate(circuits.keys()) } for qid in [self.control_qid, self.target_qid]}
        return shot_dict

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

def raw_curve_fun(x, A, B, C, period):
    return A + B*np.cos(2*np.pi*x/period) + C*np.sin(2*np.pi*x/period)

def raw_curve_fits_to_hmodel(xfit, yfit, zfit):
    omega = (1/3)*2*np.pi*(1/xfit[3] + 1/yfit[3] + 1/zfit[3])
    omega_y = xfit[2]*omega
    omega_x = -yfit[1]*omega
    comp_omega = omega_x**2 + omega_y**2
    delta = np.sqrt(abs(omega**2 - comp_omega))
    if (omega_x > 0 and xfit[1] < 0) or (omega_x < 0 and xfit[1] > 0):
        delta = -delta
        
    return delta, omega_x, omega_y


def r_diff(tomo_arr):
    """
    returns the r-value curve given the two tomographic curves
    """
    return 0.5 * np.sqrt((tomo_arr[:, 1, 0] - tomo_arr[:, 0, 0]) ** 2 +
                         (tomo_arr[:, 1, 1] - tomo_arr[:, 0, 1]) ** 2 +
                         (tomo_arr[:, 1, 2] - tomo_arr[:, 0, 2]) ** 2
                         )


def trace_difference(tomo_arr):
    """
    Calculate the trace difference between the two target qubit states
    """
    return 0.5 * (abs(tomo_arr[:, 1, 0] - tomo_arr[:, 0, 0]) +
                  abs(tomo_arr[:, 1, 1] - tomo_arr[:, 0, 1]) +
                  abs(tomo_arr[:, 1, 2] - tomo_arr[:, 0, 2])
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
    def __init__(self, cr_tomo_experiment, priors, fit_type='std'):
        self.tomo_curves = cr_tomo_experiment.curves
        self.twidth_interval = cr_tomo_experiment.time_interval
        self.has_raw_fits = False
        if fit_type == 'std':
            try:
                self.fits = self._model_fit(priors)
                c0_scaled_params = np.array(self.fits[0])*1e7
                c1_scaled_params = np.array(self.fits[1])*1e7
                self.scaled_params = (c0_scaled_params, c1_scaled_params)
                self.rates = get_interation_rates(c0_scaled_params, c1_scaled_params)
            except:
                raise RuntimeError("Could not fit the curves with the provided priors")
#         elif fit_type == 'raw':
#             self.has_raw_fits = True
#             try:
#                 self.fits = self._raw_sin_fit(priors)
#             except:
#                 raise RuntimeError("Could not fit the curves with the provided priors")
#         elif fit_type == 'interactive':
#             self.has_raw_fits = True
#             try:
#                 raw_fits = self._raw_sin_fit(priors)
#                 self.fits = raw_fits
#             except:
#                 print("Could not fit the curves with the provided priors")
#             looping = True
#             while(looping):
#                 self.plot()
#                 plt.show()
#                 x = input('Is the fit correct? Either input yes or no')
#                 if x == 'yes':
#                     looping = False
#                 else:
#                     print(raw_fits[0][0])
#                     print(raw_fits[1][0])
#                     p0 = input('new prior for prep 0')
#                     p1 = input('new prior for prep 1')
#                     try:
#                         self.fits = self._model_fit([p0, p1])
#                     except:
#                         print("Could not fit the curves with the provided priors")

#     def _raw_sin_fit(self, priors):
#         # first do a raw fit of the curves with the given prior period
#         priors_0 = priors[0]
#         priors_1 = priors[1]
#         fits0 = [curve_fit(raw_curve_fun, self.twidth_interval, self.tomo_curves[:, 0, 0], p0=priors_0),
#             curve_fit(raw_curve_fun, self.twidth_interval, self.tomo_curves[:, 0, 1], p0=priors_0),
#             curve_fit(raw_curve_fun, self.twidth_interval, self.tomo_curves[:, 0, 2], p0=priors_0)]
        
#         fits1 = [curve_fit(raw_curve_fun, self.twidth_interval, self.tomo_curves[:, 1, 0], p0=priors_1),
#             curve_fit(raw_curve_fun, self.twidth_interval, self.tomo_curves[:, 1, 1], p0=priors_1),
#             curve_fit(raw_curve_fun, self.twidth_interval, self.tomo_curves[:, 1, 2], p0=priors_1)]
#         hmodel0 = raw_curve_fits_to_hmodel(*[f[0] for f in fits0])
#         hmodel1 = raw_curve_fits_to_hmodel(*[f[0] for f in fits1])
        
#         self.raw_fits = [fits0, fits1]
                 
#         fit_ctl0 = fit_rt_evol(self.twidth_interval*1e7, self.tomo_curves[:, 0, 0],
#                                self.tomo_curves[:, 0, 1], self.tomo_curves[:, 0, 2], hmodel0)

#         fit_ctl1 = fit_rt_evol(self.twidth_interval, self.tomo_curves[:, 1, 0],
#                                self.tomo_curves[:, 1, 1], self.tomo_curves[:, 1, 2], hmodel1)
#         return [fit_ctl0, fit_ctl1]
    
    def _model_fit(self, priors):
        fit_ctl0, _ = fit_rt_evol(self.twidth_interval*1e7, self.tomo_curves[:, 0, 0],
                               self.tomo_curves[:, 0, 1], self.tomo_curves[:, 0, 2], priors[0])
        fit_ctl0 = np.array(fit_ctl0)

        fit_ctl1, _ = fit_rt_evol(self.twidth_interval*1e7, self.tomo_curves[:, 1, 0],
                               self.tomo_curves[:, 1, 1], self.tomo_curves[:, 1, 2], priors[1])
        fit_ctl1 = np.array(fit_ctl1)
        return [fit_ctl0, fit_ctl1]
    
    def plot(self, show_fits=True, fit_res=100):
        fig, axs = plt.subplots(4)

        # scatter plots of collected data
        axs[0].scatter(self.twidth_interval, self.tomo_curves[:, 0, 0], label='X0')
        axs[0].scatter(self.twidth_interval, self.tomo_curves[:, 1, 0], label='X1')
        axs[1].scatter(self.twidth_interval, self.tomo_curves[:, 0, 1], label='Y0')
        axs[1].scatter(self.twidth_interval, self.tomo_curves[:, 1, 1], label='Y1')
        axs[2].scatter(self.twidth_interval, self.tomo_curves[:, 0, 2], label='Z0')
        axs[2].scatter(self.twidth_interval, self.tomo_curves[:, 1, 2], label='Z1')
        
        if show_fits:
            # plots of predicted curves with extracted Hamiltonian rates
            interval = np.linspace(self.twidth_interval[0], self.twidth_interval[-1], fit_res)
            axs[0].plot(interval, avg_X(interval, *self.scaled_params[0]))
            axs[0].plot(interval, avg_X(interval, *self.scaled_params[1]))
            axs[1].plot(interval, avg_Y(interval, *self.scaled_params[0]))
            axs[1].plot(interval, avg_Y(interval, *self.scaled_params[1]))
            axs[2].plot(interval, avg_Z(interval, *self.scaled_params[0]))
            axs[2].plot(interval, avg_Z(interval, *self.scaled_params[1]))
#         if self.has_raw_fits:
#             # plots of predicted curves with raw oscillations
#             interval = np.linspace(self.twidth_interval[0], self.twidth_interval[-1], fit_res)
#             axs[0].plot(interval, raw_curve_fun(interval, *self.raw_fits[0][0][0]), linestyle='dashed')
#             axs[0].plot(interval, raw_curve_fun(interval, *self.raw_fits[1][0][0]), linestyle='dashed')
#             axs[1].plot(interval, raw_curve_fun(interval, *self.raw_fits[0][1][0]), linestyle='dashed')
#             axs[1].plot(interval, raw_curve_fun(interval, *self.raw_fits[1][1][0]), linestyle='dashed')
#             axs[2].plot(interval, raw_curve_fun(interval, *self.raw_fits[0][2][0]), linestyle='dashed')
#             axs[2].plot(interval, raw_curve_fun(interval, *self.raw_fits[1][2][0]), linestyle='dashed')

        axs[0].set_title('X0 and X1 response')
        axs[1].set_title('Y0 and Y1 response')
        axs[2].set_title('Z0 and Z1 response')
        td = axs[3].plot(self.twidth_interval, trace_difference(self.tomo_curves), label='td')
        rdiff = axs[3].plot(self.twidth_interval, r_diff(self.tomo_curves), label='rdiff')
        axs[3].set_xlabel('pulse width')

        for i in range(4):
            axs[i].legend()
        plt.tight_layout()