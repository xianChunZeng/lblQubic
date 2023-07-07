import matplotlib.pyplot as plt

from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
from collections import OrderedDict
from qubic.job_manager_jpm import JobManager
import numpy as np

from scipy.optimize import curve_fit

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

class CrossResonanceTomography(): #AbstractTomographyExperiment
    def __init__(self, control_qid, target_qid, pulse_width_interval, drive_amp, pulse_type='std'):
        self.target_qid = target_qid
        self.control_qid = control_qid
        self.readout_register = [control_qid, target_qid]

        self.pulse_width_interval = pulse_width_interval
        self.drive_amp = drive_amp

        if pulse_type == 'std':
            self.circuits = self._make_circuits()
        elif pulse_type == 'echoed':
            self.circuits = self._make_circuits_echoed()
        self.estimated_hamiltonian_params = None

        self.prior_fits = [
            [2e8, 2e8, 2e8],
            [2e8, 2e8, 2e8]
        ]

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment and report
        """
        self.run(jobmanager, num_shots_per_circuit, qchip)
        self.report()

        

        # fit = self._fit_data(data)
        
    def run(self, jobmanager, num_shots_per_circuit, qchip):
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        self.data = data
        # tomographic curves of the target qubit for analysis and plotting
        # stored in the format: [tomo_axis_idx, twidth_idx, amp_idx]
        tomographic_curves = np.zeros(
            (6, len(self.pulse_width_interval)))
        for twidth_idx in range(len(self.pulse_width_interval)):
            for pauli_idx in range(6):
                tomographic_curves[pauli_idx, twidth_idx] = \
                    1-2 * np.average(data[self.target_qid][twidth_idx*6+pauli_idx])
        self.tomographic_curves = tomographic_curves


        fits = self._fit_data(tomographic_curves, self.prior_fits)
        self.fits = fits
        rates = get_interation_rates(fits[0][0], fits[1][0])
        self.rates = rates
        
    def report(self):
        fig, axs = plt.subplots(4)
        tomographic_curves = self.tomographic_curves
        fits = self.fits

        # scatter plots of collected data
        axs[0].scatter(self.pulse_width_interval, tomographic_curves[0, :], label='X0')
        axs[0].scatter(self.pulse_width_interval, tomographic_curves[3, :], label='X1')
        axs[1].scatter(self.pulse_width_interval, tomographic_curves[1, :], label='Y0')
        axs[1].scatter(self.pulse_width_interval, tomographic_curves[4, :], label='Y1')
        axs[2].scatter(self.pulse_width_interval, tomographic_curves[2, :], label='Z0')
        axs[2].scatter(self.pulse_width_interval, tomographic_curves[5, :], label='Z1')

        # plots of predicted curves with extracted Hamiltonian rates
        axs[0].plot(self.pulse_width_interval,  avg_X(self.pulse_width_interval, *fits[0][0]))
        axs[0].plot(self.pulse_width_interval, avg_X(self.pulse_width_interval, *fits[1][0]))
        axs[1].plot(self.pulse_width_interval,  avg_Y(self.pulse_width_interval, *fits[0][0]))
        axs[1].plot(self.pulse_width_interval, avg_Y(self.pulse_width_interval, *fits[1][0]))
        axs[2].plot(self.pulse_width_interval,  avg_Z(self.pulse_width_interval, *fits[0][0]))
        axs[2].plot(self.pulse_width_interval, avg_Z(self.pulse_width_interval, *fits[1][0]))

        axs[0].set_title('X0 and X1 response')
        axs[1].set_title('Y0 and Y1 response')
        axs[2].set_title('Z0 and Z1 response')
        td = axs[3].plot(self.pulse_width_interval, self.trace_difference(tomographic_curves), label='td')
        rdiff = axs[3].plot(self.pulse_width_interval, self.r_diff(tomographic_curves), label='rdiff')
        axs[3].set_xlabel('pulse width')

        for i in range(4):
            axs[i].legend()
        plt.tight_layout()
            
            
            
    def _fit_data(self, tomographic_curves, prior_fits):
        fit_ctl0 = fit_rt_evol(self.pulse_width_interval, tomographic_curves[0],
                               tomographic_curves[1], tomographic_curves[2], prior_fits[0])

        fit_ctl1 = fit_rt_evol(self.pulse_width_interval, tomographic_curves[3],
                               tomographic_curves[4], tomographic_curves[5], prior_fits[1])
        return [fit_ctl0, fit_ctl1]


        
    @staticmethod
    def r_diff(tomo_arr):
        """
        returns the r-value curve given the two tomographic curves
        """
        return 0.5 * np.sqrt((tomo_arr[3, :] - tomo_arr[0, :]) ** 2 +
                             (tomo_arr[4, :] - tomo_arr[1, :]) ** 2 +
                             (tomo_arr[5, :] - tomo_arr[2, :]) ** 2
                             )
    @staticmethod
    def trace_difference(tomo_arr):
        """
        Calculate the trace difference between the two target qubit states
        """
        return 0.5 * (abs(tomo_arr[3, :] - tomo_arr[0, :]) +
                      abs(tomo_arr[4, :] - tomo_arr[1, :]) +
                      abs(tomo_arr[5, :] - tomo_arr[2, :]))

    def _make_circuits_echoed(self, crosstalk_drive_amp=0.0):
        circuits = []
        control_qubit = self.control_qid
        target_qubit = self.target_qid
        for id_twidth, twidth in enumerate(self.pulse_width_interval):
            for ctrl_state in [0, 1]:
                for axis in ('X', 'Y', 'Z'):
                    circ = [{'name': 'delay', 't': 400.e-6}]
                    if ctrl_state == 1:
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'CR_crosstalk', 'qubit': [control_qubit, target_qubit],
                                    'modi': {(0, 'twidth'): twidth,
                                             (1, 'twidth'): twidth,
                                             (0, 'amp'): self.drive_amp,
                                             (1, 'amp'): crosstalk_drive_amp}})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'CR_crosstalk', 'qubit': [control_qubit, target_qubit],
                                    'modi': {(0, 'twidth'): twidth,
                                             (1, 'twidth'): twidth,
                                             (0, 'amp'): self.drive_amp,
                                             (1, 'amp'): crosstalk_drive_amp,
                                             (0, 'pcarrier'): np.pi,
                                             (1, 'pcarrier'): np.pi}})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    if axis == 'X':
                        circ.append({'name': 'Y-90', 'qubit': [target_qubit]})

                    elif axis == 'Y':
                        circ.append({'name': 'X90', 'qubit': [target_qubit]})
                    circ.append({'name': 'read', 'qubit': [control_qubit]})
                    circ.append({'name': 'read', 'qubit': [target_qubit]})
                    circuits.append(circ)
        return circuits

    # %%

    def _make_circuits(self, crosstalk_drive_amp=0.01):
        """
        makes and returns the circuits

        circuits are stored as an OrderedDict, with keys (id_twidth, id_amp)
        that correspond to the indices in self.drive_amp_interval and self.pulse_width_interval

        all the default circuit parameters are stored in qchip
        and any changes are stored as class properties set at initialization
        note that the qchip is not stored in a calibration experiment,
        it is managed by the jobmanager, and a calibration class can update its parameters
        """

        circuits = []
        control_qubit = self.control_qid
        target_qubit = self.target_qid
        for id_twidth, twidth in enumerate(self.pulse_width_interval):
            amp_time_circuit_list = []  # experiments for a fixed time and amplitude
            for ctrl_state in [0, 1]:
                for axis in ['X', 'Y', 'Z']:
                    circ = [{'name': 'delay', 't': 400.e-6}]
                    circ = [{'name': 'delay', 't': 400.e-6}]
                    if ctrl_state == 1:
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                    circ.append({'name': 'CR_crosstalk', 'qubit': [control_qubit, target_qubit],
                                 'modi': {(0, 'twidth'): twidth,
                                          (1, 'twidth'): twidth,
                                          (0, 'amp'): self.drive_amp,
                                          (1, 'amp'): crosstalk_drive_amp}})
                    if axis == 'X':
                            circ.append({'name': 'Y-90', 'qubit': [target_qubit]})
                    elif axis == 'Y':
                        circ.append({'name': 'X90', 'qubit': [target_qubit]})
                    circ.append({'name': 'read', 'qubit': [control_qubit]})
                    circ.append({'name': 'read', 'qubit': [target_qubit]})
                    circuits.append(circ)
        return circuits

    def _collect_data(self, jobmanager: JobManager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        the GMMM and the FPGA/Channel configs and the qchip is managed
        """
        return jobmanager.collect_classified_shots(self.circuits, num_shots_per_circuit, qchip)

#====================================================================

class CrossResonanceCalibration(AbstractCalibrationExperiment):
    """
    Calibrate the cross resonance gate's time and amplitude
    """

    def __init__(self, control_qid, target_qid, pulse_width_interval, drive_amp_interval, pulse_type='std'):
        self.target_qid = target_qid
        self.control_qid = control_qid
        self.readout_register = [control_qid, target_qid]

        self.pulse_width_interval = pulse_width_interval
        self.drive_amp_interval = drive_amp_interval

        if pulse_type == 'std':
            self.circuits = self._make_circuits()
        elif pulse_type == 'echoed':
            self.circuits = self._make_circuits_echoed()
        self.optimization_parameters = ['CR.amp', 'CR.twidth']
        self.final_estimated_params = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment and report
        """
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        self.data = data
        # tomographic curves of the target qubit for analysis and plotting
        # stored in the format: [tomo_axis_idx, twidth_idx, amp_idx]
        tomographic_curves = np.zeros(
            (6, len(self.pulse_width_interval), len(self.drive_amp_interval)))
        for index_pair in data.keys():
            tomographic_curves[:, index_pair[0], index_pair[1]] = \
                1-2 * np.average(data[index_pair][self.target_qid], axis=1)[:, 0]
        self.tomographic_curves = tomographic_curves

        # fit = self._fit_data(data)

        if len(self.drive_amp_interval) == 1:
            amp = self.drive_amp_interval[0]
            fig, axs = plt.subplots(4)
            axs[0].plot(self.pulse_width_interval, tomographic_curves[0, :, 0], label='X0')
            axs[0].plot(self.pulse_width_interval, tomographic_curves[3, :, 0], label='X1')
            axs[1].plot(self.pulse_width_interval, tomographic_curves[1, :, 0], label='Y0')
            axs[1].plot(self.pulse_width_interval, tomographic_curves[4, :, 0], label='Y1')
            axs[2].plot(self.pulse_width_interval, tomographic_curves[2, :, 0], label='Z0')
            axs[2].plot(self.pulse_width_interval, tomographic_curves[5, :, 0], label='Z1')
            axs[0].set_title('X0 and X1 response')
            axs[1].set_title('Y0 and Y1 response')
            axs[2].set_title('Z0 and Z1 response')
            td = axs[3].plot(self.pulse_width_interval, trace_difference(tomographic_curves), label='td')
            rdiff = axs[3].plot(self.pulse_width_interval, r_diff(tomographic_curves), label='rdiff')
            axs[3].set_xlabel('pulse width')
            
            for i in range(4):
                axs[i].legend()
            plt.tight_layout()
            
        elif len(self.pulse_width_interval) == 1:
            twidth = self.pulse_width_interval[0]
            fig, axs = plt.subplots(4)
            axs[0].plot(self.drive_amp_interval, tomographic_curves[0, 0, :], label='X0')
            axs[0].plot(self.drive_amp_interval, tomographic_curves[3, 0, :], label='X1')
            axs[1].plot(self.drive_amp_interval, tomographic_curves[1, 0, :], label='Y0')
            axs[1].plot(self.drive_amp_interval, tomographic_curves[4, 0, :], label='Y1')
            axs[2].plot(self.drive_amp_interval, tomographic_curves[2, 0, :], label='Z0')
            axs[2].plot(self.drive_amp_interval, tomographic_curves[5, 0, :], label='Z1')
            td = axs[3].plot(self.drive_amp_interval, trace_difference(tomographic_curves)[0, :], label='td')
            rdiff = axs[3].plot(self.drive_amp_interval, r_diff(tomographic_curves)[0, :], label='rdiff')
            axs[3].set_xlabel('drive amp')
            axs[0].set_title('X0 and X1 response')
            axs[1].set_title('Y0 and Y1 response')
            axs[2].set_title('Z0 and Z1 response')
            for i in range(4):
                axs[i].legend()
            plt.tight_layout()
        else: # plot a 2d color grid
            X, Y = np.meshgrid(self.pulse_width_interval, self.drive_amp_interval, indexing='ij')
            fig, axs = plt.subplots(3)
            cmX0 = axs[0].pcolormesh(X, Y, tomographic_curves[0, :, :], label='X0', cmap='Blues', alpha=0.5)
            cmX1 = axs[0].pcolormesh(X, Y, tomographic_curves[3, :, :], label='X1', cmap='Reds', alpha=0.5)
            cmY0 = axs[1].pcolormesh(X, Y, tomographic_curves[1, :, :], label='Y0', cmap='Blues', alpha=0.5)
            cmY1 = axs[1].pcolormesh(X, Y, tomographic_curves[4, :, :], label='Y1', cmap='Reds', alpha=0.5)
            cmZ0 = axs[2].pcolormesh(X, Y, tomographic_curves[2, :, :], label='Z0', cmap='Blues', alpha=0.5)
            cmZ1 = axs[2].pcolormesh(X, Y, tomographic_curves[5, :, :], label='Z1', cmap='Reds', alpha=0.5)
            
            fig.colorbar(cmX0, ax=axs[0])
            fig.colorbar(cmX1, ax=axs[0])
            fig.colorbar(cmY0, ax=axs[1])
            fig.colorbar(cmY1, ax=axs[1])
            fig.colorbar(cmZ0, ax=axs[2])
            fig.colorbar(cmZ1, ax=axs[2])
            for i in range(3):
                axs[i].set_ylim([self.drive_amp_interval[0], self.drive_amp_interval[-1]])
            axs[0].set_title('X0 (blue) and X1 (red) response')
            axs[1].set_title('Y0 (blue) and Y1 (red) response')
            axs[2].set_title('Z0 (blue) and Z1 (red) response')
            plt.tight_layout()
                
                
            plt.show()
            
            fig, axs = plt.subplots(2)
            cmTD = axs[0].pcolormesh(X, Y, trace_difference(tomographic_curves), cmap='Greens')
            axs[0].set_ylim([self.drive_amp_interval[0], self.drive_amp_interval[-1]])
            axs[0].set_ylabel('drive amp')
            axs[0].set_xlabel('pulse width')
            axs[0].set_title('trace distance between target tomography states')
            cmRD = axs[1].pcolormesh(X, Y, r_diff(tomographic_curves), cmap='Greens')
            axs[1].set_ylim([self.drive_amp_interval[0], self.drive_amp_interval[-1]])
            axs[1].set_ylabel('drive amp')
            axs[1].set_xlabel('pulse width')
            axs[1].set_title('r-distance between target tomography states')
            fig.colorbar(cmTD, ax=axs[0])
            fig.colorbar(cmRD, ax=axs[1])
            plt.tight_layout()
            plt.show()
            

        # plt.plot(fit)
        # make report and save

        self.final_estimated_params = 1  # find the estimate based on the experiment type, the data, and the fit
        return self.final_estimated_params

    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        pass


    
    def _make_circuits_echoed(self, crosstalk_drive_amp=0.01):
        circuits = OrderedDict()
        control_qubit = self.control_qid
        target_qubit = self.target_qid
        for id_twidth, twidth in enumerate(self.pulse_width_interval):
            for id_amp, amp in enumerate(self.drive_amp_interval):
                amp_time_circuit_list = []  # experiments for a fixed time and amplitude
                for ctrl_state in [0, 1]:
                    for axis in ('X', 'Y', 'Z'):
                        circ = [{'name': 'delay', 't': 400.e-6}]
                        if ctrl_state == 1:
                            circ.append({'name': 'X90', 'qubit': [control_qubit]})
                            circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'CR_crosstalk', 'qubit': [control_qubit, target_qubit],
                                        'modi': {(0, 'twidth'): twidth,
                                                 (1, 'twidth'): twidth,
                                                 (0, 'amp'): amp,
                                                 (1, 'amp'): crosstalk_drive_amp}})
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'CR_crosstalk', 'qubit': [control_qubit, target_qubit],
                                        'modi': {(0, 'twidth'): twidth,
                                                 (1, 'twidth'): twidth,
                                                 (0, 'amp'): amp,
                                                 (1, 'amp'): crosstalk_drive_amp,
                                                 (0, 'pcarrier'): np.pi,
                                                 (1, 'pcarrier'): np.pi}})
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        if axis == 'X':
                            circ.append({'name': 'Y-90', 'qubit': [target_qubit]})

                        elif axis == 'Y':
                            circ.append({'name': 'X90', 'qubit': [target_qubit]})
                        circ.append({'name': 'read', 'qubit': [control_qubit]})
                        circ.append({'name': 'read', 'qubit': [target_qubit]})
                        amp_time_circuit_list.append(circ)
                circuits[(id_twidth, id_amp)] = amp_time_circuit_list
        return circuits

    # %%

    def _make_circuits(self, crosstalk_drive_amp=0.01):
        """
        makes and returns the circuits

        circuits are stored as an OrderedDict, with keys (id_twidth, id_amp)
        that correspond to the indices in self.drive_amp_interval and self.pulse_width_interval

        all the default circuit parameters are stored in qchip
        and any changes are stored as class properties set at initialization
        note that the qchip is not stored in a calibration experiment,
        it is managed by the jobmanager, and a calibration class can update its parameters
        """

        circuits = OrderedDict()
        control_qubit = self.control_qid
        target_qubit = self.target_qid
        for id_twidth, twidth in enumerate(self.pulse_width_interval):
            for id_amp, amp in enumerate(self.drive_amp_interval):
                amp_time_circuit_list = []  # experiments for a fixed time and amplitude
                for ctrl_state in [0, 1]:
                    for axis in ['X', 'Y', 'Z']:
                        circ = [{'name': 'delay', 't': 400.e-6}]
                        circ = [{'name': 'delay', 't': 400.e-6}]
                        if ctrl_state == 1:
                            circ.append({'name': 'X90', 'qubit': [control_qubit]})
                            circ.append({'name': 'X90', 'qubit': [control_qubit]})
                        circ.append({'name': 'CR_crosstalk', 'qubit': [control_qubit, target_qubit],
                                     'modi': {(0, 'twidth'): twidth,
                                              (1, 'twidth'): twidth,
                                              (0, 'amp'): amp,
                                              (1, 'amp'): crosstalk_drive_amp}})
                        circ.append({'name': 'read', 'qubit': [control_qubit]})
                        circ.append({'name': 'read', 'qubit': [target_qubit]})
                        amp_time_circuit_list.append(circ)
                circuits[(id_twidth, id_amp)] = amp_time_circuit_list
        return circuits

    def _collect_data(self, jobmanager: JobManager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        the GMMM and the FPGA/Channel configs and the qchip is managed
        """
        collected_shots = OrderedDict()
        for index_pairs in self.circuits.keys():
            collected_shots[index_pairs] = jobmanager.collect_classified_shots(self.circuits[index_pairs],
                                                                               num_shots_per_circuit, qchip)
        return collected_shots
