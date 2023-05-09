import matplotlib.pyplot as plt

from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
from collections import OrderedDict
from qubic.job_manager_jpm import JobManager
import numpy as np


class CrossResonanceCalibration(AbstractCalibrationExperiment):
    """
    Calibrate the cross resonance gate's time and amplitude
    """

    def __init__(self, control_qid, target_qid, pulse_width_interval, drive_amp_interval):
        self.target_qid = target_qid
        self.control_qid = control_qid
        self.readout_register = [control_qid, target_qid]

        self.pulse_width_interval = pulse_width_interval
        self.drive_amp_interval = drive_amp_interval

        self.circuits = self._make_circuits()
        self.optimization_parameters = ['CR.amp', 'CR.twidth']
        self.final_estimated_params = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment and report
        """
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)

        # tomographic curves of the target qubit for analysis and plotting
        # stored in the format: [tomo_axis_idx, twidth_idx, amp_idx]
        tomographic_curves = np.zeros(
            (6, len(self.pulse_width_interval), len(self.drive_amp_interval)))
        for index_pair in data.keys():
            for id_axis in range(6):
                tomographic_curves[id_axis, index_pair[0], index_pair[1]] = \
                    2 * np.average(data[index_pair][self.target_qid], axis=1)[id_axis] - 1

        # fit = self._fit_data(data)

        assert len(self.pulse_width_interval) != 1, "must be a non-trivial drive interval"
        if len(self.drive_amp_interval) == 1:
            amp = self.drive_amp_interval[0]
            fig, axs = plt.subplots(4)
            axs[0].plot(tomographic_curves[0, :, 0], label='X0')
            axs[0].plot(tomographic_curves[3, :, 0], label='X1')
            axs[1].plot(tomographic_curves[1, :, 0], label='Y0')
            axs[1].plot(tomographic_curves[4, :, 0], label='Y1')
            axs[2].plot(tomographic_curves[2, :, 0], label='Z0')
            axs[2].plot(tomographic_curves[5, :, 0], label='Z1')
            axs[3].plot(self.trace_difference(tomographic_curves))
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
                axs[i].set_ylim([self.drive_amp_interval[0], self.drive_amp_interval[1]])
            axs[0].set_title('X0 (blue) and X1 (red) response')
            axs[1].set_title('Y0 (blue) and Y1 (red) response')
            axs[2].set_title('Z0 (blue) and Z1 (red) response')
            plt.tight_layout()
                
                
            plt.show()
            
            fig, axs = plt.subplots()
            cmTD = axs.pcolormesh(X, Y, self.trace_difference(tomographic_curves), cmap='Greens')
            fig.colorbar(cmTD, ax=axs)
            axs.set_ylim([self.drive_amp_interval[0], self.drive_amp_interval[1]])
            axs.set_ylabel('drive amp')
            axs.set_xlabel('pulse width')
            axs.set_title('trace distance between target tomography states')
            plt.show()
            

        # plt.plot(fit)
        # make report and save

        self.final_estimated_params = 1  # find the estimate based on the experiment type, the data, and the fit
        return self.final_estimated_params

    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        pass

    def r_curve(self, tomo_arr):
        """
        returns the r-value curve given the two tomographic curves
        """
        return 0.5 * np.sqrt((tomo_arr[3, :, :] - tomo_arr[0, :, :]) ** 2 +
                             (tomo_arr[4, :, :] - tomo_arr[1, :, :]) ** 2 +
                             (tomo_arr[5, :, :] - tomo_arr[2, :, :]) ** 2
                             )

    def trace_difference(self, tomo_arr):
        """
        Calculate the trace difference between the two target qubit states
        """
        return 0.5 * (abs(tomo_arr[3, :, :] - tomo_arr[0, :, :]) +
                      abs(tomo_arr[4, :, :] - tomo_arr[1, :, :]) +
                      abs(tomo_arr[5, :, :] - tomo_arr[2, :, :])
                      )

    def _make_circuits(self):
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
        for id_twidth, twidth in enumerate(self.pulse_width_interval):
            for id_amp, amp in enumerate(self.drive_amp_interval):
                amp_time_circuit_list = []  # experiments for a fixed time and amplitude
                for control_state in [0, 1]:
                    for axis in ['X', 'Y', 'Z']:
                        circ = [{'name': 'delay', 't': 400.e-6}]
                        if control_state == 1:
                            circ.append({'name': 'X90', 'qubit': [self.control_qid]})
                            circ.append({'name': 'X90', 'qubit': [self.control_qid]})
                        circ.append({'name': 'CR', 'qubit': [self.control_qid, self.target_qid],
                                     'mode': {(0, 'twidth'): twidth, (0, 'amp'): amp}})
                        if axis == 'X':
                            circ.append({'name': 'Y-90', 'qubit': [self.target_qid]})
                        elif axis == 'Y':
                            circ.append({'name': 'X90', 'qubit': [self.target_qid]})
                        circ.append({'name': 'barrier', 'qubit': self.readout_register})
                        for qid in self.readout_register:
                            circ.append({'name': 'read', 'qubit': [qid]})
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
