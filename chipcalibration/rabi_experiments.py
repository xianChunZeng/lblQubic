import matplotlib.pyplot as plt
import numpy as np

from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
from qubic.state_disc import GMMManager
import warnings
from scipy.optimize import curve_fit

from collections import OrderedDict

class GMMRabi(AbstractCalibrationExperiment):
    """
    Time Rabi experiments to construct GMM Managers for the register
    """
    def __init__(self, target_qubits, target_amplitude, pulse_width_interval,
                 channel_configs):
        self.target_register = target_qubits
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in self.target_register}

        self.drive_amplitude = target_amplitude
        self.pulse_widths = pulse_width_interval

        self.circuits = self._make_circuits()

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip, plotting=True):
        """
        run the time Rabi experiments and make a GMM Manager
        """
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        self.raw_iq_shots = data
        self.gmm_manager.fit(data)
        self.shots = self.gmm_manager.predict(data)

        if plotting:
            fig, axs = plt.subplots(len(self.target_register))
            if len(self.target_register) == 1:
                axs = [axs]
            for idx, qid in enumerate(self.target_register):
                axs[idx].set_title(f'{qid}')
                rchan = str(self.readout_chanmap[qid])
                axs[idx].scatter(data[rchan].real, data[rchan].imag, c=self.shots[qid], cmap='viridis')
            plt.tight_layout()
            plt.show()

           
    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        pass

    def _make_circuits(self):
        """
        Make list of circuits used for rabi measurement. and the list of pulse width. So there will be a total of
        1 circuits, each of which contains len(pulse_widths) measurements. A 400 us
        delay is inserted between each measurement.
        """
        circuits = []
        for qid in self.target_register:
            for twidth in self.pulse_widths:
                cur_circ = []
                cur_circ.append({'name': 'delay', 't': 400.e-6})
                cur_circ.append({'name': 'rabi', 'qubit': [qid],
                                 'modi': {(0, 'twidth'): twidth}, (0, 'amp'): self.drive_amplitude})
                cur_circ.append({'name': 'barrier', 'qubit': self.target_register})
                for qid_read in self.target_register:
                    cur_circ.append({'name': 'read', 'qubit': [qid_read]})
                circuits.append(cur_circ)
        return circuits

    def _collect_data(self, jobmanager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        returns raw IQ shots
        """
        return jobmanager.collect_raw_IQ(self.circuits, num_shots_per_circuit, qchip)

class TimeRabi(AbstractCalibrationExperiment):
    """
    Time Rabi experiment based off standard Experiment class
    """

    def __init__(self, target_qubit, readout_register, target_amplitude, pulse_width_interval, chanmap_or_channel_configs=None):
        if chanmap_or_channel_configs is None:
            warnings.warn("""This class does not have a chanmap_or_channelconfigs, and so cannot fit GMMManagers
                          please set the chanmap with set_channel_info() before running""")
        if type(target_qubit) is not list:
            target_qubit = [target_qubit]
        if len(target_qubit) > 1:
            raise ValueError("TimeRabi can only target 1 qubit at a time")
        self.target_register = target_qubit
        self.readout_register = readout_register
        self.chanmap_or_channel_configs = chanmap_or_channel_configs
        self.gmm_manager = None
        self.drive_amplitude = target_amplitude

        self.pulse_widths = pulse_width_interval

        self.circuits = self._make_circuits()
        self.optimization_parameters = ['X90.twidth', 'X90.amp']
        self.final_estimated_params = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip, plotting=True):
        """
        run the time Rabi experiment

        will also create a GMM Manager along the way
        """
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        self.gmm_manager = self._fit_gmmm(data)
        self.shots = self.gmm_manager.predict(data)
        fits = self._fit_data()

        if plotting:
            fig, axs = plt.subplots(len(self.readout_register))
            if len(self.readout_register) == 1:
                axs = [axs]
            for idx, qid in enumerate(self.readout_register):
                axs[idx].plot(self.pulse_widths, np.average(self.shots[qid], axis=1))
                axs[idx].plot(self.pulse_widths, self._cos(self.pulse_widths, *fits[qid][0]), c='red')
            plt.tight_layout()
            plt.show()

        self.final_estimated_params = 1 # find the estimate based on the experiment type, the data, and the fit
        return self.final_estimated_params

    def set_channel_info(self, chanmap_or_channel_configs):
        self.chanmap_or_channel_configs = chanmap_or_channel_configs

    def _cos(self, x, A, B, drive_period, phi):
        return A*np.cos(2*np.pi*x/drive_period - phi) + B
    def _fit_data(self, data, fit_routine='fft', prior_fit_params=None):
        """
        fit the count data to a cosine
        """
        fits = dict()
        if prior_fit_params is None:
            prior_fit_params = [-0.5, 0.5, 1, 0]

        for qid in self.target_register:
            average_response = np.array([np.average(self.shots[qid][qid][i]) for i in range(len(self.pulse_widths))])
            if fit_routine == 'ff':
                try:
                        # this is "frequency" in terms of the rabi amplitude oscillation period
                        freq_ind_max = np.argmax(np.abs(np.fft.rfft(average_response)[1:])) + 1
                        freq_max = np.fft.rfftfreq(len(average_response), np.diff(self.pulse_widths)[0])[freq_ind_max]
                        prior_fit_params[qid][2] = 1 / freq_max
                        fits[qid] = curve_fit(self._cos, self.pulse_widths, average_response, prior_fit_params[qid])
                except:
                    print(f'Could not fit {qid}')
        return fits

    def _fit_gmmm(self, data):
        gmmm = GMMManager(chanmap_or_chan_cfgs=self.chanmap_or_channel_configs)
        if data is not None:
            gmmm.fit(data)
        else:
            raise ValueError("Must collect data before fitting a GMM Manager")
        return gmmm
    def _make_circuits(self):
        """
        Make list of circuits used for rabi measurement. and the list of pulse width. So there will be a total of
        1 circuits, each of which contains len(pulse_widths) measurements. A 400 us
        delay is inserted between each measurement.
        """
        circuits = []
        for twidth in self.pulse_widths:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6})
            cur_circ.append({'name': 'rabi', 'qubit': self.target_register,
                             'modi': {(0, 'twidth'): twidth}, (0, 'amp'): self.drive_amplitude})
            cur_circ.append({'name': 'barrier', 'qubit': self.readout_register})
            for qid in self.readout_register:
                cur_circ.append({'name': 'read', 'qubit': [qid]})
            circuits.append(cur_circ)
        return circuits


    def _collect_data(self, jobmanager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        returns raw IQ shots
        """
        return jobmanager.collect_raw_IQ(self.circuits, num_shots_per_circuit, qchip)
