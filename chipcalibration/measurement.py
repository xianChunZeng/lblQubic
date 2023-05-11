import matplotlib.pyplot as plt

from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
import numpy as np
from collections import OrderedDict
from qubic.job_manager_jpm import JobManager
from qubic.state_disc import GMMManager

class MeasurementCalibration(AbstractCalibrationExperiment):
    """
    Measurement calibration -- can target one qubit

    uses a jobmanager without a GMM Manager
    """

    def __init__(self, readout_qid, amp_interval, freq_interval, channel_config):
        if isinstance(readout_qid, list):
            if len(readout_qid) > 1:
                raise ValueError("Can only target one qubit at a time")
            else:
                self.target_register = readout_qid
        else:
            self.target_register = [readout_qid]


        self.amp_interval = amp_interval
        self.freq_interval = freq_interval
        self.channel_config = channel_config

        self.circuits = self._make_circuits()
        self.optimization_parameters = ['MEAS.amp', 'MEAS.freq']
        self.final_estimated_params = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip, center_amp=None, center_freq=None):
        """
        run the experiment

        makes a color mesh grid of the GMM separation as a function of amp and freq
        """
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        separations, gmmms = self._fit_data(data)
        self.fitted_gmmms = gmmms
        self.separations = separations
        colormesh = plt.pcolormesh(self.amp_interval, self.freq_interval, separations.T, cmap='cividis')
        plt.colorbar(colormesh)
        if center_amp is not None and center_freq is not None:
            plt.scatter(center_amp, center_freq, marker='X', s=400, c='red')
        plt.title("Sigma separation as a function of amp and freq")
        plt.xlabel("amp")
        plt.ylabel("freq")
        plt.xlim([self.amp_interval[0], self.amp_interval[-1]])
        plt.ylim([self.freq_interval[0], self.freq_interval[-1]])
        plt.show()


    @staticmethod
    def sigma_separation(gmm):
        mean_i = gmm.means_[0]
        covar_i = gmm.covariances_[0]
        mean_j = gmm.means_[1]
        covar_j = gmm.covariances_[1]
        return np.sqrt((mean_i - mean_j) @ np.linalg.inv((covar_i + covar_j) / 2) @ (mean_i - mean_j).T)

    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        separations = np.zeros((len(self.amp_interval), len(self.freq_interval)))
        gmmms = OrderedDict()
        for id_amp, amp, in enumerate(self.amp_interval):
            for id_freq, freq in enumerate(self.freq_interval):
                gmm_manager = GMMManager(chanmap_or_chan_cfgs=self.channel_config)
                gmm_manager.fit(data[(id_amp, id_freq)])
                gmmms[(id_amp, id_freq)] = gmm_manager
                separations[(id_amp, id_freq)] = self.sigma_separation(gmm_manager.gmm_dict[self.target_register[0]].gmmfit)
        return separations, gmmms




    def _make_circuits(self):
        """
        makes the circuits,


        all the default circuit parameters are stored in qchip
        and any changes are stored as class properties set at initialization
        note that the qchip is not stored in a calibration experiment,
        it is managed by the jobmanager, and a calibration class can update its parameters
        """
        circuits = OrderedDict()
        for id_amp, amp, in enumerate(self.amp_interval):
            for id_freq, freq in enumerate(self.freq_interval):
                circ_instruction = [
                    {'name': 'delay', 't': 400.e-6},
                    {'name': 'X90', 'qubit': self.target_register},
                    {'name': 'read', 'qubit': self.target_register,
                     'modi': {(0, 'amp'): amp, (0, 'fcarrier'):
                        freq, (1, 'fcarrier'): freq}},
                ]
                circuits[(id_amp, id_freq)] = circ_instruction
        return circuits

    def _collect_data(self, jobmanager: JobManager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        the GMMM and the FPGA/Channel configs and the qchip is managed
        """
        data = OrderedDict()
        amp_range = range(len(self.amp_interval))
        freq_range = range(len(self.freq_interval))
        
        for id_amp, amp, in enumerate(self.amp_interval):
            for id_freq, freq in enumerate(self.freq_interval):
                # TODO: convert to a single batch using np.ravel or some such
                data[(id_amp, id_freq)] = jobmanager.collect_raw_IQ([self.circuits[(id_amp, id_freq)]], num_shots_per_circuit, qchip=qchip)
        return data