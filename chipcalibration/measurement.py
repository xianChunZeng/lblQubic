import matplotlib.pyplot as plt

from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
import numpy as np
from collections import OrderedDict
from qubic.job_manager import JobManager
from qubic.state_disc import GMMManager

class MeasurementCalibration(AbstractCalibrationExperiment):
    """
    Measurement calibration -- can target one qubit

    uses a jobmanager without a GMM Manager
    """

    def __init__(self, readout_qid, delta_amps_db, delta_freqs, channel_config):
        if isinstance(readout_qid, list):
            if len(readout_qid) > 1:
                raise ValueError("Can only target one qubit at a time")
            else:
                self.target_register = readout_qid
        else:
            self.target_register = [readout_qid]

        self.delta_amps_db = delta_amps_db
        self.delta_freqs = delta_freqs
        self.channel_config = channel_config

        self.optimization_parameters = ['MEAS.amp', 'MEAS.freq']
        #self.final_estimated_params = None
        self._results = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment

        makes a color mesh grid of the GMM separation as a function of amp and freq
        """
        self.circuits = self._make_circuits(qchip)
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        separations, gmmms = self._fit_data(data)
        self.fitted_gmmms = gmmms
        self.separations = separations
        opt_amp_ind, opt_freq_ind = np.unravel_index(np.argmax(separations), separations.shape)
        self.opt_amp = self.amps[opt_amp_ind]
        self.opt_freq = self.freqs[opt_freq_ind]
        self._results = {'opt_amp': self.opt_amp, 'opt_freq': self.opt_freq}
        

    def plot_results(self, fig):
        ax = fig.add_subplot(111)
        colormesh = ax.pcolormesh(self.amps, self.freqs, self.separations.T, cmap='cividis')
        ax.colorbar(colormesh)
        fig.title("Sigma separation as a function of amp and freq")
        ax.set_xlabel("amp")
        ax.set_ylabel("freq")
        ax.set_xlim([self.amps[0], self.amps[-1]])
        ax.set_ylim([self.freqs[0], self.freqs[-1]])
        plt.show()


    @staticmethod
    def sigma_separation(gmm):
        mean_i = gmm.means_[0]
        covar_i = gmm.covariances_[0]
        mean_j = gmm.means_[1]
        covar_j = gmm.covariances_[1]
        return np.sqrt((mean_i - mean_j) @ np.linalg.inv((covar_i + covar_j) / 2) @ (mean_i - mean_j).T)

    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        separations = np.zeros((len(self.amps), len(self.freqs)))
        gmmms = OrderedDict()
        for id_amp, amp, in enumerate(self.amps):
            for id_freq, freq in enumerate(self.freqs):
                gmm_manager = GMMManager(chanmap_or_chan_cfgs=self.channel_config)
                gmm_manager.fit(data[(id_amp, id_freq)])
                gmmms[(id_amp, id_freq)] = gmm_manager
                separations[(id_amp, id_freq)] = self.sigma_separation(gmm_manager.gmm_dict[self.target_register[0]].gmmfit)
        return separations, gmmms

    @property
    def results(self):
        return self._results

    def _make_circuits(self, qchip):
        """
        makes the circuits,


        all the default circuit parameters are stored in qchip
        and any changes are stored as class properties set at initialization
        note that the qchip is not stored in a calibration experiment,
        it is managed by the jobmanager, and a calibration class can update its parameters
        """
        circuits = OrderedDict()
        self.amps = 10**(self.delta_amps_db/20)*qchip.gates[self.readgate_name].contents[0].amp
        self.freqs = self.delta_freqs + qchip.qubits[self.target_register[0]].readfreq
        for id_amp, amp, in enumerate(self.amps):
            for id_freq, freq in enumerate(self.freqs):
                circ_instruction = [
                    {'name': 'delay', 't': 400.e-6},
                    {'name': 'X90', 'qubit': self.target_register},
                    {'name': 'read', 'qubit': self.target_register,
                     'modi': {(0, 'amp'): amp, (0, 'freq'):
                        freq, (1, 'freq'): freq}},
                ]
                circuits[(id_amp, id_freq)] = circ_instruction
        return circuits

    def _collect_data(self, jobmanager: JobManager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        the GMMM and the FPGA/Channel configs and the qchip is managed
        """
        data = OrderedDict()
        amp_range = range(len(self.amps))
        freq_range = range(len(self.freqs))
        
        for id_amp, amp, in enumerate(self.amps):
            for id_freq, freq in enumerate(self.freqs):
                # TODO: convert to a single batch using np.ravel or some such
                data[(id_amp, id_freq)] = jobmanager.collect_raw_IQ([self.circuits[(id_amp, id_freq)]], num_shots_per_circuit, qchip=qchip)
        return data

    def update_gmm_manager(self, gmm_manager):
        pass

    @property
    def readgate_name(self):
        return '{}read'.format(self.target_register[0])

    def update_qchip(self, qchip):
        qchip.qubits[self.target_register[0]].readfreq = self.opt_freq
        qchip.gates[self.readgate_name].contents[0].amp = self.opt_amp
