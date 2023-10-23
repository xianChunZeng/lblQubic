import matplotlib.pyplot as plt
import numpy as np
import pdb
from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
from qubic.state_disc import GMMManager
import warnings
from scipy.optimize import curve_fit
import ipdb

from collections import OrderedDict

class GMMRabi(AbstractCalibrationExperiment):
    """
    Time Rabi experiments to construct GMM Managers for the register
    """
    def __init__(self, target_qubits, pulse_width_interval,
                 channel_configs, drive_amplitude=None,rabigate='rabi'):
        self.target_register = target_qubits
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in self.target_register}

        self.drive_amplitude = drive_amplitude
        self.pulse_widths = pulse_width_interval
        self.rabigate=rabigate

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip, plotting=True):
        """
        run the time Rabi experiments and make a GMM Manager
        """
        self.circuits = self._make_circuits(qchip)
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip=qchip)
        self.raw_iq_shots = data.copy()
        self.gmm_manager.fit(data)
        self.gmm_manager.set_labels_maxtomin(self.first_batch, [0, 1])
        self.shots = self.gmm_manager.predict(data)

        if plotting:
            fig, axs = plt.subplots(len(self.target_register))
            if len(self.target_register) == 1:
                axs = [axs]
            for idx, qid in enumerate(self.target_register):
                axs[idx].set_title(f'{qid}')
                rchan = str(self.readout_chanmap[qid])
                axs[idx].scatter(data[rchan].real, data[rchan].imag, c=self.shots[qid], cmap='viridis')
            #plt.tight_layout()
            plt.show()
            fig, axs = plt.subplots(len(self.target_register))
            if len(self.target_register) == 1:
                axs = [axs]
            for idx, qid in enumerate(self.target_register):
                axs[idx].set_title(f'{qid}')
                rchan = str(self.readout_chanmap[qid])
                axs[idx].hexbin(data[rchan].real, data[rchan].imag)
            #plt.tight_layout()
            plt.show()

           
    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        pass

    def _make_circuits(self, qchip):
        """
        Make list of circuits used for rabi measurement. and the list of pulse width. So there will be a total of
        1 circuits, each of which contains len(pulse_widths) measurements. A 400 us
        delay is inserted between each measurement.
        """
        circuits = []
        for qid in self.target_register:
            if self.drive_amplitude is None:
                drive_amplitude = qchip.gates['{}'.format(qid)+self.rabigate].contents[0].amp
            else:
                drive_amplitude = self.drive_amplitude
            for twidth in self.pulse_widths:
                cur_circ = []
                cur_circ.append({'name': 'delay', 't': 400.e-6})
                if twidth == 0: 
                    pass
                else: 
                    cur_circ.append({'name': self.rabigate, 'qubit': [qid],
                                 'modi': {(0, 'twidth'): twidth, (0, 'amp'): drive_amplitude}})
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
        self.first_batch = jobmanager.collect_raw_IQ(self.circuits, num_shots_per_circuit, qchip=qchip)
        return jobmanager.collect_raw_IQ(self.circuits, num_shots_per_circuit, qchip=qchip)

    @property
    def results(self):
        return None

    def plot_results(self, fig):
        pass

    def update_gmm_manager(self, gmm_manager):
        for qubit in self.target_register:
            gmm_manager.gmm_dict[qubit] = self.gmm_manager.gmm_dict[qubit]

    def update_qchip(self, qchip):
        pass

class TimeRabi(AbstractCalibrationExperiment):
    """
    Time Rabi experiment based off standard Experiment class
    """

    def __init__(self, target_qubit, readout_register, pulse_width_interval, drive_amplitude=None, rabigate='rabi'):
        if type(target_qubit) is not list:
            target_qubit = [target_qubit]
        if len(target_qubit) > 1:
            raise ValueError("TimeRabi can only target 1 qubit at a time")
        self.target_register = target_qubit
        self.readout_register = readout_register
        self.drive_amplitude = drive_amplitude

        self.pulse_widths = pulse_width_interval

        self.rabigate=rabigate
        self.optimization_parameters = ['X90.twidth', 'X90.amp']
        self.final_estimated_params = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip, plotting=True, fit_type='fft', period=None):
        """
        run the time Rabi experiment

        will also create a GMM Manager along the way
        """
        self.circuits = self._make_circuits(qchip=qchip)
        self.shots = self._collect_data(jobmanager, num_shots_per_circuit, qchip=qchip)
        self.raw_IQ = jobmanager.collect_raw_IQ(self.circuits, num_shots_per_circuit, qchip=qchip)
        fit = self._fit_data(self.shots[self.target_register[0]], fit_type, period)
        self.fitted_rabi_period = fit[0][2]
        
        if plotting:
            fig, axs = plt.subplots(len(self.readout_register))
            if len(self.readout_register) == 1:
                axs = [axs]
            for idx, qid in enumerate(self.readout_register):
                axs[idx].plot(self.pulse_widths, np.average(self.shots[qid], axis=1))
                axs[idx].set_title(qid)
                if qid == self.target_register[0]:
                    axs[idx].plot(self.pulse_widths, self._cos(self.pulse_widths, *fit[0]), c='red')
            plt.tight_layout()
            plt.show()
        self.final_estimated_params = [self.fitted_rabi_period/4, self.drive_amplitude] # find the estimate based on the experiment type, the data, and the fit
        return self.final_estimated_params

    def set_channel_info(self, chanmap_or_channel_configs):
        self.chanmap_or_channel_configs = chanmap_or_channel_configs

    def _cos(self, x, A, B, drive_period, phi):
        return A*np.cos(2*np.pi*x/drive_period - phi) + B
    
    def _fit_data(self, data, fit_routine='fft', period=None):
        """
        fit the count data to a cosine
        """
        prior_fit_params = [1, 0.5, 100.e-9, 0]

        average_response = np.squeeze(np.average(data, axis=1))
        if fit_routine == 'fft':
            try:
                # this is "frequency" in terms of the rabi amplitude oscillation period
                freq_ind_max = np.argmax(np.abs(np.fft.rfft(average_response)[1:])) + 1
                freq_max = np.fft.rfftfreq(len(average_response), np.diff(self.pulse_widths)[0])[freq_ind_max]
                prior_fit_params[2] = 1 / freq_max
                #ipdb.set_trace()
                fit = curve_fit(self._cos, self.pulse_widths, average_response, prior_fit_params)
                return fit
            except Exception as e:
                print('Could not fit with FFT, try again with a user defined period')
                raise e
                
        if fit_routine == 'period':
            if period is None:
                fig, axs = plt.subplots()
                axs.plot(self.pulse_widths, np.average(self.shots[self.target_register[0]], axis=1))
                axs.set_title(self.target_register[0])
                plt.show()
                period = float(input("Please input the approximate period of the oscillation for fitting"))
                prior_fit_params[2] = period
            return curve_fit(self._cos, self.pulse_widths, average_response, prior_fit_params)
        
    def update_qchip(self, qchip,x90gate='X90'):
        if self.final_estimated_params is None:
            raise ValueError('Please run a Rabi experiment before updating the qchip')
        qchip.gates[self.target_register[0]+x90gate].contents[0].amp = self.final_estimated_params[1]
        qchip.gates[self.target_register[0]+x90gate].contents[0].twidth = self.final_estimated_params[0]
        
        

    def _make_circuits(self, qchip):
        """
        Make list of circuits used for rabi measurement. and the list of pulse width. So there will be a total of
        1 circuits, each of which contains len(pulse_widths) measurements. A 400 us
        delay is inserted between each measurement.
        """
        circuits = []
        for twidth in self.pulse_widths:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6})
            if twidth == 0: 
                pass
            else: 
                if self.drive_amplitude is not None:
                    cur_circ.append({'name': self.rabigate, 'qubit': self.target_register,
                                    'modi': {(0, 'twidth'): twidth, (0, 'amp'): self.drive_amplitude}})
                else:
                    cur_circ.append({'name': self.rabigate, 'qubit': self.target_register,
                                    'modi': {(0, 'twidth'): twidth}})
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
        return jobmanager.collect_classified_shots(self.circuits, num_shots_per_circuit, qchip=qchip)

    @property
    def results(self):
        return None

    def plot_results(self):
        pass

    def update_gmm_manager(self, gmm_manager):
        pass

    def update_qchip(self, qchip):
        pass
#==========================================================================================
    
class AmpRabi(AbstractCalibrationExperiment):
    """
    Amplitude Rabi experiment based off standard Experiment class
    """

    def __init__(self, target_qubit, readout_register, amp_interval, target_pulse_width,rabigate='rabi'):
        if type(target_qubit) is not list:
            target_qubit = [target_qubit]
        if len(target_qubit) > 1:
            raise ValueError("TimeRabi can only target 1 qubit at a time")
        self.target_register = target_qubit
        self.readout_register = readout_register
        self.amp_interval = amp_interval

        self.target_pulse_width = target_pulse_width

        self.rabigate=rabigate
        self.circuits = self._make_circuits()
        self.optimization_parameters = ['X90.twidth', 'X90.amp']
        self.final_estimated_params = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip, plotting=True, fit_type='fft', period=None):
        """
        run the time Rabi experiment

        will also create a GMM Manager along the way
        """
        self.shots = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        fit = self._fit_data(self.shots[self.target_register[0]], fit_type, period)
        self.fitted_rabi_period = fit[0][2]
        
        if plotting:
            fig, axs = plt.subplots(len(self.readout_register))
            if len(self.readout_register) == 1:
                axs = [axs]
            for idx, qid in enumerate(self.readout_register):
                axs[idx].plot(self.amp_interval, np.average(self.shots[qid], axis=1))
                axs[idx].set_title(qid)
                if qid == self.target_register[0]:
                    axs[idx].plot(self.amp_interval, self._cos(self.amp_interval, *fit[0]), c='red')
            plt.tight_layout()
            plt.show()
        #self.final_estimated_params = [self.fitted_rabi_period/4, self.drive_amplitude] # find the estimate based on the experiment type, the data, and the fit
        return self.final_estimated_params

    def set_channel_info(self, chanmap_or_channel_configs):
        self.chanmap_or_channel_configs = chanmap_or_channel_configs

    def _cos(self, x, A, B, drive_period, phi):
        return A*np.cos(2*np.pi*x/drive_period - phi) + B
    
    def _fit_data(self, data, fit_routine='fft', period=None):
        """
        fit the count data to a cosine
        """
        prior_fit_params = [-0.5, 0.5, period, 0]

        average_response = np.average(data, axis=1)
        if fit_routine == 'fft':
            try:
                # this is "frequency" in terms of the rabi amplitude oscillation period
                freq_ind_max = np.argmax(np.abs(np.fft.rfft(average_response)[1:])) + 1
                freq_max = np.fft.rfftfreq(len(average_response), np.diff(self.pulse_widths)[0])[freq_ind_max]
                prior_fit_params[2] = 1 / freq_max
                fit = curve_fit(self._cos, self.amp_interval, average_response[:, 0], prior_fit_params)
                return fit
            except:
                print('Could not fit with FFT, try again with a user defined period')
        if fit_routine == 'period':
            if period is None:
                fig, axs = plt.subplots()
                axs.plot(self.amp_interval, np.average(self.shots[self.target_register[0]], axis=1))
                axs.set_title(self.target_register[0])
                plt.show()
                period = float(input("Please input the approximate period of the oscillation for fitting"))
                prior_fit_params[2] = period
            return curve_fit(self._cos, self.amp_interval, average_response[:, 0], prior_fit_params)
        
    def update_qchip(self, qchip,x90gate='X90'):
        if self.final_estimated_params is None:
            raise ValueError('Please run a Rabi experiment before updating the qchip')
        qchip.gates[self.target_register[0]+x90gate].contents[0].amp = self.final_estimated_params[1]
        qchip.gates[self.target_register[0]+x90gate].contents[0].twidth = self.final_estimated_params[0]
        
        

    def _make_circuits(self):
        """
        Make list of circuits used for rabi measurement. and the list of pulse width. So there will be a total of
        1 circuits, each of which contains len(pulse_widths) measurements. A 400 us
        delay is inserted between each measurement.
        """
        circuits = []
        for amp in self.amp_interval:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6})
            if amp == 0:
                pass
            else:
                cur_circ.append({'name': self.rabigate, 'qubit': self.target_register,
                             'modi': {(0, 'twidth'): self.target_pulse_width}, (0, 'amp'): amp})
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
        return jobmanager.collect_classified_shots(self.circuits, num_shots_per_circuit, qchip=qchip)

    @property
    def results(self):
        return None

    def plot_results(self):
        pass

    def update_gmm_manager(self, gmm_manager):
        pass

