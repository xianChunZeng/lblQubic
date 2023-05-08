import matplotlib.pyplot as plt
import numpy as np
import pdb
from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
from qubic.state_disc import GMMManager
import warnings
from scipy.optimize import curve_fit

from collections import OrderedDict

class RamseyExperiment(AbstractCalibrationExperiment):
    def __init__(self, target_register, readout_register, delay_interval, drive_frequency):
        if len(target_register) > 1:
            raise ValueError('This class only supports Ramsey experiments on a single target qubit')
        if type(target_register) is not list:
            target_register = [target_register]
        self.target_register = target_register
        self.readout_register = readout_register
        self.delay_interval = delay_interval
        self.initial_drive_frequency = drive_frequency


        self.optimization_parameters = ['Qid.freq']
        self.estimated_drive_frequency = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        initial_shots = self.run_ramsey(self.initial_drive_frequency, self.delay_interval, jobmanager, num_shots_per_circuit, qchip)
        try:
            fit = self._fit_data(initial_shots)
        except:
            print("Could not fit the initial reading")

        estimated_frequency_differential = fit[1]

        pos_shots = self.run_ramsey(self.initial_drive_frequency+estimated_frequency_differential, num_shots_per_circuit, qchip)
        neg_shots = self.run_ramsey(self.initial_drive_frequency-estimated_frequency_differential, num_shots_per_circuit, qchip)

        try:
            pos_fit = self._fit_data(pos_shots)
            neg_fit = self._fit_data(neg_shots)
        except:
            print("Could not fit one of the differential readings")

        pos_probs = np.average(pos_shots[self.target_register[0]], axis=1)
        neg_probs = np.average(neg_shots[self.target_register[0]], axis=1)
        plt.plot(self.delay_interval, pos_probs)
        plt.figure()
        plt.plot(self.delay_interval, neg_probs)

        pos_or_neg = input("Was the differential pos or neg? Please input +1 or -1")
        if pos_or_neg == +1:
            self.estimated_drive_frequency = self.initial_drive_frequency+estimated_frequency_differential
        elif pos_or_neg == -1:
            self.estimated_drive_frequency = self.initial_drive_frequency-estimated_frequency_differential

        return self.estimated_drive_frequency

    @staticmethod
    def _cos_exp(x, scale, offset, drive_freq, phi, exp_decay):
        return scale*np.exp(-x/exp_decay)*np.cos(2*np.pi*x*drive_freq - phi) + offset

    def run_ramsey(self, drive_frequency, delay_interval, jobmanager, num_shots_per_circuit, qchip):
        circuits = self.ramsey_circuits(drive_frequency, delay_interval)
        return jobmanager.collect_classified_shots(circuits, num_shots_per_circuit, qchip,
                                                    reads_per_shot=len(self.delay_interval),
                                                    delay_per_shot=(500.e-6*len(self.delay_interval))+sum(self.delay_interval) )

    def _fit_data(self, observed_probabilities, fit_routine=None, prior_estimates=None):
        """
        cosine decaying exponentially with offset
        params are [A, B, drive_freq, phi, exp_decay]
        """
        self.fit_params = {}
        prior_fit_params = copy.deepcopy(prior_fit_params)
        for qubit in self.qubits:
            if use_fft:
                freq_ind_max = np.argmax(np.abs(np.fft.rfft(self.ones_frac[qubit])[1:])) + 1
                freq_max = np.fft.rfftfreq(len(self.delaytime), np.diff(self.delaytime)[0])[freq_ind_max]
                prior_fit_params[qubit][2] = freq_max
            try:
                self.fit_params[qubit] = curve_fit(self._cos_exp, self.delaytime[1:], self.ones_frac[qubit][1:],
                                                   prior_fit_params[qubit])
            except RuntimeError:
                print('{} could not be fit')

    def ramsey_circuits(self, drive_freq, delay_interval):
        circuits = []
        for dtime in delay_interval:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6, 'qubit': self.target_register})
            cur_circ.append(
                {'name': 'X90', 'qubit': self.target_register, 'modi': {(0, 'fcarrier'): drive_freq}})
            cur_circ.append({'name': 'delay', 't': dtime, 'qubit': self.target_register})
            cur_circ.append(
                {'name': 'X90', 'qubit': self.target_register, 'modi': {(0, 'fcarrier'): drive_freq}})
            cur_circ.append({'name': 'read', 'qubit': self.target_register})
            circuits.append(cur_circ)
        return circuits

    def _make_circuits(self):
        """
        Make a circuit used for ramsey measurement with the list of delaytime. So there will be a total of
        1 circuits, each of which contains len(delaytime) measurements. A 400 us
        delay is inserted between each measurement.
        Limitations:
        1. Length of the delaytime list is restricted to 1024 because of memory
        2. Due to command buffer limitation we can have 341 of x90->delay->X90
        """
        return self.ramsey_circuits(self.drive_frequency, self.delay_interval)

    def _collect_data(self, jobmanager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        the GMMM and the FPGA/Channel configs and the qchip is managed
        """
        pass