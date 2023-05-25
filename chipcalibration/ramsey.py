import matplotlib.pyplot as plt
import numpy as np
import pdb
from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
from qubic.state_disc import GMMManager
import warnings
from scipy.optimize import curve_fit
import copy
import logging
from chipcalibration.graph import CalibrationGraph
from chipcalibration.rabi_experiments import GMMRabi


from collections import OrderedDict

class RamseyExperiment(AbstractCalibrationExperiment):

    def __init__(self, target_register, readout_register, delay_interval, drive_frequency):
        if len(target_register) > 1:
            raise ValueError('This class only supports Ramsey experiments on a single target qubit')
        if type(target_register) is not list:
            target_register = [target_register]

        super().__init__(target_register, readout_register)
        self.delay_interval = delay_interval
        self.drive_frequency = drive_frequency
        self.circuits = self._make_ramsey_circuits()
        self.prior_fit_params = [0.720, 0.5, 500e3, 0, 1.5e-5]
        self._results = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        shots = jobmanager.collect_classified_shots(self.circuits, num_shots_per_circuit, qchip=qchip)
        observed_probabilities = np.average(shots[self.target_register[0]], axis=1).flatten()

        cos_exp_fit = self._fit_cos_exp(observed_probabilities, prior_estimates=self.prior_fit_params)
        exp_fit = self._fit_exp(observed_probabilities, prior_estimates=[self.prior_fit_params[0], 
                                                        self.prior_fit_params[1], self.prior_fit_params[-1]])

        self._results = {}
        if cos_exp_fit is not None:
            self._results['cos_exp'] = {'scale': cos_exp_fit[0][0], 'offset': cos_exp_fit[0][1], 
                                       'freq': cos_exp_fit[0][2], 'phase': cos_exp_fit[0][3], 'decay_time': cos_exp_fit[0][4],
                                        'total_var': np.diag(cos_exp_fit[1]), 
                                        'fit_resid': self._calc_fit_res_sq(self.delay_interval, observed_probabilities, cos_exp_fit[0], self._cos_exp)}
        else:
            self._results['cos_exp'] = None

        if exp_fit is not None:
            self._results['exp'] = {'scale': exp_fit[0][0], 'offset': exp_fit[0][1], 'decay_time': exp_fit[0][2],
                                    'total_var': np.diag(exp_fit[1]),
                                    'fit_resid': self._calc_fit_res_sq(self.delay_interval, observed_probabilities, exp_fit[0], self._exp)}
        else:
            self._results['exp'] = None

        self.cos_exp_fit = cos_exp_fit
        self.exp_fit = exp_fit
        self.observed_probabilities = observed_probabilities
        self.shots = shots

    @staticmethod
    def _calc_fit_res_sq(delay_interval, observed_probabilities, fit_params, fit_function):
        return np.sum((observed_probabilities.flatten() - fit_function(delay_interval, *fit_params))**2)


    def plot_results(self, fig):
        ax = fig.add_subplot(111)
        ax.plot(self.delay_interval, self.observed_probabilities, '.-', label='results')
        if self.cos_exp_fit is not None:
            ax.plot(self.delay_interval, self._cos_exp(self.delay_interval, *self.cos_exp_fit[0]), label='cos exp fit')
        if self.exp_fit is not None:
            ax.plot(self.delay_interval, self._exp(self.delay_interval, *self.exp_fit[0]), label='exp_fit')
        ax.set_xlabel('ramsey delay (s)')
        ax.set_ylabel('|1> state population')
        ax.legend()


    def _make_ramsey_circuits(self):
        circuits = []
        for dtime in self.delay_interval:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6, 'qubit': self.target_register})
            cur_circ.append(
                {'name': 'X90', 'qubit': self.target_register, 'modi': {(0, 'fcarrier'): self.drive_frequency}})
            cur_circ.append({'name': 'delay', 't': dtime, 'qubit': self.target_register})
            cur_circ.append(
                {'name': 'X90', 'qubit': self.target_register, 'modi': {(0, 'fcarrier'): self.drive_frequency}})
            cur_circ.append({'name': 'read', 'qubit': self.target_register})
            circuits.append(cur_circ)
        return circuits

    def _fit_cos_exp(self, observed_probabilities, use_initial_fft=True, prior_estimates=None):
        """
        cosine decaying exponentially with offset
        params are [A, B, drive_freq, phi, exp_decay]
        """
        self.fit_params = {}
        prior_fit_params = copy.deepcopy(prior_estimates)
        qid = self.target_register[0]
        if use_initial_fft:
            freq_ind_max = np.argmax(np.abs(np.fft.rfft(observed_probabilities)[1:])) + 1
            freq_max = np.fft.rfftfreq(len(self.delay_interval), 
                                       np.diff(self.delay_interval)[0])[min(freq_ind_max, len(observed_probabilities))]
            prior_fit_params[2] = freq_max
        try:
            fit = curve_fit(self._cos_exp, self.delay_interval[1:], observed_probabilities[1:].flatten(),
                                               prior_fit_params)
        except RuntimeError:
            fit = None
            logging.getLogger(__name__).warning('cos_exp could not be fit')
        return fit

    def _fit_exp(self, observed_probabilities, prior_estimates=None):
        """
        cosine decaying exponentially with offset
        params are [A, B, drive_freq, phi, exp_decay]
        """
        self.fit_params = {}
        prior_fit_params = copy.deepcopy(prior_estimates)
        qid = self.target_register[0]
        try:
            fit = curve_fit(self._exp, self.delay_interval[1:], observed_probabilities[1:].flatten(),
                                               prior_fit_params)
        except RuntimeError:
            fit = None
            logging.getLogger(__name__).warning('exp could not be fit')
        return fit

    @staticmethod
    def _cos_exp(x, scale, offset, drive_freq, phi, exp_decay):
        return scale*np.exp(-x/exp_decay)*np.cos(2*np.pi*x*drive_freq - phi) + offset

    @staticmethod
    def _exp(x, scale, offset, exp_decay):
        return scale*np.exp(-x/exp_decay) + offset

    @property
    def results(self):
        return self._results

    def update_gmm_manager(self, gmm_manager):
        pass

    def update_qchip(self, qchip):
        pass

class RamseyOptimize(AbstractCalibrationExperiment):

    def __init__(self, target_register, readout_register, delay_stepsize, n_delay_steps):
        self.delay_interval = np.arange(0, n_delay_steps*delay_stepsize, delay_stepsize)

        super().__init__(target_register, readout_register)

        self._results = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        self.initial_ramsey = RamseyExperiment(self.target_register, self.readout_register, 
                                               self.delay_interval, qchip.qubits[self.target_register[0]].freq)
        self.initial_ramsey.run_and_report(jobmanager, num_shots_per_circuit, qchip)

        self._results = {'initial': self.initial_ramsey.results}
        logging.getLogger(__name__).info('done initial ramsey; results: {}'.format(self._results))

        if self.initial_ramsey.results['exp'] is not None and self.initial_ramsey.results['cos_exp'] is not None:
            if self.initial_ramsey.results['exp']['fit_resid'] > self.initial_ramsey.results['cos_exp']['fit_resid']:
                self.delta_freq = self.initial_ramsey.results['cos_exp']['freq']
            else:
                self.delta_freq = 0
                logging.getLogger(__name__).info('found deltaf = 0')

        elif self.initial_ramsey.results['exp'] is None and self.initial_ramsey.results['cos_exp'] is None:
            raise RuntimeError('All initial Ramsey fits failed!')

        elif self.initial_ramsey.results['exp'] is None:
            self.delta_freq = self.initial_ramsey.results['cos_exp']['freq']

        else:
            self.delta_freq = 0
            self._results['delta_freq'] = 0
            logging.getLogger(__name__).info('found deltaf = 0')


        if self.delta_freq != 0:
            self.pos_ramsey = RamseyExperiment(self.target_register, self.readout_register, 
                                               self.delay_interval, qchip.qubits[self.target_register[0]].freq + self.delta_freq)
            self.neg_ramsey = RamseyExperiment(self.target_register, self.readout_register, 
                                               self.delay_interval, qchip.qubits[self.target_register[0]].freq - self.delta_freq)

            self.pos_ramsey.run_and_report(jobmanager, num_shots_per_circuit, qchip)
            logging.getLogger(__name__).info('done + ramsey; results: {}'.format(self.pos_ramsey.results))
            self.neg_ramsey.run_and_report(jobmanager, num_shots_per_circuit, qchip)
            logging.getLogger(__name__).info('done - ramsey; results: {}'.format(self.neg_ramsey.results))

            pos_ramsey_valid = self.pos_ramsey.results['exp'] is not None and (self.pos_ramsey.results['cos_exp'] is None
                                            or (self.pos_ramsey.results['exp']['fit_resid'] < self.pos_ramsey.results['cos_exp']['fit_resid'])\
                                            or (self.pos_ramsey.results['cos_exp']['freq'] < 50))
            
            neg_ramsey_valid = self.neg_ramsey.results['exp'] is not None and (self.neg_ramsey.results['cos_exp'] is None
                                            or (self.neg_ramsey.results['exp']['fit_resid'] < self.neg_ramsey.results['cos_exp']['fit_resid'])\
                                            or (self.neg_ramsey.results['cos_exp']['freq'] < 50))

            if pos_ramsey_valid and neg_ramsey_valid: #neither show oscillations, meaning something went wrong
                raise RuntimeError('+/- freq could not be determined!')

            if not (pos_ramsey_valid or neg_ramsey_valid): #both show significant oscillations
                self.delta_freq = 0

            if neg_ramsey_valid:
                self.delta_freq *= -1

            self._results['pos'] = self.pos_ramsey.results
            self._results['neg'] = self.neg_ramsey.results
            self._results['delta_freq'] = self.delta_freq

    @property
    def results(self):
        return self._results

    def update_gmm_manager(self, gmm_manager):
        pass

    def update_qchip(self, qchip):
        qchip.qubits[self.target_register[0]].freq += self.delta_freq

    def plot_results(self, fig):
        figs = fig.subfigures(3,1, hspace=1)
        self.initial_ramsey.plot_results(figs[0])
        self.neg_ramsey.plot_results(figs[1])
        self.pos_ramsey.plot_results(figs[2])

def get_refinement_graph(qubits, dt_steps=[2.e-9, 30.e-9, 500.e-9], n_delay_times=100, shots_per_circuit=1000):
    cal_graph = CalibrationGraph()
    for qubit in qubits:
        pred_nodes = None

        for i, step in enumerate(dt_steps):
            cal_obj = RamseyOptimize([qubit], [qubit], step, n_delay_times)
            step_name = '{}_ramsey_opt_{}'.format(qubit, i)
            cal_graph.add_calibration_step(step_name, cal_obj, [qubit], pred_nodes, shots_per_circuit)
            pred_nodes = ['{}_ramsey_opt_{}'.format(qubit, i)]

    return cal_graph


