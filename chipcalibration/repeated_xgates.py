import pdb

import matplotlib.pyplot as plt
import numpy as np
from abc import ABC
from collections import OrderedDict
from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
import logging


class XGateRepetition(AbstractCalibrationExperiment):
    """
    Gang's X-gate repetition trick
    """

    def __init__(self, qubits, delta_amp_frac, n_amps=30, n_full_rotations=1, gateX90='X90'):
        if len(qubits) > 1:
            raise ValueError("XGateRepetition can only target 1 qubit")
        self.qubits = qubits

        self.delta_amp_frac = delta_amp_frac
        self.n_amps = n_amps
        self.n_full_rotations = n_full_rotations
        self.gateX90=gateX90
        self.optimization_parameters = [f'{qubits}X90.amp']
        self._results = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment

        """
        self.circuits = self._make_circuits(qchip)

        shots = jobmanager.collect_classified_shots(self.circuits, num_shots_per_circuit, qchip=qchip)
        self.ones_population = np.average(shots[self.qubits[0]], axis=1).flatten()
        self.ones_population = self.ones_population.astype(np.float64)
        self._fit_data(self.ones_population)

    def _fit_data(self, ones_pop, fit_routine=None, prior_estimates=None):
        """
        fit ones_pop vs amp to parabola
        """
        fit_params, cov = np.polyfit(self.amplitudes, ones_pop, deg=2, cov=True)
        argmax = -fit_params[1]/(2*fit_params[0])
        assert self.amplitudes[0] < argmax < self.amplitudes[-1]
        self._results = {'opt_amplitude': argmax,
                         'quadratic_fit_params': fit_params,
                         'fit_cov': cov}
        self.opt_amplitude = argmax
        logging.getLogger(__name__).info('finished {} repeatgate. results: {}'.format(self.qubits[0], self._results))

    def plot_results(self, fig):
        ax = fig.add_subplot(111)
        ax.plot(self.amplitudes, self.ones_population)
        ax.axvline(self.opt_amplitude, 0, 1, color='r')
        ax.axvline(self.amplitudes[len(self.amplitudes)//2], 0, 1, color='r', ls=':')
        ax.set_xlabel('amplitude')
        ax.set_ylabel('|1> population')

    def _make_circuits(self, qchip):
        """
        makes the circuits,
        circuits consist of initial delay followed by
        4n+2 X-pi/2 rotations
        and measurement of all the readout register
        """
        circuits = []
        tmp1='{}'+self.gateX90
        qchip_amp = qchip.gates[tmp1.format(self.qubits[0])].contents[0].amp
        self.amplitudes = qchip_amp + np.linspace(-qchip_amp*self.delta_amp_frac, qchip_amp*self.delta_amp_frac, self.n_amps)
        self.amplitudes = self.amplitudes[self.amplitudes < 1] #clip/delete amps
        for amp in self.amplitudes:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6})
            for _ in range(4*self.n_full_rotations+2):
                cur_circ.append({'name': self.gateX90, 'qubit': self.qubits, 'modi': {(0, 'amp'): amp}})

            cur_circ.append({'name': 'read', 'qubit': self.qubits}) 

            circuits.append(cur_circ)

        return circuits

    @property
    def results(self):
        return self._results

    def update_qchip(self, qchip):
        tmp='{}'+self.gateX90
        qchip.gates[tmp.format(self.qubits[0])].contents[0].amp = self.opt_amplitude

    def update_gmm_manager(self, gmm_manager):
        pass
