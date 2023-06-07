import pdb

import matplotlib.pyplot as plt
import numpy as np
from abc import ABC
from qubic.job_manager_jpm import JobManager
from pygsti.circuits import Circuit
from collections import OrderedDict
from abstract_calibration import AbstractCalibrationExperiment


class XGateRepetition(AbstractCalibrationExperiment):
    """
    Gang's X-gate repetition trick
    """

    def __init__(self, qubits, delta_amp_range, n_full_rotations=1):
        if len(qubits) > 1:
            raise ValueError("XGateRepetition can only target 1 qubit")
        self.qubits = qubits

        self.delta_amp_range = delta_amp_range
        self.n_full_rotations = n_full_rotations

        self.optimization_parameters = [f'{qubits}X90.amp']
        self._results = None

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment

        """
        self.circuits = self._make_circuits(qchip)

        shots = jobmanager.collect_classified_shots(self.circuits, num_shots_per_circuit, qchip)
        self.ones_population = np.average(shots[self.qubits[0]], axis=1).flatten()

    def _fit_data(self, ones_pop, fit_routine=None, prior_estimates=None):
        """
        fit ones_pop vs amp to parabola
        """
        fit_params, cov = np.polyfit(self.amplitudes, ones_pop, deg=2)
        argmax = -fit_params[1]/(2*fit_params[2])
        assert self.amplitudes[0] < argmax < self.amplitudes[1]
        self._results = {'opt_amplitude': argmax,
                         'quadratic_fit_params': fit_params,
                         'fit_cov': cov}
        self.opt_amplitude = argmax

    def _make_circuits(self, qchip):
        """
        makes the circuits,
        circuits consist of initial delay followed by
        4n+2 X-pi/2 rotations
        and measurement of all the readout register
        """
        circuits = []
        self.amplitudes = self.delta_amp_range + qchip.gates['{}X90'.format(self.qubits)].contents[0].amp
        for amp in self.amplitudes:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6})
            for i in range(self.n_full_rotations):
                for _ in range(4*self.n_full_rotations+2):
                    cur_circ.append({'name': 'X90', 'qubit': self.qubits, 'modi': {(0, 'amp'): amp}})

                cur_circ.append({'name': 'read', 'qubit': self.qubits}) 

            circuits.append(cur_circ)

        return circuits

    @property
    def results(self):
        return self._results

    def update_qchip(self, qchip):
        qchip.gates['{}X90'.format(self.qubits[0])].contents[0].amp = self.opt_amplitude

    def update_gmm_manager(self, gmm_manager):
        pass
