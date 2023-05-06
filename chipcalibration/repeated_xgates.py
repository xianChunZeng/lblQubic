import abc
import pdb

import matplotlib.pyplot as plt
import numpy as np
from abc import ABC
from qubic.job_manager_jpm import JobManager
from pygsti.circuits import Circuit
from collections import OrderedDict



class AbstractCalibrationExperiment(ABC):
    """
    Abstract calibration class

    These classes are design to be used once as follows:
    1) the class is defined with any required arguments
    2) the class constructs circuits upon initialization
    3) the circuits are run using a provided jobmanager
    4) the results of the circuits are processed and analyzed and graphed
    5) changes to the qchip file can be made after the experiment

    =====================
    Public functions:
    =====================
    run_and_report(num_shots_per_circuit):
        runs circuits, analyzes results, and reports and or plots results

    update_chip(qchip):
        writes the parameters that the experiment found to be optimal to the qchip file

    optimization_parameters():
        returns a list of the labels of the parameters that the class optimizes

    =====================
    Member functions:
    =====================
    _make_circuits():
        returns a list of circuits that define the experiment

    _collect_data(jobmanager, num_shots_per_circuit, data_format):
        runs the circuits and returns the data in the data_format

    _fit_data(fit_routine, prior_params):
        fits the data according to the fit routine with the given priors
    """

    @abc.abstractmethod
    def __init__(self, target_register, readout_register, experimental_parameters):
        self.target_register = target_register
        self.readout_register = readout_register

        self.params = experimental_parameters # pulse widths, amplitudes, frequencies, etc...

        self.circuits = self._make_circuits()
        self.optimization_parameters = ['param_label_1', 'param_label_2', '...']
        self.final_estimated_params = None

    @abc.abstractmethod
    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment

        """
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        fit = self._fit_data(data)
        #plt.plot(fit)
        #make report and save

        self.final_estimated_params = 1 # find the estimate based on the experiment type, the data, and the fit
        return self.final_estimated_params

    @abc.abstractmethod
    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        pass



    @abc.abstractmethod
    def _make_circuits(self):
        """
        makes the circuits,
        all the default circuit parameters are stored in qchip
        and any changes are stored as class properties set at initialization
        note that the qchip is not stored in a calibration experiment,
        it is managed by the jobmanager, and a calibration class can update its parameters
        """
        pass

    @abc.abstractmethod
    def _collect_data(self, jobmanager, num_shots_per_circuit, qchip):
        """
        runs the circuits using the jabmanager
        the GMMM and the FPGA/Channel configs and the qchip is managed
        """
        pass


class XGateRepetition(AbstractCalibrationExperiment):
    """
    Gang's X-gate repetition trick
    """

    def __init__(self, target_register, readout_register, center_amp, sigma_amp, num_repetitions=10, num_partitions=20):
        if len(target_register) > 1:
            raise ValueError("XGateRepetition can only target 1 qubit")
        self.target_register = target_register
        self.readout_register = readout_register

        self.amplitudes = np.linspace(center_amp-sigma_amp, center_amp+sigma_amp, num_partitions)
        self.num_repetitions = num_repetitions

        self.circuits = self._make_circuits()

        self.optimization_parameters = [f'{target_register}X90.amp']
        self.final_estimated_params = None


    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment

        """
        data = self._collect_data(jobmanager, num_shots_per_circuit, qchip)
        self.data = data
        outcome_frequencies = [{
            qid : np.average(data[i][qid], axis=1) for qid in self.readout_register
        } for i in range(self.num_repetitions)]
        #fit = self._fit_data(data)
        #plt.plot(fit)
        #make report and save

        fig, axs = plt.subplots(len(self.readout_register))
        if len(self.readout_register) == 1:
            for i in range(self.num_repetitions):
                axs.plot(self.amplitudes, outcome_frequencies[i][self.readout_register[0]][:, 0])
        else:
            for idx, qid in enumerate(self.readout_register):
                for i in range(self.num_repetitions):
                    axs[idx].plot(self.amplitudes, outcome_frequencies[i][qid][:, 0])
        fig.show()

        self.final_estimated_params = 1 # find the estimate based on the experiment type, the data, and the fit
        return self.final_estimated_params

    def _fit_data(self, data, fit_routine=None, prior_estimates=None):
        pass



    def _make_circuits(self):
        """
        makes the circuits,
        circuits consist of initial delay followed by
        4n+2 X-pi/2 rotations
        and measurement of all the readout register
        """
        circuits = [OrderedDict()]*self.num_repetitions
        for i in range(self.num_repetitions):
            for amp in self.amplitudes:
                cur_circ = []
                cur_circ.append({'name': 'delay', 't': 400.e-6})
                for _ in range(4*i+2):
                    cur_circ.append({'name': 'X90', 'qubit': self.target_register, 'modi': {(0, 'amp'): amp}})
                cur_circ.append({'name': 'barrier', 'qubit': self.readout_register})
                for qubit in self.readout_register:
                    cur_circ.append({'name': 'read', 'qubit': [qubit]})
                pygsti_circuit = Circuit([('Gxpi2', self.target_register[0])]*(4*i+1))
                circuits[i][pygsti_circuit] = cur_circ
        return circuits

    def _collect_data(self, jobmanager: JobManager, num_shots_per_circuit, qchip):
        """
        """
        data = []
        for i in range(self.num_repetitions):
            data.append(jobmanager.collect_classified_shots(self.circuits[i], num_shots_per_circuit, qchip))
        return data
