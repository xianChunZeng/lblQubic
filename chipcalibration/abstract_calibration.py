import abc
from abc import ABC

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
    def __init__(self, target_register, readout_register):
        self.target_register = target_register
        self.readout_register = readout_register

    @abc.abstractmethod
    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        """
        run the experiment

        TODO: this function should really be split into to, and only have the body
        self.run(jobmanager, num_shots_per_circuit, qchip)
        return self.report()
        """
        pass

    @abc.abstractproperty
    def results(self):
        pass

    @abc.abstractmethod
    def plot_results(self, fig):
        pass

    @abc.abstractmethod
    def update_qchip(self, qchip):
        pass

    @abc.abstractmethod
    def update_gmm_manager(self, gmm_manager):
        pass
