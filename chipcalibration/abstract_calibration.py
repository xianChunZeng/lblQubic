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