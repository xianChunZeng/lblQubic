import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic.toolchain as tc
from chipcalibration.abstract_calibration import AbstractCalibrationExperiment
import qubitconfig.qchip as qc
import qubic.job_manager as jm
import pdb
from scipy.signal import detrend

RINGDOWN_TIME = 40.e-6

class Vna(AbstractCalibrationExperiment):
    def __init__(self, drive_amp: float, freqs: np.ndarray, nsamples: int = 100, integration_time: float = 2.e-6, 
                 ringdown_time: float = RINGDOWN_TIME, qubit: str|list = None):
        if qubit is None:
            qubit = 'Q0'
        self.drive_amp = drive_amp
        self.freqs = freqs
        self.nsamples = nsamples
        self.integration_time = integration_time
        super().__init__(qubit, qubit)
        self.qubit = qubit
        self._make_circuits()

    def _make_circuits(self):
        self._circuits = []
        for freq in self.freqs:
            circuit = [{'name': 'delay', 't': RINGDOWN_TIME},
                       #{'name': 'barrier', 'scope': [f'{self.qubit}.rdrv', f'{self.qubit}.rdlo']},
                       {'name': 'pulse', 'amp': self.drive_amp, 'freq': freq, 'phase': 0, 'dest': f'{self.qubit}.rdrv',
                        'env': {
                                "env_func": "cos_edge_square",
                                "paradict": {"ramp_fraction": 0.25}}, 'twidth': self.integration_time},
                       {'name': 'delay', 't': 600.e-9, 'scope': [f'{self.qubit}.rdlo']},
                       {'name': 'pulse', 'amp': 1, 'freq': freq, 'phase': 0, 'dest': f'{self.qubit}.rdlo',
                        'env': {
                            "env_func": "square",
                            "paradict": { "phase": 0.0, "amplitude": 1.0}}, 'twidth': self.integration_time}]
            self._circuits.append(circuit)


    def run_and_report(self, jobmanager: jm.JobManager, num_shots_per_circuit: int = None, qchip: qc.QChip = None):
        if num_shots_per_circuit is not None:
            nshots = num_shots_per_circuit
        else:
            nshots = self.nsamples

        s11 = jobmanager.build_and_run_circuits(self._circuits, nshots, outputs=['s11'], reload_env=False)['s11']
        s11 = list(s11.values())[0]
        s11_avg = np.squeeze(np.average(s11, axis=1))

        phase = detrend(np.unwrap(np.angle(s11_avg)))
        amp = np.abs(s11_avg)

        self._results = {'amp': amp, 'phase': phase, 's11': s11_avg}

    @property
    def results(self):
        return self._results

    def plot_results(self, fig):
        pass

    def update_qchip(self, qchip):
        pass

    def update_gmm_manager(self, gmm_manager):
        pass



