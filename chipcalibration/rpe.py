import matplotlib.pyplot as plt
import pygsti
import numpy as np

from pygsti.modelpacks import smq1Q_XZ


from pyrpe.src.quapack.pyRPE.quantum import Q as RPE_Experiment
from pyrpe.src.quapack.pyRPE import RobustPhaseEstimation

import qubic.pygsti.qupig as _qpig


class RpeX90:
    """
    Class that can:
        1) make RPE circuit experiment designs
        2) fit curves to RPE results
        3) plot curves

    """
    sin_circs: dict

    def __init__(self, processor_spec, target_qubit, max_max_depth=7):
        self.qid = target_qubit
        self.max_depths = [2 ** i for i in range(max_max_depth)]
        self.processor_spec = processor_spec
        self.target_model = pygsti.models.modelconstruction.create_explicit_model(processor_spec)
        self._make_pygsti_circuits()
        self.circuits = self.rpe_circuits()

    # def update_drive_amp(self, drive_amp):
    #     self.qchip.update(('Gates', f'{self.qid}X90', 0, 'amp'), drive_amp)
    #
    # def update_frequency(self, frequency):
    #     self.qchip.update(('Qubits', self.qid, 'freq'), frequency)

    def run_and_report(self,  pygsti_jobmanager, num_shots_per_circuit, qchip):
        ds = self._collect_data(pygsti_jobmanager, num_shots_per_circuit, qchip)
        self.rpe_results = self.process_rpe(ds)
        self.plot_rpe_verbose(ds, num_shots_per_circuit, self.rpe_results)
        last_good_estimate_generation = self.rpe_results.check_unif_local(historical=True)
        print(f'Last good generation: {last_good_estimate_generation}')
        print(f'Estimated phase: {self.rpe_results.angle_estimates[last_good_estimate_generation]}')
        return self.rpe_results

    def rpe_circuits(self):
        return {circ: _qpig.qubic_instructions_from_pygsti_circuit(circ, [self.qid])
                for circ in self.all_circuits_needing_data}

    @property
    def all_circuits_needing_data(self):
        all_circs = list(self.cos_circs.values()) + [c for c in self.sin_circs.values() if
                                                     c not in self.cos_circs.values()]
        return all_circs

    def make_cos_circuit(self, power: int):
        return pygsti.circuits.Circuit([[('Gxpi2', self.qid)]])*power

    def make_sin_circuit(self, power: int):
        return pygsti.circuits.Circuit([[('Gxpi2', self.qid)]])*power + pygsti.circuits.Circuit(
            [[('Gxpi2', self.qid)]])

    def _make_pygsti_circuits(self):
        self.sin_circs = {i: self.make_sin_circuit(i) for i in self.max_depths}
        self.cos_circs = {i: self.make_cos_circuit(i) for i in self.max_depths}

    def simulate_target(self, num_samples):
        """
        Simulate the pygsti circuits at the target model
        :return:
        """
        target_ds = pygsti.data.simulate_data(self.target_model, self.all_circuits_needing_data,
                                              num_samples=num_samples, sample_error='multinomial', seed=None)
        return target_ds

    def process_rpe(self, dataset):
        # Post-process the RPE data from the pyGSTi dataset
        the_experiment = RPE_Experiment()
        for i in self.max_depths:
            the_experiment.process_sin(i, (int(dataset[self.sin_circs[i]]['0']), int(dataset[self.sin_circs[i]]['1'])))
            the_experiment.process_cos(i, (int(dataset[self.cos_circs[i]]['0']), int(dataset[self.cos_circs[i]]['1'])))
        return RobustPhaseEstimation(the_experiment)

    def plot_rpe_verbose(self, dataset, num_shots, rpe_results):
        # 1st process the dataset -- probably should design so that I don't have to process every time
        fit = rpe_results
        target_ds = self.simulate_target(num_shots)
        target_fit = self.process_rpe(target_ds)
        fig, axs = plt.subplots(3, 1)
        axs[0].semilogx(self.max_depths, [dataset[c]['1'] for c in self.sin_circs.values()])
        axs[0].semilogx(self.max_depths, [target_ds[c]['1'] for c in self.sin_circs.values()])
        axs[0].set_title("1 counts on sin circuits")
        axs[1].semilogx(self.max_depths, [dataset[c]['1'] for c in self.cos_circs.values()])
        axs[1].plot(self.max_depths, [target_ds[c]['1'] for c in self.cos_circs.values()])
        axs[1].set_title("1 counts on cos circuits")
        axs[2].semilogx(self.max_depths, fit.angle_estimates)
        axs[2].plot(self.max_depths, target_fit.angle_estimates, c='orange')
        axs[2].plot((self.max_depths[0], self.max_depths[-1]), (np.pi / 2, np.pi / 2), c='red')
        num_estimates = len(self.max_depths)
        axs[2].fill_between(self.max_depths,
                            [fit.angle_estimates[i] + 1 / (i + 1) for i in range(num_estimates)],
                            [fit.angle_estimates[i] - 1 / (i + 1) for i in range(num_estimates)], color='blue')
        axs[2].set_xlabel("Circuit depth")
        plt.tight_layout()
        plt.show()

    def _collect_data(self, pygsti_jobmanager, num_shots_per_circuit, qchip):
        """

        """
        return pygsti_jobmanager.collect_dataset(self.circuits, num_shots_per_circuit, qchip)
