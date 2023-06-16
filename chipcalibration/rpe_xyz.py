import matplotlib.pyplot as plt
import pygsti
import numpy as np
from pygsti.processors import QubitProcessorSpec


from pygsti.modelpacks import smq1Q_XYZI

from pyrpe.src.quapack.pyRPE.quantum import Q as RPE_Experiment
from pyrpe.src.quapack.pyRPE import RobustPhaseEstimation

import qubic.pygsti.qupig as _qpig

def make_type1_cos_sequences(x_gate_label, rep_interval):
    """
    produces gate sequences for measuring X- or Y-type rpe phases
    """
    sequences = []
    for r in rep_interval:
        sequences.append([x_gate_label]*r)
    return sequences

def make_type1_sin_sequences(x_gate_label, rep_interval):
    """
    produces gate sequences for measuring X- or Y-type rpe phases
    """
    sequences = []
    for r in rep_interval:
        sequences.append([x_gate_label]*(r+1))
    return sequences

def make_z_cos_sequences(x_gate_label, y_gate_label, rep_interval):
    """
    produces gate sequences for measuring compiled Z rpe phases
    """
    sequences = []
    for r in rep_interval:
        s = [y_gate_label]
        s = s + [y_gate_label, x_gate_label, y_gate_label, y_gate_label, y_gate_label]*r
        s = s + [y_gate_label]*3
        sequences.append(s)
    return sequences

def make_z_sin_sequences(x_gate_label, y_gate_label, rep_interval):
    """
    produces gate sequences for measuring compiled Z rpe phases
    """
    sequences = []
    for r in rep_interval:
        s = [x_gate_label]
        s = s + [y_gate_label, x_gate_label, y_gate_label, y_gate_label, y_gate_label]*r
        s = s + [y_gate_label]*3
        sequences.append(s)
    return sequences

def convert_sequence_to_pygsti(sequence):
    circ = pygsti.circuits.Circuit([sequence[0]])
    for s in sequence[1:]:
        circ += pygsti.circuits.Circuit([s])
    return circ


class RPE_XYZ_Experiment:
    """
    RPE XYZ experiment

    collects data for postprocessing and optimization with RPEAnalyzer and RPEOptimizer
    """

    def __init__(self, target_qid, calibration, max_max_depth=7):
        self.qid = target_qid
        self.max_depths = [2 ** i for i in range(max_max_depth)]
        pspec = QubitProcessorSpec(num_qubits=1, qubit_labels=[target_qid],
                                   gate_names=['Gxpi2', 'Gzpi2'], geometry='line')
        self.processor_spec = pspec
        self.target_model = pygsti.models.modelconstruction.create_explicit_model(pspec)
        self.calibration = calibration
        self.x_instructions = self.make_x_instructions(calibration)
        self.y_instructions = self.make_y_instructions(calibration)
        self.z_instructions = self.make_z_instructions(calibration)

    # def update_drive_amp(self, drive_amp):
    #     self.qchip.update(('Gates', f'{self.qid}X90', 0, 'amp'), drive_amp)
    #
    # def update_frequency(self, frequency):
    #     self.qchip.update(('Qubits', self.qid, 'freq'), frequency)

    def run(self, pygsti_jobmanager, num_shots_per_circuit, qchip):
        datasets = dict()
        datasets['X'] = self.collect_x_dataset(pygsti_jobmanager, num_shots_per_circuit, qchip)
        datasets['Y'] = self.collect_y_dataset(pygsti_jobmanager, num_shots_per_circuit, qchip)
        datasets['Z'] = self.collect_z_dataset(pygsti_jobmanager, num_shots_per_circuit, qchip)
        return datasets

    def collect_x_dataset(self, pygsti_jobmanager, num_shots_per_circuit, qchip):
        return pygsti_jobmanager.collect_dataset(self.x_instructions, num_shots_per_circuit, qchip)

    def collect_y_dataset(self, pygsti_jobmanager, num_shots_per_circuit, qchip):
        return pygsti_jobmanager.collect_dataset(self.y_instructions, num_shots_per_circuit, qchip)

    def collect_z_dataset(self, pygsti_jobmanager, num_shots_per_circuit, qchip):
        return pygsti_jobmanager.collect_dataset(self.z_instructions, num_shots_per_circuit, qchip)



    def make_x_circuits(self):
        x_cos_dict = {i: convert_sequence_to_pygsti(s)
                      for i, s in enumerate(make_type1_cos_sequences(('Gxpi2', self.qid), self.max_depths))}
        x_sin_dict = {i: convert_sequence_to_pygsti(s)
                      for i, s in enumerate(make_type1_sin_sequences(('Gxpi2', self.qid), self.max_depths))}
        return x_cos_dict, x_sin_dict

    def make_y_circuits(self):
        y_cos_dict = {i: convert_sequence_to_pygsti(s)
                      for i, s in enumerate(make_type1_cos_sequences(('Gypi2', self.qid), self.max_depths))}
        y_sin_dict = {i: convert_sequence_to_pygsti(s)
                      for i, s in enumerate(make_type1_sin_sequences(('Gypi2', self.qid), self.max_depths))}
        return y_cos_dict, y_sin_dict

    def make_z_circuits(self):
        z_cos_dict = {i: convert_sequence_to_pygsti(s)
                      for i, s in enumerate(make_z_cos_sequences(('Gxpi2', self.qid), ('Gypi2', self.qid), self.max_depths))}
        z_sin_dict = {i: convert_sequence_to_pygsti(s)
                      for i, s in enumerate(make_z_sin_sequences(('Gxpi2', self.qid), ('Gypi2', self.qid), self.max_depths))}
        return z_cos_dict, z_sin_dict


    @property
    def all_circuits_needing_data(self):
        all_circs = list(self.cos_circs.values()) + [c for c in self.sin_circs.values() if
                                                     c not in self.cos_circs.values()]
        return all_circs

    def make_cos_circuit(self, power: int):
        return (pygsti.circuits.Circuit([[('Gxpi2', self.register[0])]]) * 2 +
                pygsti.circuits.Circuit([[('Gcr', self.register[0], self.register[1])]]) * power +
                pygsti.circuits.Circuit([[('Gxpi2', self.register[0])]]) * 2)

    def make_sin_circuit(self, power: int):
        return (pygsti.circuits.Circuit([[('Gxpi2', self.register[0])]]) * 2 +
                pygsti.circuits.Circuit([[('Gcr', self.register[0], self.register[1])]]) * power +
                pygsti.circuits.Circuit([[('Gxpi2', self.register[0])]]) * 2 +
                pygsti.circuits.Circuit([[('Gzpi2', self.register[1])]]) * 2 +
                pygsti.circuits.Circuit([[('Gxpi2', self.register[1])]]) * 3)

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

class RPEAnalyzer:
    def __init__(self, rpe_xyz_experiment):
        self.datasets = rpe_xyz_experiment.datasets


    def process_rpe(self, dataset):
        # Post-process the RPE data from the pyGSTi dataset
        the_experiment = RPE_Experiment()
        for i in self.max_depths:
            the_experiment.process_sin(i,
                                       (int(dataset[self.sin_circs[i]]['00']), int(dataset[self.sin_circs[i]]['01'])))
            the_experiment.process_cos(i,
                                       (int(dataset[self.cos_circs[i]]['00']), int(dataset[self.cos_circs[i]]['01'])))
        return RobustPhaseEstimation(the_experiment)

    def plot_rpe_verbose(self, dataset, num_shots, rpe_results):
        # 1st process the dataset -- probably should design so that I don't have to process every time
        fit = rpe_results
        target_ds = self.simulate_target(num_shots)
        self.target_ds = target_ds
        target_fit = self.process_rpe(target_ds)
        fig, axs = plt.subplots(3, 1)
        axs[0].semilogx(self.max_depths, [dataset[c]['00'] for c in self.sin_circs.values()])
        axs[0].semilogx(self.max_depths, [target_ds[c]['00'] for c in self.sin_circs.values()])
        axs[0].set_title("00 counts on sin circuits")
        axs[1].semilogx(self.max_depths, [dataset[c]['01'] for c in self.cos_circs.values()])
        axs[1].plot(self.max_depths, [target_ds[c]['01'] for c in self.cos_circs.values()])
        axs[1].set_title("01 counts on cos circuits")
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
        return pygsti_jobmanager.collect_dataset(self.circuits, num_shots_per_circuit, qchip=qchip)
