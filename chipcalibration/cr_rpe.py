import matplotlib.pyplot as plt
import pygsti
import numpy as np


from pyrpe.src.quapack.pyRPE.quantum import Q as RPE_Experiment
from pyrpe.src.quapack.pyRPE import RobustPhaseEstimation

import qubic.pygsti.qupig as _qpig

class RPE_XR_Experiment:
    """
    RPE XR experiment class.  
    This class is used to generate the circuits for the RPE XR experiment and to
    run them on the QPU.
    """

    def __init__(self, control_qubit, target_qubit, max_max_depth=7):
        self.register = [control_qubit, target_qubit]
        if type(max_max_depth) is int:
            self.max_depths = [2 ** i for i in range(max_max_depth)]
        else:
            self.max_depths = max_max_depth
        self.state_pair_lookup = {(0,2):{('cos','+'):'01', ('cos','-'):'00',
                            ('sin','+'):'00', ('sin','-'):'01'},
                     (1,2):{('cos','+'):'00', ('cos','-'):'10',
                            ('sin','+'):'00', ('sin','-'):'10'},
                     (1,3):{('cos','+'):'11', ('cos','-'):'10',
                            ('sin','+'):'00', ('sin','-'):'01'}}
        self._make_pygsti_circuits()
        self.circuits = self.all_circuits_needing_data()

    def run_and_report(self,  pygsti_jobmanager, num_shots_per_circuit, qchip):
        ds = self._collect_data(pygsti_jobmanager, num_shots_per_circuit, qchip)
        self.ds = ds
    
    def rpe_circuits(self):
        return {circ: _qpig.qubic_instructions_from_pygsti_circuit(circ, self.register)
                for circ in self.all_circuits_needing_data}

    def all_circuits_needing_data(self):
        all_circs = []
        for trig_dict in [self.sin_dict, self.cos_dict]:
            for state_pair in self.state_pair_lookup.keys():
                all_circs += list(trig_dict[state_pair].values())
        return pygsti.remove_duplicates(all_circs)

    def make_cos_circuit(self, k, state_pair):
        if state_pair in [(0,2), (2,0)]:  #<01| for (1+cos(1/2 (theta_ix + theta_zx))
                                        #<00| for (1-cos(1/2 (theta_ix + theta_zx))
            circ = (pygsti.circuits.Circuit([[('Gxpi2',self.register[1])],
                                            ('Gxpi2',self.register[1])],
                                            line_labels=self.register)+
                    pygsti.circuits.Circuit([[('Gcr',self.register[0],self.register[1])]])*k)
        elif state_pair in [(1,2), (2,1)]:#<00| for (1+cos(1/2 (theta_zi - theta_zx))
                                        #<10| for (1+cos(1/2 (theta_zi - theta_zx))
            circ = (pygsti.circuits.Circuit([[('Gypi2',self.register[0])]],line_labels=self.register)+
                    pygsti.circuits.Circuit([[('Gypi2',self.register[1])]],line_labels=self.register)*3+
                    pygsti.circuits.Circuit([[('Gcr',self.register[0],self.register[1])]])*k+
                    pygsti.circuits.Circuit([[('Gypi2',self.register[0])]],line_labels=self.register)*3+
                    pygsti.circuits.Circuit([[('Gypi2',self.register[1])]],line_labels=self.register))
        elif state_pair in [(1,3), (3,1)]:#<11| for (1+cos(1/2 (theta_ix - theta_zx))
                                        #<10| for (1-cos(1/2 (theta_ix - theta_zx))
            circ = (pygsti.circuits.Circuit([[('Gxpi2',self.register[0])]],line_labels=self.register)*2+
                    pygsti.circuits.Circuit([[('Gxpi2',self.register[1])]],line_labels=self.register)*2+
                    pygsti.circuits.Circuit([[('Gcr',self.register[0],self.register[1])]])*k)
        else:
            assert False, "state_pair must be in [(0,2), (2,0), (1,2), (2,1), (1,3), (3,1)]"
        return circ

    def make_sin_circuit(self, k, state_pair):
        if state_pair in [(0,2), (2,0)]:  #<00| for (1+sin(1/2 (theta_ix + theta_zx))
                                        #<01| for (1+sin(1/2 (theta_ix - theta_zx))
            circ = (pygsti.circuits.Circuit([[('Gxpi2',self.register[1])],
                                            ('Gxpi2',self.register[1])],
                                            line_labels=self.register)+
                    pygsti.circuits.Circuit([[('Gcr',self.register[0],self.register[1])]])*k+
                    pygsti.circuits.Circuit([[('Gxpi2',self.register[1])]]))
        elif state_pair in [(1,2),(2,1)]:#<00| for (1+sin(1/2 (theta_zi - theta_zx))
                                        #<10| for (1+sin(1/2 (theta_zi - theta_zx))
            circ = (pygsti.circuits.Circuit([[('Gypi2',self.register[0])]],line_labels=self.register)+
                    pygsti.circuits.Circuit([[('Gypi2',self.register[1])]],line_labels=self.register)*3+
                    pygsti.circuits.Circuit([[('Gcr',self.register[0],self.register[1])]])*k+
                    pygsti.circuits.Circuit([[('Gxpi2',self.register[0])],
                                            [('Gypi2',self.register[1])]]))
        elif state_pair in [(1,3), (3,1)]:#<0-| for (1+sin(1/2 (theta_ix - theta_zx))
                                        #<01| for (1-sin(1/2 (theta_ix - theta_zx))
            circ = (pygsti.circuits.Circuit([[('Gxpi2',self.register[0])]],line_labels=self.register)*2+
                    pygsti.circuits.Circuit([[('Gxpi2',self.register[1])]],line_labels=self.register)*2+
                    pygsti.circuits.Circuit([[('Gcr',self.register[0],self.register[1])]])*k+
                    pygsti.circuits.Circuit([[('Gxpi2',self.register[0])]])*2+
                    pygsti.circuits.Circuit([[('Gxpi2',self.register[1])]]))
        else:
            assert False, "state_pair must be in [(0,2), (2,0), (1,2), (2,1), (1,3), (3,1)]"
        return circ

    def _make_pygsti_circuits(self):
        state_pairs = [(0,2), (1,2), (1,3)]
        self.sin_dict = {state_pair: {i: self.make_sin_circuit(i,state_pair, self.register) for i in self.max_depths} for state_pair in state_pairs}
        self.cos_dict = {state_pair: {i: self.make_cos_circuit(i,state_pair, self.register) for i in self.max_depths} for state_pair in state_pairs}



    def _collect_data(self, pygsti_jobmanager, num_shots_per_circuit, qchip):
        return pygsti_jobmanager.collect_dataset(self.circuits, num_shots_per_circuit, qchip=qchip)

 
 class CR_RPE_Analyzer:
    """
    Class for analyzing the results of a CR-RPE experiment
    """
    def __init__(self, max_depths, sin_circs, cos_circs):
        pass

    def process_rpe(self, dataset):
        # Post-process the RPE data from the pyGSTi dataset
        the_experiment = RPE_Experiment()
        for i in self.max_depths:
            the_experiment.process_sin(i, (int(dataset[self.sin_circs[i]]['00']), int(dataset[self.sin_circs[i]]['01'])))
            the_experiment.process_cos(i, (int(dataset[self.cos_circs[i]]['00']), int(dataset[self.cos_circs[i]]['01'])))
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

