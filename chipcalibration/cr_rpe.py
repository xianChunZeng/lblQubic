import matplotlib.pyplot as plt
import pygsti
import numpy as np



from pyrpe.src.quapack.pyRPE.quantum import Q as RPE_Experiment
from pyrpe.src.quapack.pyRPE import RobustPhaseEstimation

import qubic.pygsti.qupig as _qpig

def rectify_angle(theta):
    if theta > np.pi:
        theta -= 2*np.pi
    return theta

class CR_RPE_Experiment:
    """
    RPE CR experiment class.  
    This class is used to generate the circuits for the RPE XR experiment and to
    run them on the QPU.
    """

    def __init__(self, control_qubit, target_qubit, calibration, max_max_depth=7):
        self.calibration = calibration
        self.target_qid = target_qubit
        self.control_qid = control_qubit
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
        
        self.state_pairs = list(self.state_pair_lookup.keys())
        self._make_pygsti_circuits()
        self.pygsti_circuits = self.all_pygsti_circuits_needing_data()
        self.qubic_instructions = {circ : self.transpile_qubic_instructions(circ, calibration) for circ in self.pygsti_circuits}

    def run(self,  pygsti_jobmanager, num_shots_per_circuit, qchip):
        self.ds = pygsti_jobmanager.collect_dataset(self.qubic_instructions, num_shots_per_circuit, qchip)
                                                    
    def rpe_circuits(self):
        return {circ: _qpig.qubic_instructions_from_pygsti_circuit(circ, self.register)
                for circ in self.all_circuits_needing_data}

    def all_pygsti_circuits_needing_data(self):
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
        self.sin_dict = {state_pair: {i: self.make_sin_circuit(i,state_pair) for i in self.max_depths} for state_pair in state_pairs}
        self.cos_dict = {state_pair: {i: self.make_cos_circuit(i,state_pair) for i in self.max_depths} for state_pair in state_pairs}
        
    def make_cr_pulse(self, calibration, phase_offset=0):
        return [
            {
            'name': 'pulse',
            "freq": f"{self.target_qid}.freq",
            "dest": f"{self.control_qid}.qdrv",
            "twidth": calibration['twidth']/2.0,
            "amp": calibration['c_amp'],
            "phase": calibration['c_phase'] + phase_offset,
            "env": [
                {
                    "env_func": "cos_edge_square",
                    "paradict": {
                        "ramp_fraction": 0.25,
                    }
                }
            ]
        }, {
            'name': 'pulse',
            "freq": f"{self.target_qid}.freq",
            "dest": f"{self.target_qid}.qdrv",
            "twidth": calibration['twidth']/2.0,
            "amp": calibration['t_amp'],
            "phase": calibration['t_phase'] + phase_offset, 
            "env": [
                {
                    "env_func": "cos_edge_square",
                    "paradict": {
                        "ramp_fraction": 0.25,
                    }
                }
            ]
            }, {'name': 'barrier'}
        ]
    
    def make_echoed_cr_pulse_schedule(self, calibration):
        circ = self.make_cr_pulse(calibration)
        circ.append({'name': 'delay', 't': 10.e-9})
        circ.append({'name': 'X90', 'qubit': [self.control_qid]})
        circ.append({'name': 'X90', 'qubit': [self.control_qid]})
        circ.append({'name': 'barrier'})
        circ.append({'name': 'delay', 't': 10.e-9})
        circ.extend(self.make_cr_pulse(calibration, np.pi))
        circ.append({'name': 'X90', 'qubit': [self.control_qid]})
        circ.append({'name': 'X90', 'qubit': [self.control_qid]})
        circ.append({'name': 'barrier'})
        return circ
        
    def transpile_qubic_instructions(self, pygsti_circuit, cr_calibration):
        """
        converts the pygsti circuit into qubic instructions
        
        """
        qubic_circuit = list()
        qubic_circuit.append({'name': 'delay', 't': 400.e-6})
        for layer in pygsti_circuit:
            if layer.name == 'Gcr':
                qubic_circuit.extend( self.make_echoed_cr_pulse_schedule(cr_calibration))
            else:
                qubic_circuit.extend(_qpig.parse_layer(layer))
            qubic_circuit.append({'name': 'barrier', 'qubit': list(self.register)})
        for qid in self.register:
            qubic_circuit.append({'name': 'read', 'qubit': [qid]}, )
        return qubic_circuit
            



    def _collect_data(self, pygsti_jobmanager, num_shots_per_circuit, qchip):
        return pygsti_jobmanager.collect_dataset(self.circuits, num_shots_per_circuit, qchip=qchip)

class CR_RPE_Analyzer:
    """
    Class for analyzing the results of a CR-RPE experiment
    """
    def __init__(self, cr_rpe_exp: CR_RPE_Experiment):
        self.experiment = cr_rpe_exp
        self.state_pairs = cr_rpe_exp.state_pairs
        self.state_pair_lookup = cr_rpe_exp.state_pair_lookup
        self.max_depths = cr_rpe_exp.max_depths
        self.sin_dict = cr_rpe_exp.sin_dict
        self.cos_dict = cr_rpe_exp.cos_dict


    def process_rpe(self, dataset):
        #Post-process the RPE data from the pyGSTi dataset
        the_experiments = {}
        analyses = {}
        for state_pair in self.state_pairs:
            the_experiments[state_pair] = RPE_Experiment()

        for state_pair in self.state_pairs:
            cos_plus = self.state_pair_lookup[state_pair]['cos','+']
            cos_minus = self.state_pair_lookup[state_pair]['cos','-']
            sin_plus = self.state_pair_lookup[state_pair]['sin','+']
            sin_minus = self.state_pair_lookup[state_pair]['sin','-']
            for i in self.max_depths:
                the_experiments[state_pair].process_sin(i,(int(dataset[self.sin_dict[state_pair][i]][sin_plus]),
                                                        int(dataset[self.sin_dict[state_pair][i]][sin_minus])))
                the_experiments[state_pair].process_cos(i,(int(dataset[self.cos_dict[state_pair][i]][cos_plus]),
                                                        int(dataset[self.cos_dict[state_pair][i]][cos_minus])))
        for state_pair in self.state_pairs:
            analyses[state_pair] = RobustPhaseEstimation(the_experiments[state_pair])

        last_good_estimate_generations = {}

        for state_pair in self.state_pairs:
            last_good_estimate_generations[state_pair] = analyses[(state_pair)].check_unif_local(historical=True)

        last_good_estimate_generation = min(list(last_good_estimate_generations.values()))
        print('Generation of last good estimate: ', last_good_estimate_generation)

        #Now shift angles into correct principle range.
        for state_pair in self.state_pairs:
            if state_pair in [(1,2),(1,3)]:#If you're getting funny results for your (0,2) pair estimate, try including (0,2) here as well.
                analyses[state_pair].angle_estimates = [rectify_angle(theta) for theta in analyses[state_pair].angle_estimates]#rectify_angle(analyses[state_pair].angle_estimates)

        #Turn lin. comb. estimates into direct phase estimates.
        theta_ix_estimates = 0.5 * (np.array(analyses[(0,2)].angle_estimates) + np.array(analyses[(1,3)].angle_estimates))
        theta_zx_estimates = 0.5 * (np.array(analyses[(0,2)].angle_estimates) - np.array(analyses[(1,3)].angle_estimates))
        theta_zi_estimates = np.array(analyses[(1,2)].angle_estimates)+theta_zx_estimates

        #Get final trusted estimates:
        #What are the actual phase estimates (in radians) for all depths?
        print(f'Trusted theta_ix : {theta_ix_estimates[last_good_estimate_generation]}')
        print(f'Trusted theta_zx : {theta_zx_estimates[last_good_estimate_generation]}')
        print(f'Trusted theta_zi : {theta_zi_estimates[last_good_estimate_generation]}')


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

