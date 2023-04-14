import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic.toolchain as tc
import qubic.run as rc
from qubic.state_disc import GMMManager
from scipy.optimize import curve_fit

ACC_BUFSIZE = 1000

class Ramsey:
    """
    Define circuits, take data, and plot Ramsey patterns
    """
    def __init__(self, qubits, delaytime, qchip, fpga_config, channel_configs):
        """
        Create rabi circuits according to input parameters, then compile to asm binaries.
        """
        self.circuits = self._make_ramsey_circuit(qubits, delaytime, qchip)
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in qubits}
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        compiled_progs = tc.run_compile_stage(self.circuits, fpga_config, qchip)
        self.raw_asm_progs = tc.run_assemble_stage(compiled_progs, channel_configs)
        self.qubits = qubits
        
    def _make_ramsey_circuit(self, qubits, delaytime, qchip):
        """
        Make a circuit used for ramsey measurement with the list of delaytime. So there will be a total of 
        1 circuits, each of which contains len(delaytime) measurements. A 400 us 
        delay is inserted between each measurement.
        Limitations: 
        1. Length of the delaytime list is restricted to 1024 because of memory
        2. Due to command buffer limitation we can have 341 of x90->delay->X90
        
        """
        circuit = []
        self.delaytime = delaytime
#        cur_circ = []
        for dtime in delaytime:
            for qubit in qubits:
                circuit.append({'name': 'X90', 'qubit': [qubit]})
                circuit.append({'name': 'delay', 't': dtime, 'qubit':[qubit]})
                circuit.append({'name': 'X90', 'qubit': [qubit]})
                circuit.append({'name': 'read', 'qubit': [qubit]})
                circuit.append({'name': 'delay', 't': 400.e-6, 'qubit': [qubit]})
        return circuit
        
        
    def run(self, circuit_runner, nsamples):
        """
        Run the ramsey circuit with nsamples shots.
    
        Parameters
        ----------
            circuit_runner : CircuitRunner object
            nsamples : int
        """
        circuit_runner.load_circuit(self.raw_asm_progs)
        #Total delay is sum of 500 us for every measurment and all the delays
        
        self.s11 = circuit_runner.run_circuit(nsamples,len(self.delaytime), delay_per_shot=(500.e-6*len(self.delaytime))+sum(self.delaytime))

    def _fit_gmm(self):
        self.gmm_manager.fit(self.s11)
        self.gmm_manager.set_labels_maxtomin({chan: data[:,0] for chan, data in self.s11.items()},labels_maxtomin = [1,0])
        self.state_disc_shots = self.gmm_manager.predict(self.s11)
        self.ones_frac = {qubit: np.sum(self.state_disc_shots[qubit], axis=0) for qubit in self.state_disc_shots.keys()}
        self.zeros_frac = {qubit: np.sum(self.state_disc_shots[qubit] == 0, axis=0) for qubit in self.state_disc_shots.keys()}

    def _cos_exp(self, x, A, B, drive_freq, phi,exp_decay):
        return A*np.exp(-x/exp_decay)*np.cos(2*np.pi*x*drive_freq - phi) + B
        

    def fit_ramsey_freq(self,prior_fit_params):
        #cosine decaying exponentially with offset
        self.param = {}
        for qubit in self.qubits:
            self.param[qubit] = curve_fit(self._cos_exp, self.delaytime, self.ones_frac[qubit], prior_fit_params[qubit])
        
	
