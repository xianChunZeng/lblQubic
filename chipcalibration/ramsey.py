import matplotlib.pyplot as plt
import numpy as np
import qubic.toolchain as tc
from qubic.state_disc import GMMManager
from scipy.optimize import curve_fit
import copy

ACC_BUFSIZE = 1000

class Ramsey:
    """
    Define circuits, take data, and plot Ramsey patterns

    TODO: add FFT for coarse period finding
    """
    def __init__(self, qubits, delaytime, qchip, fpga_config, channel_configs, freq_offs_dict=None):
        """
        Create rabi circuits according to input parameters, then compile to asm binaries.
        """
        self.qubits = qubits
        self.circuits = self._make_ramsey_circuit(qubits, delaytime, qchip, freq_offs_dict)
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in qubits}
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        compiled_progs = tc.run_compile_stage(self.circuits, fpga_config, qchip)
        self.raw_asm_progs = tc.run_assemble_stage(compiled_progs, channel_configs)
        
    def _make_ramsey_circuit(self, qubits, delaytime, qchip, freq_offs_dict=None):
        """
        Make a circuit used for ramsey measurement with the list of delaytime. So there will be a total of 
        1 circuits, each of which contains len(delaytime) measurements. A 400 us 
        delay is inserted between each measurement.
        Limitations: 
        1. Length of the delaytime list is restricted to 1024 because of memory
        2. Due to command buffer limitation we can have 341 of x90->delay->X90

        TODO: have global class attribute circuit that gets set
        
        """
        circuit = []
        self.delaytime = delaytime
        if freq_offs_dict is None:
            freq_offs_dict = {q: 0 for q in qubits}
        for dtime in delaytime:
            for qubit in qubits:
                drivefreq = qchip.gates[qubit + 'X90'].contents[0].fcarrier + freq_offs_dict[qubit]
                circuit.append({'name': 'X90', 'qubit': [qubit], 'modi':{(0, 'fcarrier'): drivefreq}})
                circuit.append({'name': 'delay', 't': dtime, 'qubit':[qubit]})
                circuit.append({'name': 'X90', 'qubit': [qubit], 'modi':{(0, 'fcarrier'): drivefreq}})
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
        self._fit_gmm()

    def _fit_gmm(self):
        self.gmm_manager.fit(self.s11)
        self.gmm_manager.set_labels_maxtomin({chan: data[:,0] for chan, data in self.s11.items()},labels_maxtomin = [1,0])
        self.state_disc_shots = self.gmm_manager.predict(self.s11)
        self.ones_frac = {qubit: np.sum(self.state_disc_shots[qubit], axis=0) for qubit in self.state_disc_shots.keys()}
        self.zeros_frac = {qubit: np.sum(self.state_disc_shots[qubit] == 0, axis=0) for qubit in self.state_disc_shots.keys()}

    def plot_fits(self, qubit):
        plt.plot(self.delaytime, self.ones_frac[qubit], '.')
        try:
            plt.plot(self.delaytime, self._cos_exp(self.delaytime, *self.fit_params[qubit][0]))
        except KeyError:
            pass


    def _cos_exp(self, x, A, B, drive_freq, phi,exp_decay):
        return A*np.exp(-x/exp_decay)*np.cos(2*np.pi*x*drive_freq - phi) + B
        

    def fit_ramsey_freq(self, prior_fit_params, use_fft=True):
        """
        cosine decaying exponentially with offset
        params are [A, B, drive_freq, phi, exp_decay]
        """
        self.fit_params = {}
        prior_fit_params = copy.deepcopy(prior_fit_params)
        for qubit in self.qubits:
            if use_fft:
                freq_ind_max = np.argmax(np.abs(np.fft.rfft(self.ones_frac[qubit])[1:])) + 1
                freq_max = np.fft.rfftfreq(len(self.delaytime), np.diff(self.delaytime)[0])[freq_ind_max]
                prior_fit_params[qubit][2] = freq_max
            try:
                self.fit_params[qubit] = curve_fit(self._cos_exp, self.delaytime, self.ones_frac[qubit], prior_fit_params[qubit])
            except RuntimeError:
                print('{} could not be fit')



class RamseyOptimize:

    def __init__(self, qubits, delaytime, qchip, fpga_config, channel_configs):
        self.qubits = qubits
        self.delaytime = delaytime
        self.qchip = qchip
        self.fpga_config = fpga_config
        self.channel_configs = channel_configs
        self.initial_ramsey = Ramsey(qubits, delaytime, qchip, fpga_config, channel_configs)

    def run_optimize_step(self, circuit_runner, nsamples, prior_fit_params):
        #run initial sweep
        self.initial_ramsey.run(circuit_runner, nsamples)
        self.initial_ramsey.fit_ramsey_freq(prior_fit_params)
        initial_params = copy.deepcopy(self.initial_ramsey.fit_params)
        for qubit in self.qubits:
            self.initial_ramsey.plot_fits(qubit)
            plt.show()

        #run positive freq
        freq_offs_dict = {q: initial_params[q][0][2] for q in initial_params.keys()}
        pos_ramsey = Ramsey(self.qubits, self.delaytime, self.qchip, self.fpga_config, 
                            self.channel_configs, freq_offs_dict)
        pos_ramsey.run(circuit_runner, nsamples)
        pos_ramsey.fit_ramsey_freq({q: initial_params[q][0] for q in self.qubits})
        for qubit in self.qubits:
            pos_ramsey.plot_fits(qubit)
            plt.title(qubit + ' + {} Hz'.format(freq_offs_dict[qubit]))
            plt.show()

        #run negative freq
        freq_offs_dict = {q: -initial_params[q][0][2] for q in initial_params.keys()}
        neg_ramsey = Ramsey(self.qubits, self.delaytime, self.qchip, self.fpga_config, 
                            self.channel_configs, freq_offs_dict)
        neg_ramsey.run(circuit_runner, nsamples)
        neg_ramsey.fit_ramsey_freq({q: initial_params[q][0] for q in self.qubits})
        for qubit in self.qubits:
            neg_ramsey.plot_fits(qubit)
            plt.title(qubit + ' - {} Hz'.format(freq_offs_dict[qubit]))
            plt.show()

        return pos_ramsey.fit_params, neg_ramsey.fit_params

    def update_qubit_freq(self, qubit, sign, qchip):
        """
        update the qubit frequency in the qchip object

        sign : int
            if sign is +1, add the initial fit freq, if -1 subtract
        """

        if np.abs(sign) != 1:
            raise Exception('sign must be +/-1')

        qchip.qubits[qubit].freq += sign * self.initial_ramsey.fit_params[qubit][0][2]
