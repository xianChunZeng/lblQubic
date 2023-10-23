import matplotlib.pyplot as plt
import numpy as np
import qubic.toolchain as tc
from qubic.state_disc import GMMManager
import logging

ACC_BUFSIZE = 1000

class Chevron:
    """
    Define circuits, take data, and plot Chevron patterns
    """

    def __init__(self, qubits, freqspan, nfreq, pulse_widths, qchip, fpga_config, 
                 channel_configs, centerfreq=None, ef=False, ef_rotate_back=False
                 ,rabigate='rabi'
                 ):
        """
        Create chevron circuits according to input parameters, then compile to asm binaries.
        """
        self.rabigate=rabigate
        self.circuits = self._make_chevron_circuits(qubits, freqspan, nfreq, pulse_widths, 
                                                    qchip, centerfreq, ef, ef_rotate_back)
        self.qubits = qubits
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in qubits}
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        compiled_progs = tc.run_compile_stage(self.circuits, fpga_config, qchip)
        self.raw_asm_progs = tc.run_assemble_stage(compiled_progs, channel_configs)

    def _make_chevron_circuits(self, qubits, freqspan, nfreq, pulse_widths, qchip, initial_centerfreq, ef, ef_rotate_back):
        """
        Make list of circuits used for chevron pattern measurement. Each circuit covers the 
        full range of frequencies, and a single pulse width. So there will be a total of 
        len(pulse_widths) circuits, each of which contains nfreq measurements. A 400 us 
        delay is inserted between each measurement.
        """
        circuits = []
        self.centerfreq={}
        self.freqoffsets = np.linspace(-freqspan/2, freqspan/2, nfreq)
        self.pulse_widths = pulse_widths
        for twidth in pulse_widths:
            cur_circ = []
            for freqoffset in self.freqoffsets:
                cur_circ.append({'name': 'delay', 't': 500.e-6})
                cur_circ.append({'name': 'barrier'})
                for qubit in qubits:
                    if initial_centerfreq is None:
                        centerfreq = qchip.qubits[qubit].freq_ef if ef else qchip.qubits[qubit].freq
                    else:
                        centerfreq = initial_centerfreq
                        if len(qubits) > 0:
                            logging.getLogger(__name__).warning(f'All qubits will be set to f = {centerfreq}')
                    self.centerfreq[qubit]=centerfreq
                    if ef:
                        cur_circ.append({'name': 'X90', 'qubit': [qubit]})
                        cur_circ.append({'name': 'X90', 'qubit': [qubit]})
                    freq = centerfreq + freqoffset
                    cur_circ.append({'name': self.rabigate, 'qubit': [qubit], 
                                     'modi': {(0, 'twidth'): twidth, (0, 'freq'): freq}})
                    if ef_rotate_back:
                        cur_circ.append({'name': 'X90', 'qubit': [qubit]})
                        cur_circ.append({'name': 'X90', 'qubit': [qubit]})
                    cur_circ.append({'name': 'read', 'qubit': [qubit]})
            circuits.append(cur_circ)
        return circuits

    def run(self, circuit_runner, nsamples):
        """
        Run the chevron circuit with nsamples shots per freq x pulse_width point.

        Parameters
        ----------
            circuit_runner : CircuitRunner object
            nsamples : int
        """
        self.s11 = circuit_runner.run_circuit_batch(self.raw_asm_progs, nsamples, len(self.freqoffsets), 
                                                    delay_per_shot=0)
        for chan in self.s11.keys():
            self.s11[chan] = np.transpose(self.s11[chan], (0, 2, 1)) # shape of s11 is pulse_widths x n_freq x nsamples
        self._fit_gmm()

    def _fit_gmm(self):
        self.gmm_manager.fit(self.s11)
        self.gmm_manager.set_labels_maxtomin({chan: shots[np.argmin(self.pulse_widths)].flatten() 
                                              for chan, shots in self.s11.items()}, [0, 1])
        self.state_disc_shots = self.gmm_manager.predict(self.s11)
        self.ones_frac = {qubit: np.sum(self.state_disc_shots[qubit], axis=2) for qubit in self.state_disc_shots.keys()}
        self.zeros_frac = {qubit: np.sum(self.state_disc_shots[qubit] == 0, axis=2) for qubit in self.state_disc_shots.keys()}

    def plot(self, qubit):
        plt.pcolormesh(self.pulse_widths, self.freqoffsets, self.ones_frac[qubit].T)
        plt.show()
