import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic.toolchain as tc
import qubic.run as rc
from qubic.state_disc import GMMManager

ACC_BUFSIZE = 1000

class Chevron:

    def __init__(self, qubits, freqspan, nfreq, pulse_widths, qchip, fpga_config, channel_configs):
        self.circuits = self._make_chevron_circuits(qubits, freqspan, nfreq, pulse_widths, qchip)
        self.qubits = qubits
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in qubits}
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        compiled_progs = tc.run_compile_stage(self.circuits, fpga_config, qchip)
        self.raw_asm_progs = tc.run_assemble_stage(compiled_progs, channel_configs)

    def _make_chevron_circuits(self, qubits, freqspan, nfreq, pulse_widths, qchip):
        circuits = []
        self.freqoffsets = np.linspace(-freqspan/2, freqspan/2, nfreq)
        self.pulse_widths = pulse_widths
        for twidth in pulse_widths:
            cur_circ = []
            for freqoffset in self.freqoffsets:
                for qubit in qubits:
                    freq = qchip.qubits[qubit].freq + freqoffset
                    cur_circ.append({'name': 'rabi', 'qubit': [qubit], 
                                     'modi': {(0, 'twidth'): twidth, (0, 'fcarrier'): freq}})
                    cur_circ.append({'name': 'read', 'qubit': [qubit]})
                cur_circ.append({'name': 'delay', 't': 400.e-6})
            circuits.append(cur_circ)
        return circuits

    def run(self, circuit_runner, nsamples):
        reads_per_shot = len(self.freqoffsets)
        nshot = ACC_BUFSIZE//reads_per_shot
        navg = int(np.ceil(nsamples/nshot))
        self.s11 = {chan: np.zeros((len(self.pulse_widths), len(self.freqoffsets), nshot*navg), dtype=np.complex128) 
                    for chan in self.readout_chanmap.values()}

        for i, raw_asm in enumerate(self.raw_asm_progs):
            circuit_runner.load_circuit(raw_asm)
            shots = circuit_runner.run_circuit(nshot, navg, reads_per_shot, delay=0.5)
            for chan, iqdata in shots:
                self.s11[chan][i] = np.reshape(iqdata, (nshot*navg, len(self.freqoffsets)))

        self._fit_gmm()

    def _fit_gmm(self):
        self.gmm_manager.fit(self.s11)
        self.gmm_manager.set_labels_maxtomin({chan: shots[np.argmin(self.pulse_widths)].flatten() 
                                              for chan, shots in self.s11.items()}, [0, 1])
        self.state_disc_shots = self.gmm_manager.predict(self.s11)
        self.ones_frac = {qubit: np.sum(self.state_disc_shots[qubit], axis=2) for qubit in self.state_disc_shots.keys()}
        self.zeros_frac = {qubit: np.sum(self.state_disc_shots[qubit] == 0, axis=2) for qubit in self.state_disc_shots.keys()}
