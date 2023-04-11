import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic.toolchain as tc
import qubic.run as rc
from qubic.state_disc import GMMManager

ACC_BUFSIZE = 1000

class Rabi:
    """
    Define circuits, take data, and plot Rabi patterns
    """

    def __init__(self, drvqubit, readqubits, pulse_widths, qchip, fpga_config, channel_configs):
        """
        Create rabi circuits according to input parameters, then compile to asm binaries.
        """
        self.circuits = self._make_rabi_circuits(drvqubit,readqubits, pulse_widths, qchip)
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in readqubits}
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        compiled_progs = tc.run_compile_stage(self.circuits, fpga_config, qchip)
        self.raw_asm_progs = tc.run_assemble_stage(compiled_progs, channel_configs)

    def _make_rabi_circuits(self, drvqubit, readqubits, pulse_widths, qchip):
        """
        Make list of circuits used for rabi measurement. and the list of pulse width. So there will be a total of 
        1 circuits, each of which contains len(pulse_widths) measurements. A 400 us 
        delay is inserted between each measurement.
        """
        circuits = []
        self.pulse_widths = pulse_widths
        for twidth in pulse_widths:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6})
            cur_circ.append({'name': 'rabi', 'qubit': drvqubit, 'modi': {(0, 'twidth'): twidth }})
            for qubit in readqubits:
                cur_circ.append({'name': 'read', 'qubit': [qubit]})
            circuits.append(cur_circ)
        return circuits

    def run(self, circuit_runner, nsamples):
        """
        Run the rabi circuit with nsamples shots x pulse_width point.

        Parameters
        ----------
            circuit_runner : CircuitRunner object
            nsamples : int
        """
        reads_per_shot = 1 
        nshot = min(ACC_BUFSIZE//reads_per_shot,nsamples)
        navg = int(np.ceil(nsamples/nshot))
        self.s11 = {chan: np.zeros((len(self.pulse_widths), nshot*navg), dtype=np.complex128)
                    for chan in self.readout_chanmap.values()}
        print((len(self.pulse_widths), nshot*navg))        
        print('nshot',nshot,'navg',navg)        

        for i, raw_asm in enumerate(self.raw_asm_progs):
            print(i)
            circuit_runner.load_circuit(raw_asm)
            shots = circuit_runner.run_circuit(nshot, navg, reads_per_shot, delay=0.5)
            for chan, iqdata in shots.items():
                self.s11[chan][i] = iqdata #np.reshape(iqdata, (nshot*navg, len(self.pulse_widths))).T

        self._fit_gmm()

    def _fit_gmm(self):
        self.gmm_manager.fit(self.s11)
        self.gmm_manager.set_labels_maxtomin({chan: shots[np.argmin(self.pulse_widths)].flatten() 
                                              for chan, shots in self.s11.items()}, [0, 1])
        self.state_disc_shots = self.gmm_manager.predict(self.s11)
        self.ones_frac = {qubit: np.sum(self.state_disc_shots[qubit],axis=1) for qubit in self.state_disc_shots.keys()}
        self.zeros_frac = {qubit: np.sum(self.state_disc_shots[qubit] == 0,axis=1) for qubit in self.state_disc_shots.keys()}
