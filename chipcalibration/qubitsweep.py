import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic.qcvv.vnaqubit as vq

N_SAMPLES = 100
MAX_READ_QUBITS = 3

def run_qubit_sweep(drive_qubits, read_qubits, qchip, inst_cfg, drive_sweep_freqs, n_samples=N_SAMPLES):
    separation = np.zeros((len(drive_qubits), len(read_qubits), len(drive_sweep_freqs)))
    for i, qubit in enumerate(drive_qubits):
        vqubit = vq.c_vnaqubit(drive_qubit, qchip, inst_cfg)
        vqubit.seqs(drive_sweep_freqs, readqubitid=read_qubits, ef=False)
        vqubit.run(n_samples)
        separation[i] = np.asarray(vqubit.separation)

    return separation

