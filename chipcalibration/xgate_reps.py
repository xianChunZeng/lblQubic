import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


class XGateReps:
    """
    repeated X-pi/2 gates with different amps with 
        4n+2 reps, and measure
    """
    def __init__(self, target_qid, center_amp, sigma_amp, num_partitions=100, num_reps=5):
        self.qid = target_qid
        self.amplitude_partitions = np.linspace(center_amp-sigma_amp, center_amp+sigma_amp, num_partitions)
        self.circuits = self._make_circuits()
        
    def _make_circuits(self):
        circuits = []
        for amp in self.amplitude_partitions:
            for reps in range(self.num_reps):
                cur_circ = []
                cur_circ.append({'name': 'delay', 't': 400e-6})
                for _ in range(4*reps+2):
                    cur_circ.append({'name': 'X90', 'qubit': [self.qid], 'modi': {(0, 'amp'): amp}})
                cur_circ.append({'name': 'read', 'qubit': [self.qid]})
                circuits.append(cur_circ)
        return circuits
                    
        