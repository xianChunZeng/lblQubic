from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import sys


class dragalpha:
    """
    Define circuits, take data, and plot dragalpha
    """

    def __init__(self, qubit, alphas):
        """
        Create dragalpha circuits according to input parameters, then compile to asm binaries.
        """
        self.qubit = qubit
        self.alphas = alphas
        self.circuits = self._make_dragalpha_circuits()

    def _make_dragalpha_circuits(self,delaybeforecircuit=600e-6):
        """
        Make list of circuits used for drag alpha measurement. 
        """

        circuits=[]
        for alpha in self.alphas:
            x90alpha=[{'name': 'X90', 'qubit': self.qubit, 'modi' : {(0,"env",0,"paradict","alpha"):alpha}}]
            circuit=[{'name': 'delay', 't': delaybeforecircuit}]
            circuit.extend(x90alpha)
            circuit.append({'name': 'virtual_z', 'qubit': self.qubit, 'phase': -np.pi/2})
            circuit.extend(x90alpha)
            circuit.extend(x90alpha)
            circuit.append({'name': 'virtual_z', 'qubit': self.qubit, 'phase': np.pi/2})
            circuit.append({'name': 'read', 'qubit': self.qubit})
            circuits.append(circuit)
            circuit=[{'name': 'delay', 't': delaybeforecircuit}]
            circuit.append({'name': 'virtual_z', 'qubit': self.qubit, 'phase': -np.pi/2})
            circuit.extend(x90alpha)
            circuit.append({'name': 'virtual_z', 'qubit': self.qubit, 'phase': np.pi/2})
            circuit.extend(x90alpha)
            circuit.extend(x90alpha)
            circuit.append({'name': 'read', 'qubit': self.qubit})
            circuits.append(circuit)

        return circuits

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        self.output_dict= jobmanager.collect_all(self.circuits, num_shots_per_circuit, qchip=qchip)
        self.x90y180=self.output_dict['counts'].bitstring_dict['0'][0::2].flatten()
        self.y90x180=self.output_dict['counts'].bitstring_dict['0'][1::2].flatten()        

    @staticmethod
    def _line(x,a,b):
        return a*x+b

    def _fit_line(self,x,y):
        fit = curve_fit(self._line, x,y)
        return fit
            


    def fit(self):
        fitpara1=self._fit_line(x=self.alphas,y=self.x90y180)
        popt1=fitpara1[0]

        fitpara2=self._fit_line(x=self.alphas,y=self.y90x180)
        popt2=fitpara2[0]
        self.alphaoptm=(popt2[1]-popt1[1])/(popt1[0]-popt2[0])
        return self.alphaoptm
        
    def update_qchip(self, qchip,x90gate='X90'):
        qchip.gates[self.qubit + x90gate].contents[0].env[0]['paradict']['alpha'] = self.alphaoptm
