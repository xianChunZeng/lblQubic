from scipy.optimize import curve_fit
from matplotlib import pyplot
import numpy 
import sys
import chipcalibration.fit as fit


class rabi:
    """
    Define circuits, take data, and plot rabi
    """

    def __init__(self, qubit,qchip,rabigate='rabi3_cos'):
        """
        Create rabi circuits according to input parameters, then compile to asm binaries.
        """
        self.qubit = qubit
        self.qchip=qchip
        self.rabigate=rabigate
    def rabiamp(self, amps):
        self.circuits=[]
        circuit=[]
        for amp in amps:
            circuit.extend(self.rabi_circuit(amp=amp))
        self.circuits.append(circuit)
        self.x=amps
        self.reads_per_shot=len(amps)

    def rabitime(self, twidths):        
        self.circuits=[]
        circuit=[]
        for twidth in twidths:
            circuit.extend(self.rabi_circuit(twidth=twidth))
        self.circuits.append(circuit)
        self.x=twidths
        self.reads_per_shot=len(twidths)

    def rabi_circuit(self,amp=None,twidth=None,delaybeforecircuit=600e-6):
        """
        Make list of circuits used for drag alpha measurement. 
        """
        rabimodi={}
        if amp is not None:
            rabimodi.update({(0, 'amp'): amp,(1,'amp'):amp,(2,'amp'):amp})
        if twidth is not None:
            rising=self.qchip.gates['%s%s'%(self.qubit,self.rabigate)].contents[0].twidth
            faling=self.qchip.gates['%s%s'%(self.qubit,self.rabigate)].contents[2].twidth
            rabimodi.update({(2,'t0'):twidth-faling})
            print(twidth-faling)

        circuit=[]
        circuit.append({'name': 'delay', 't': delaybeforecircuit, 'qubit': self.qubit})
        circuit.append({'name': 'barrier', 'qubit': self.qubit})
        circuit.append({'name': self.rabigate, 'qubit': self.qubit ,'modi':rabimodi})
        circuit.append({'name': 'barrier', 'qubit': self.qubit})
        circuit.append({'name': 'read', 'qubit': self.qubit})

        return circuit

    def run_and_report(self, jobmanager, num_shots_per_circuit):
        self.output_dict= jobmanager.collect_all(program_list=self.circuits, num_shots_per_circuit=num_shots_per_circuit, reads_per_shot=self.reads_per_shot, qchip=self.qchip)
        self.raw_IQ=self.output_dict['s11']
        predict=jobmanager.gmm_manager.predict(self.raw_IQ)
        self.p1=numpy.mean(predict[self.qubit],axis=1).flatten()
        return self.p1


    def fitexpsin(self):
        estsine=fit.sinestimate(x=self.delay_interval,y=self.p1)
        est=list(estsine)
        est.append(self.delay_interval[-1])
        expsinfit=fit.fitexpsin(x=self.delay_interval,y=self.p1,p0=est,bounds=(([0,-numpy.inf,0,-1,1e-9],[2,numpy.inf,2*numpy.pi,2,numpy.inf])))
        return expsinfit

    

    def update_qchip(self, qchip, plusorminus=None):
        if plusorminus is None:
            plusorminus=input('use the plus(p) or minus(m) or original(o) in the qubitcfg.json? or quit(q)')
        while plusorminus.lower() not in ['p','m','plus','minus','o','orig','original','q','quit']:
            plusorminus=input('use the plus(p) or minus(m) or original(o) in the qubitcfg.json? or quit(q)')
        if plusorminus.lower() in ['p','plus']:            
            qchip.qubits[self.qubit].freq += self.ffit
        elif plusorminus.lower() in ['m','minus']:            
            qchip.qubits[self.qubit].freq -= self.ffit
        else:
            pass
