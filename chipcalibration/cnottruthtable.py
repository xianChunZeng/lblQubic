from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy
import sys


class cnottruthtable:
    """
    Define circuits, take data, and plot cnottruthtable
    """

    def __init__(self, control_qubit, target_qubit,delaybeforecircuit=600e-6,normorder='fro'):
        """
        Create cnottruthtable circuits according to input parameters, then compile to asm binaries.
        """
        self.control_qubit=control_qubit
        self.target_qubit=target_qubit
        self.delaybeforecircuit=delaybeforecircuit
        self.normorder=normorder

    def base_circuit(self,repeat,cnotgate='CNOT',twidth=None, amp=None, pztgt=None,pzctl=None,axtgt=None):
        self.cnotgate=cnotgate
        self.repeat=repeat
        cnotmodi={}
        if twidth is not None:
            cnotmodi.update({(1, 'twidth'): twidth})
            cnotmodi.update({(4, 't0'): twidth})
        if amp is not None:
            cnotmodi.update({(1, 'amp'): amp})
        if pztgt is not None:
            cnotmodi.update({(0, 'phase'): pztgt})
            cnotmodi.update({(2, 'phase'): -pztgt})
        if pzctl is not None:
            cnotmodi.update({(3, 'phase'): pzctl})
        if axtgt is not None:
            cnotmodi.update({(4, 'amp'): axtgt})

        qubitid=[self.control_qubit, self.target_qubit]
        circuit=[]
        for ctlinit in [0,1]:
            for tgtinit in [0,1]:
                circuit.append({'name': 'delay', 't': self.delaybeforecircuit})
                if ctlinit==1:
                    circuit.append({'name': 'barrier', 'qubit': qubitid})
                    circuit.append({'name': 'X90', 'qubit': qubitid[0]})
                    circuit.append({'name': 'X90', 'qubit': qubitid[0]})
                if tgtinit==1:
                    circuit.append({'name': 'barrier', 'qubit': qubitid})
                    circuit.append({'name': 'X90', 'qubit': qubitid[1]})
                    circuit.append({'name': 'X90', 'qubit': qubitid[1]})
                for i in range(repeat):
                    circuit.append({'name': 'barrier', 'qubit': qubitid})
                    circuit.append({'name': self.cnotgate, 'qubit': qubitid, 'modi':cnotmodi})
                    circuit.append({'name': 'barrier', 'qubit': qubitid})
                circuit.extend([{'name': 'read', 'qubit': qubit} for qubit in qubitid])
        return circuit
    def truthcircuit(self,repeat=1):
        self.circuits=[]
        self.circuits.append(self.base_circuit(repeat=repeat))
        self.x=[]
        self.reads_per_shot=4
    def scanaxtgt(self,axtgts,repeat=1):
        self.circuits=[]
        self.reads_per_shot=0
        self.x=axtgts
        circuit=[]
        for axtgt in axtgts:
            circuit.extend(self.base_circuit(repeat=repeat,axtgt=axtgt))
            self.reads_per_shot+=4
        self.circuits.append(circuit)
    def scanpzctl(self,pzctls,repeat=1):
        self.circuits=[]
        self.reads_per_shot=0
        self.x=pzctls
        circuit=[]
        for pzctl in pzctls:
            circuit.extend(self.base_circuit(repeat=repeat,pzctl=pzctl))
            self.reads_per_shot+=4
        self.circuits.append(circuit)
    def scanpztgt(self,pztgts,repeat=1):
        self.circuits=[]
        self.reads_per_shot=0
        self.x=pztgts
        circuit=[]
        for pztgt in pztgts:
            circuit.extend(self.base_circuit(repeat=repeat,pztgt=pztgt))
            self.reads_per_shot+=4
        self.circuits.append(circuit)
    def scanamp(self,amps,repeat=1):
        self.circuits=[]
        self.reads_per_shot=0
        self.x=amps
        circuit=[]
        for amp in amps:
            circuit.extend(self.base_circuit(repeat=repeat,amp=amp))
            self.reads_per_shot+=4
        self.circuits.append(circuit)

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        self.output_dict= jobmanager.collect_all(self.circuits, num_shots_per_circuit, qchip=qchip,reads_per_shot=self.reads_per_shot)
        qubits=self.output_dict['counts'].qubits
        ttlist=[]
        for ctlinit in [0,1]:
            for tgtinit in [0,1]:
                ttlist.append({self.control_qubit:str(ctlinit),self.target_qubit:str(tgtinit)})
        kklist=[]
        for tt in ttlist:
            k=''
            for q in qubits:
                k=k+tt[q]
            kklist.append(k)
#        print(kklist)                
        resultarray=[]
        for i in range(len(self.circuits)):
            countarray=numpy.array([self.output_dict['counts'].bitstring_dict[k][i] for k in kklist])
            total=numpy.sum(countarray,axis=0)
#            print(total)
            resultarray.append((countarray/total).T)
        self.tt=numpy.array(resultarray).reshape([-1,4,4])
        return self.tt
    def error(self,normorder=2):
        cnot=numpy.array([[1,0,0,0]
            ,[0,1,0,0]
            ,[0,0,0,1]
            ,[0,0,1,0]])
        iden=numpy.array([[1,0,0,0]
            ,[0,1,0,0]
            ,[0,0,1,0]
            ,[0,0,0,1]
            ])
        norm=[]
        for tt in self.tt:
            norm.append(numpy.linalg.norm((tt-(cnot if self.repeat%2==1 else iden)),ord=normorder ))
        return self.x,norm


    def fit(self):
        pass

    def update_qchip(self, qchip,twidth=None,amp=None,pztgt=None,pzctl=None,axtgt=None):
        if amp is not None:
            qchip.gates['%s%s%s'%(self.control_qubit,self.target_qubit,self.cnotgate)].contents[1].amp=amp
        pass
