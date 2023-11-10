from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy
import sys


class tcr:
    """
    Define circuits, take data, and plot tcr
    """

    def __init__(self, control_qubit, target_qubit, delaybeforecircuit=600e-6):
        
        """
        Create tcr circuits according to input parameters, then compile to asm binaries.
        """
        self.control_qubit=control_qubit
        self.target_qubit=target_qubit
        self.delaybeforecircuit=delaybeforecircuit

    def cr_amp_sweep(self,amps,repeat=1,crcnotgate='CR'):
        self.circuits=[]
        circuit=[]
        for amp in amps:
            circuit.extend(self.base_circuit(amp=amp,repeat=repeat,crcnotgate=crcnotgate))
        self.circuits.append(circuit)
        self.reads_per_shot=len(amps)*6
        self.x=amps
    def cr_amp_width_sweep(self,twidths, amps,repeat=1,crcnotgate='CR'):
        self.circuits=[]
        for twidth in twidths:
            circuit=[]
            for amp in amps:
                circuit.extend(self.base_circuit(twidth=twidth, amp=amp,repeat=repeat,crcnotgate=crcnotgate))
            self.circuits.append(circuit)                
        self.reads_per_shot=len(amps)*6
        self.x=twidths
    def cnot_sweep_pztgt(self,pztgts,repeat=1,crcnotgate='CNOT'):
        self.circuits=[]
        circuit=[]
        for pztgt in pztgts:
            circuit.extend(self.base_circuit(pztgt= pztgt,repeat=repeat,crcnotgate=crcnotgate))
        self.circuits.append(circuit)        
        self.reads_per_shot=len(pztgts)*6
        self.x=pztgts
        
#plt.legend()    

    def base_circuit(self,repeat,crcnotgate,twidth=None, amp=None, pztgt=None,pzctl=None,axtgt=None):
        circuit = []
        index=0
        for ctrl_state in [0,1]:
            for axis in ('X', 'Y', 'Z'):
                crcnotmodi={}
                if crcnotgate=='CR':
                    if twidth is not None:
                        crcnotmodi.update({(0, 'twidth'): twidth})
                    if amp is not None:
                        crcnotmodi.update({(0, 'amp'): amp})
                elif crcnotgate=='CNOT':
                    if twidth is not None:
                        crcnotmodi.update({(1, 'twidth'): twidth})
                        crcnotmodi.update({(4, 't0'): twidth})
                    if amp is not None:
                        crcnotmodi.update({(1, 'amp'): amp})
                    if pztgt is not None:
                        crcnotmodi.update({(0, 'phase'): pztgt})
                        crcnotmodi.update({(2, 'phase'): -pztgt})
                    if pzctl is not None:
                        crcnotmodi.update({(3, 'phase'): pzctl})
                    if axtgt is not None:
                        crcnotmodi.update({(4, 'amp'): axtgt})
                        
                circuit.append({'name': 'delay', 't': self.delaybeforecircuit})
                if ctrl_state == 1:
                    circuit.append({'name': 'X90', 'qubit':[self.control_qubit]})
                    circuit.append({'name': 'X90', 'qubit':[self.control_qubit]})
                circuit.append({'name': 'barrier'})
                for i in range(repeat):
                    circuit.append({'name': crcnotgate, 'qubit': [self.control_qubit, self.target_qubit],
                               'modi':crcnotmodi})
                if axis == 'X':
                    circuit.append({'name': 'Y-90', 'qubit': [self.target_qubit]})
                elif axis == 'Y':
                    circuit.append({'name': 'X90', 'qubit': [self.target_qubit]})
                circuit.append({'name': 'barrier'})
                circuit.append({'name': 'read', 'qubit': [self.control_qubit]})
                circuit.append({'name': 'read', 'qubit': [self.target_qubit]})
        return circuit

    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        self.output_dict= jobmanager.collect_all(self.circuits, num_shots_per_circuit, qchip=qchip,reads_per_shot=self.reads_per_shot)
        self.c0={}
        self.c1={}
        self.r={}
        for k,v in self.output_dict['shots'].items():
            v=v.astype(numpy.float64)
            vavr=numpy.average(v,axis=1)
            ss=vavr.reshape((len(self.circuits),self.reads_per_shot//6,6))
            self.c0[k]=ss[:,:,:3]*2-1
            self.c1[k]=ss[:,:,3:]*2-1
            self.r[k]=0.5*numpy.sqrt(numpy.sum(numpy.square(self.c0[k]-self.c1[k]),axis=2))
#            print('v',v.shape)
#            print('vavr',vavr.shape)
#            print('ss',ss.shape)
#            print('c0',self.c0[k].shape)
#            print('c1',self.c1[k].shape)
#            print('r',self.r[k].shape)
            
        return self.c0,self.c1,self.r
    def plotxyzr(self,figure):
        ax=figure.subplots(4)
        for itomo,tomo in enumerate(['x','y','z']):
            for iamp,amp in enumerate(amps):
                ax[itomo].plot(self.x,c0['Q0'][:,iamp,itomo],label='%s %8.2f'%('c0',amp))
                ax[itomo].plot(self.x,c1['Q0'][:,iamp,itomo],label='%s %8.2f'%('c1',amp))
                ax[itomo].set_ylabel(tomo)
                ax[itomo].set_ylim((-1,1))
                #ax[itomo].legend()
        for iamp,amp in enumerate(amps):
            ax[3].plot(twidths,r['Q0'][:,iamp],'*-',label='%s %8.2f'%('r',amp))
#ax[3].legend()
        ax[3].set_ylim((0,1))



    def fit(self):
        fitpara1=self._fit_line(x=self.alphas,y=self.x90y180)
        popt1=fitpara1[0]

        fitpara2=self._fit_line(x=self.alphas,y=self.y90x180)
        popt2=fitpara2[0]
        self.alphaoptm=(popt2[1]-popt1[1])/(popt1[0]-popt2[0])
        return self.alphaoptm
        
    def update_qchip(self, qchip,x90gate='X90'):
        qchip.gates[self.qubit + x90gate].contents[0].env[0]['paradict']['alpha'] = self.alphaoptm
