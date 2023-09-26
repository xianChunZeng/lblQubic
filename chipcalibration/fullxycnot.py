from scipy.optimize import curve_fit
from matplotlib import pyplot 
import numpy
import sys
import chipcalibration.crfit as crfit
from scipy.optimize import curve_fit

def zxzxz(qubitid,zprep):
    if isinstance(qubitid,str):
        qubitid=[qubitid]
    onequbit=len(qubitid)==1
    seq=[]
    for iz,z in enumerate(zprep):
        if iz!=0:
            seq.extend([{'name': 'X90', 'qubit': qubit} for qubit in qubitid])
        seq.extend([{'name': 'virtual_z', 'qubit': qubit, 'phase': z if onequbit else z[iq]} for iq,qubit in enumerate(qubitid)])
    return seq

def baseseqs(self,qubitid,xyrot,zprep,**kwargs):
        #    ops=dict(elementlength=80,elementstep=4e-9)
        #ops.update(**kwargs)
        seqs=[]
        for rot in xyrot:
            #        for tomo in [tomox,tomoy,tomoz]:
            for tomo in ['X90']:#,'Y-90',None]:
                seq=[]
                seq.extend(zxzxz.zxzxz(qubitid=qubitid,zprep=zprep))
                seq.append({'name': 'CNOT', 'qubit': qubitid})
                seq.extend([{'name': 'virtualz', 'qubit': qubit, 'para': {'phase': rot}} for qubit in qubitid])
                seq.append({'name': 'barrier', 'qubit': qubitid})
                if tomo is not None:
                    seq.extend([{'name': tomo, 'qubit': qubit} for qubit in qubitid])
                seq.extend([{'name': 'read', 'qubit': qubit} for qubit in qubitid])
                seqs.append(seq)
        return seqs



class fullxycnot:
    """
    Define circuits, take data, and plot fullxycnot fit for cnot finder, implement the 2 qubit calibration protocol as described in https://dl.acm.org/doi/full/10.1145/3529397
    """

    def __init__(self, control_qubit, target_qubit, xyrot,zprep,axes=('X', 'Y', 'Z')):
        """
        Create fullxycnot circuits according to input parameters, then compile to asm binaries.
        """
        self.control_qubit = control_qubit
        self.target_qubit= target_qubit
        self.xyrot=xyrot
        self.zprep=zprep
        self.axes=axes
        self.circuits = self._make_fullxycnot_circuits()

    def _make_fullxycnot_circuits(self,delaybeforecircuit=600e-6):
        """
        Make list of circuits used for drag alpha measurement. 
        """
        qubitid=[self.control_qubit, self.target_qubit]
        circuits = []
        circuit=[]
        for rot in self.xyrot:
            for axis in self.axes:
                circuit.append({'name': 'delay', 't': delaybeforecircuit})
                circuit.extend(zxzxz(qubitid=qubitid,zprep=self.zprep))
                circuit.append({'name': 'barrier'})
                
                circuit.append({'name': 'CNOT', 'qubit': qubitid})
                circuit.extend([{'name': 'virtual_z', 'qubit': qubit, 'phase': rot} for qubit in qubitid])
                
                circuit.append({'name': 'barrier'})
                if axis == 'X':
                    circuit.append({'name': 'Y-90', 'qubit': [self.control_qubit]})
                    circuit.append({'name': 'Y-90', 'qubit': [self.target_qubit]})    
                elif axis == 'Y':
                    circuit.append({'name': 'X90', 'qubit': [self.control_qubit]})
                    circuit.append({'name': 'X90', 'qubit': [self.target_qubit]})
                    
                circuit.append({'name': 'barrier'})
                circuit.extend([{'name': 'read', 'qubit': qubit} for qubit in qubitid])
        circuits.append(circuit)            
        
        return circuits    


    def run_and_report(self, jobmanager, num_shots_per_circuit, qchip):
        self.results_out = jobmanager.build_and_run_circuits(self.circuits, num_shots_per_circuit, outputs=['shots', 'counts'],qchip=qchip,reads_per_shot=len(self.xyrot))
        self.result_dict=self.results_out['counts'].bitstring_dict

    def plot(self):
        for bitstr in ['00','01','10','11']:
            pyplot.plot(numpy.rad2deg(xyrot),self.result_dict.flatten(),label=bitstr)
        pyplot.legend()

    def fit(self,plot=True):
        xymeas4=numpy.zeros((4,len(self.xyrot)))
        total=numpy.sum(numpy.array([v for k,v in self.result_dict.items()]),axis=0)
        total=total[0][0]
        for i,k in enumerate(['00','10','01','11']):
            xymeas4[i,:]=self.result_dict[k][0]
            if plot:
                pyplot.figure('raw')
                pyplot.plot(xymeas4[i,:],label=k)
        ydata=xymeas4.T.flatten('F')/total
        xdata=numpy.tile(self.xyrot.T,4)#/numpy.pi*180
        fun_to_fit = lambda xyscan, pzctl, pztgt, pxtgt: crfit.fitfun(xyscan, pzctl, pztgt, pxtgt,zzdeglist=self.zprep*180/numpy.pi)
        popt,pcov=curve_fit(fun_to_fit,xdata=xdata,ydata=ydata,p0=(0,0,0),bounds=(([-2*numpy.pi,-2*numpy.pi,-2*numpy.pi],[2*numpy.pi,2*numpy.pi,2*numpy.pi])))
        xyfit=crfit.fitfun(xdata,*popt)
        pyplot.figure('x4fit')
        pyplot.plot(ydata,'*')
        pyplot.plot(xyfit)
        self.pzctl=popt[0]
        self.pztgt=popt[1]
        self.pxtgt=popt[2]
        return dict(pzctl=popt[0],pztgt=popt[1],pxtgt=popt[2])
        
    def update_qchip(self, qchip):
        qchip.gates['Q1Q0CNOT'].contents[0]._phase+=self.pztgt
        qchip.gates['Q1Q0CNOT'].contents[2]._phase-=self.pztgt
        qchip.gates['Q1Q0CNOT'].contents[3]._phase+=self.pzctl
        x90amp=qchip.gates['Q0X90'].contents[0].amp
        origamp=qchip.gates['Q1Q0CNOT'].contents[4].amp
        ampperrad=x90amp*2/numpy.pi
        newamp=((origamp/ampperrad)+self.pxtgt)%(numpy.pi*2)*ampperrad
        qchip.gates['Q1Q0CNOT'].contents[4].amp=newamp
        return qchip.gates['Q1Q0CNOT'].cfg_dict




