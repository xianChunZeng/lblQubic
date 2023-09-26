from scipy.optimize import curve_fit
from matplotlib import pyplot
import numpy 
import sys
import chipcalibration.fit as fit


class ramsey:
    """
    Define circuits, take data, and plot ramsey
    """

    def __init__(self, qubit,delay_interval, qchip):
        """
        Create ramsey circuits according to input parameters, then compile to asm binaries.
        """
        self.qubit = qubit
        self.delay_interval=delay_interval
        self.qchip=qchip
    def ramsey(self, framsey_offset=None):        
        self.circuits = self._make_ramsey_circuits(framsey_offset=framsey_offset)

    def _make_ramsey_circuits(self,framsey_offset,delaybeforecircuit=600e-6):
        """
        Make list of circuits used for drag alpha measurement. 
        """
        if framsey_offset is None:
            fqdrv=self.qchip.qubits[self.qubit].freq 
        else:
            fqdrv=self.qchip.qubits[self.qubit].freq + framsey_offset

        circuits= []
        circuit=[]
        for tdelay in self.delay_interval:
            circuit.append({'name': 'delay', 't': delaybeforecircuit, 'qubit': self.qubit})
            circuit.append({'name': 'X90', 'qubit': self.qubit ,'modi':{(0,'freq'):fqdrv}})
            circuit.append({'name': 'delay', 't': tdelay, 'qubit': self.qubit})
            circuit.append({'name': 'X90', 'qubit': self.qubit,'modi':{(0,'freq'):fqdrv}})
            circuit.append({'name': 'read', 'qubit': self.qubit})
        circuits.append(circuit)

        return circuits

    def run_and_report(self, jobmanager, num_shots_per_circuit):
        self.output_dict= jobmanager.collect_all(program_list=self.circuits, num_shots_per_circuit=num_shots_per_circuit, reads_per_shot=len(self.delay_interval), qchip=self.qchip)
        self.raw_IQ=self.output_dict['s11']
        predict=jobmanager.gmm_manager.predict(self.raw_IQ)
        self.p1=numpy.mean(predict[self.qubit],axis=1).flatten()
        return self.p1

    def scanoffsets(self,framsey_offsets, jobmanager, num_shots_per_circuit):
        p1offset={}
        ffitoffset={}
        for framsey_offset in framsey_offsets:
            self.circuits = self._make_ramsey_circuits(framsey_offset=framsey_offset)
            p1offset[framsey_offset]=self.run_and_report(jobmanager, num_shots_per_circuit)
            expsinfit=self.fitexpsin()
            ffitoffset[framsey_offset]=expsinfit['popt'][1]
        return dict(p1offset=p1offset,ffitoffset=ffitoffset)

    def plusminus(self,jobmanager, num_shots_per_circuit,plot=True):
        self.circuits = self._make_ramsey_circuits(framsey_offset=0)
        p1orig=self.run_and_report(jobmanager, num_shots_per_circuit)
        expsinfit=self.fitexpsin()
        self.ffit=expsinfit['popt'][1]
        pmdict=self.scanoffsets(framsey_offsets=[self.ffit,-self.ffit], jobmanager=jobmanager, num_shots_per_circuit=num_shots_per_circuit)
        if plot:
            pyplot.plot(self.delay_interval,p1orig,label='original')
            for k,v in pmdict['p1offset'].items():
                pyplot.plot(self.delay_interval,v,label='plus' if k==self.ffit else 'minus')
            pyplot.ylim((0,1))
            pyplot.xlabel('time(s)')
            pyplot.ylabel('P<1>')
            pyplot.title(self.qubit)
            pyplot.legend()





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


def absx(x,x0):
    return abs(x-x0)
def ramsey_optimize(qubit,delay_interval, qchip, framsey_offsets,jobmanager,num_shots_per_circuit,plot=True):
    iramsey=ramsey(qubit=qubit,delay_interval=delay_interval,qchip=qchip)
    dscan=iramsey.scanoffsets(framsey_offsets=framsey_offsets,jobmanager=jobmanager,num_shots_per_circuit=num_shots_per_circuit)
    x=[]
    y=[]
    for k,v in dscan['ffitoffset'].items():
        x.append(k)
        y.append(v)
    absfit=curve_fit(absx,x,y)        
    return dict(x=x,y=y,absfit=absfit)

