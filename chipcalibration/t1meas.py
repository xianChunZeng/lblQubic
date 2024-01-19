import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic.toolchain as tc
from qubic.state_disc import GMMManager
import chipcalibration.fit as fit
from scipy.optimize import curve_fit


ACC_BUFSIZE = 1000

class t1meas:
    """
    Define circuits, take data, and plot t1meas
    """

    def __init__(self, qubits, elementlength, elementstep, qchip, fpga_config, channel_configs):
        """
        Create t1meas circuits according to input parameters, then compile to asm binaries.
        """
        self.circuits = self._make_t1meas_circuits(qubits, elementlength, elementstep, qchip)
        self.qubits = qubits
        self.elementlength=elementlength
        self.elementstep=elementstep
        self.qchip=qchip
        self.readout_chanmap = {qubit: channel_configs[qubit + '.rdlo'].core_ind for qubit in qubits}
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        compiled_progs = tc.run_compile_stage(self.circuits, fpga_config, qchip)
        self.raw_asm_progs = tc.run_assemble_stage(compiled_progs, channel_configs)

    def _make_t1meas_circuits(self, qubits, elementlength, elementstep, qchip):
        """
        Make list of circuits used for chevron pattern measurement. Each circuit covers the 
        full range of frequencies, and a single pulse width. So there will be a total of 
        len(pulse_widths) circuits, each of which contains nfreq measurements. A 400 us 
        delay is inserted between each measurement.
        """
        circuits = []
        cur_circ = []
        qubit=qubits[0]
        
        for i in range(elementlength):
            #for qubit in qubits:
            cur_circ.append({'name': 'delay', 't': 500.e-6})
            cur_circ.append({'name': 'X90', 'qubit': [qubit]})
            cur_circ.append({'name': 'X90', 'qubit': [qubit]})
            cur_circ.append({'name': 'delay', 't': elementstep*(i)})
            cur_circ.append({'name': 'read', 'qubit': [qubit]})
        circuits.append(cur_circ)
        return circuits
    
    def run_and_report(self, jobmanager, num_shots_per_circuit):
        self.reads_per_shot = self.elementlength
        self.output_dict= jobmanager.collect_all(program_list=self.circuits, num_shots_per_circuit=num_shots_per_circuit, reads_per_shot=self.reads_per_shot, qchip=self.qchip)
        self.raw_IQ=self.output_dict['s11']
        predict=jobmanager.gmm_manager.predict(self.raw_IQ)
        self.p1=np.mean(predict[self.qubits[0]],axis=1).flatten()
        return self.p1

    def run(self, circuit_runner, nsamples):
        """
        Run the chevron circuit with nsamples shots per freq x pulse_width point.

        Parameters
        ----------
            circuit_runner : CircuitRunner object
            nsamples : int
        """
        reads_per_shot = self.elementlength
        nshot = ACC_BUFSIZE//reads_per_shot
        navg = int(np.ceil(nsamples/nshot))
        self.s11 = {chan: np.zeros((self.elementlength, nshot*navg), dtype=np.complex128) 
                    for chan in self.readout_chanmap.values()}

        print('nshot',nshot)        
        print('navg',navg)        
        for k,v in self.s11.items():
            print(k,v.shape)
        
        # circuit_runner.load_circuit(self.raw_asm_progs[0])
        shots = circuit_runner.run_circuit_batch(self.raw_asm_progs, nsamples, 1, delay_per_shot=500.e-6)#circuit_runner.run_circuit_(nshot, navg, reads_per_shot, delay=0.5)
                
        for chan, iqdata in shots.items():
            self.s11[chan] = iqdata.reshape((nshot*navg, self.elementlength)).T
            
#        for i, raw_asm in enumerate(self.raw_asm_progs):
#            circuit_runner.load_circuit(raw_asm)
#            shots = circuit_runner.run_circuit(nshot, navg, reads_per_shot, delay=0.5)
#            for chan, iqdata in shots.items():
#                #self.s11[chan][i] = np.reshape(iqdata, (nshot*navg, len(self.freqoffsets))).T
#                self.s11[chan][i] = np.reshape(iqdata, (nshot*navg, self.elementlength)).T

        self._fit_gmm()

    def _fit_gmm(self):
        self.gmm_manager.fit(self.s11)
        #self.gmm_manager.set_labels_maxtomin({chan: shots[np.argmin(self.elementlength)].flatten() 
        #                                      for chan, shots in self.s11.items()}, [0, 1])
        self.state_disc_shots = self.gmm_manager.predict(self.s11)
        self.ones_frac = {qubit: np.sum(self.state_disc_shots[qubit], axis=1) for qubit in self.state_disc_shots.keys()}
        self.zeros_frac = {qubit: np.sum(self.state_disc_shots[qubit] == 0, axis=1) for qubit in self.state_disc_shots.keys()}

    def fitT1(self):
        # estsine=fit.sinestimate(x=self.delay_interval,y=self.p1)
        # est=list(estsine)
        # est.append(self.delay_interval[-1])
        t=np.arange((self.elementlength))*self.elementstep
        expfit=fit.fitexp(x=t,y=self.p1, bounds=((0,np.inf)))
        return expfit
    

    # def fitt1meas(self,dt,data,plot=True,figname='fitt1meas'):
    #     tdata = dt * np.arange(len(data)) # unit: seconds
    #     popt, pcov = curve_fit(expfit, tdata, data, bounds=((0,np.inf)))
    #     t1fit = 1./popt[1]
    #     fiterr= np.sqrt(np.diag(pcov))
    #     if plot:
    #         pyplot.figure(figname)
    #         pyplot.plot(tdata,expfit(tdata,*popt))
    #         pyplot.plot(tdata,data,'*')
    #         pyplot.ylim(0,1)
    #         pyplot.savefig(figname+'.pdf')
    #     return [t1fit,fiterr]