from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic
from chipcalibration import blob
import itertools


class blob3:
    """
    Define circuits, take data, and plot blob3
    """

    def __init__(self, qubit, qchip):
        """
        Create blob3 circuits according to input parameters, then compile to asm binaries.
        """
        self.qubit = qubit
        self.qchip=qchip
    def blobamp(self,amps):
        self.circuits=[]
        circuit=[]
        for amp in amps:
            circuit.extend(self.base_circuit(amp=amp))
        self.circuits.append(circuit)
        self.x=amps
        self.reads_per_shot=3*len(self.x)
    def blobt0(self,t0s):
        self.circuits=[]
        circuit=[]
        for t0 in t0s:
            circuit.extend(self.base_circuit(t0=t0))
        self.circuits.append(circuit)
        self.x=t0s
        self.reads_per_shot=3*len(self.x)


    def blobfreq(self,freads=None,dfreads=None):
        if freads is None:
            if dfreads is not None:
                freads = self.qchip.qubits[self.qubit].readfreq+dfreads
            else:
                freads = [self.qchip.qubits[self.qubit].readfreq]
        else:
            freads = freads
        self.x=freads            
        self.reads_per_shot=3*len(self.x)
        self.circuits=[]
        circuit=[]
        for fread in freads:
            circuit.extend(self.base_circuit(fread=fread))
        self.circuits.append(circuit)

    def base_circuit(self,amp=None,t0=None,fread=None,delaybeforecircuit=600e-6):
        """
        Make list of circuits used for read freq measurement. 
        """
        modidict={}
        if fread is not None:
            modidict.update({(0,'freq'):fread,(1,'freq'):fread})
        if t0 is not None:            
            modidict.update({(1,'t0'):t0})
        if amp is not None:            
            modidict.update({(0,'amp'):amp})

        circuit=[]
        circuit.append({'name': 'delay', 't': delaybeforecircuit})
        circuit.append({'name': 'read', 'qubit': self.qubit, 'modi':modidict})
        circuit.append({'name': 'delay', 't': delaybeforecircuit})
        circuit.append({'name': 'X90', 'qubit': self.qubit })
        circuit.append({'name': 'X90', 'qubit': self.qubit })
        circuit.append({'name': 'read', 'qubit': self.qubit, 'modi':modidict})
        circuit.append({'name': 'delay', 't': delaybeforecircuit})
        circuit.append({'name': 'X90', 'qubit': self.qubit })
        circuit.append({'name': 'X90', 'qubit': self.qubit })
        circuit.append({'name': 'X90_ef', 'qubit': self.qubit })
        circuit.append({'name': 'X90_ef', 'qubit': self.qubit })
        circuit.append({'name': 'read', 'qubit': self.qubit, 'modi':modidict})


        return circuit

    def run_and_report(self, jobmanager, num_shots_per_circuit):
        self.raw_IQ= jobmanager.collect_raw_IQ(program_list=self.circuits, num_shots_per_circuit=num_shots_per_circuit
            ,reads_per_shot=self.reads_per_shot, delay_per_shot=0)
        self.dists=[]
        for chan,data in self.raw_IQ.items():
            for i,fread in enumerate(self.x):
                #                gmm = qubic.state_disc.GMMManager(chanmap_or_chan_cfgs=jobmanager.channel_configs)
                #gmm.fit({str(chan):data[:,:,i]})
                iq0=data[0,:,3*i+0]
                iq1=data[0,:,3*i+1]
                iq2=data[0,:,3*i+2]
                gmm = qubic.state_disc.GMMStateDiscriminator(n_states=3)
                gmm.fit(iqdata=data[:,:,i].flatten())
                gmm.set_labels_maxtomin(iqdata=iq0,labels_maxtomin=[0,1,2])
#                gmm.set_labels_maxtomin(iqdata=iq1,labels_maxtomin=['1'])
#                gmm.set_none_label('2')


                blob3={}
                distij={}
                for i in range(gmm.gmmfit.n_components):
                    blob3[i]=blob.c_blob(mean=gmm.gmmfit.means_[i],covar=gmm.gmmfit.covariances_[i])
                for i,j in itertools.combinations(blob3,2): 
                    distij[tuple(sorted([gmm.labels[i],gmm.labels[j]]))]=blob3[i].dist(blob3[j])
                self.dists.append(distij)                        
        return dict(x=self.x,dists=self.dists)

    def update_qchip(self, qchip,x90gate='X90'):
        #qchip.gates[self.qubit + x90gate].contents[0].env[0]['paradict']['alpha'] = self.alphaoptm
        pass
