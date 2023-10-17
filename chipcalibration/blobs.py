from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
import sys
import qubic
from chipcalibration import blob
import itertools


class blobfreq:
    """
    Define circuits, take data, and plot blobfreq
    """

    def __init__(self, qubit, qchip, freads=None,dfreads=None):
        """
        Create blobfreq circuits according to input parameters, then compile to asm binaries.
        """
        self.qubit = qubit
        self.qchip=qchip
        if freads is None:
            if dfreads is not None:
                self.freads = qchip.qubits[qubit].readfreq+dfreads
            else:
                self.freads = [qchip.qubits[qubit].readfreq]
        else:
            self.freads = freads
        self.circuits = self._make_blobfreq_circuits()

    def _make_blobfreq_circuits(self,delaybeforecircuit=600e-6):
        """
        Make list of circuits used for read freq measurement. 
        """

        circuits=[]
        circuit=[]
        for fread in self.freads:
            circuit.append({'name': 'delay', 't': delaybeforecircuit})
            circuit.append({'name': 'X90', 'qubit': self.qubit })
            circuit.append({'name': 'read', 'qubit': self.qubit, 'modi':{(0,'freq'):fread,(1,'freq'):fread}})
        circuits.append(circuit)

        return circuits

    def run_and_report(self, jobmanager, num_shots_per_circuit):
        self.raw_IQ= jobmanager.collect_raw_IQ(program_list=self.circuits, num_shots_per_circuit=num_shots_per_circuit
            ,reads_per_shot=len(self.freads), delay_per_shot=0, qchip=self.qchip)
        self.dists=[]
        for chan,data in self.raw_IQ.items():
            for i,fread in enumerate(self.freads):
                gmm = qubic.state_disc.GMMManager(chanmap_or_chan_cfgs=jobmanager.channel_configs)
                gmm.fit({str(chan):data[:,:,i]})
                blobs={}
                distij={}
                gmmfit=gmm.gmm_dict[self.qubit].gmmfit
                for i in range(gmmfit.n_components):
                    blobs[i]=blob.c_blob(mean=gmmfit.means_[i],covar=gmmfit.covariances_[i])
                for i,j in itertools.combinations(blobs,2): 
                    distij[tuple(sorted([i,j]))]=blobs[i].dist(blobs[j])
                self.dists.append(distij)                        
        return dict(freads=self.freads,dists=self.dists)

    def update_qchip(self, qchip,x90gate='X90'):
        qchip.gates[self.qubit + x90gate].contents[0].env[0]['paradict']['alpha'] = self.alphaoptm
