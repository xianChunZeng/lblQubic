"""
Script for "punching out" qubit; i.e. initial estimates of readout 
resonator drive attenuation and frequency.

TODO:
    - integrate with cfg
    - see which params should be CL args
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
from qubic.qcvv.punchout import c_punchout

FBW = 6e6 
N_FREQ = 200
ATTEN_START = -35
ATTEN_STOP = 0.2
ATTEN_STEP = 5.0
N_SAMPLES = 100

class PunchoutGUI:
    """
    Implements clickable GUI for selecting resonator power/freq
    """

    def __init__(self, punchout, sweep_index=None, qubitid=''):
        """
        Parameters
        ----------
            punchout : c_punchout object
                object containing punchout data to be plotted
            qubitid : str
                qubit identifier (todo: get this from punchout)
        """
        self.fig1 = plt.figure(figsize=(10, 10))
        self.sub = self.fig1.subplots(2, 1)
        self.fig1.suptitle(qubitid)
        punchout.plotamp(self.sub[0], sweep_index)
        punchout.plotang(self.sub[1], sweep_index)
        self.fig1.canvas.mpl_connect('button_press_event', self.onClick)
        print('Click any plot to select desired resonator attenuation and frequency')
        plt.show()

    def onClick(self, event):
        self.freq = event.xdata
        self.atten = event.ydata
        print('Selected resonator frequency {} and attenutation {}'.format(self.freq, self.atten))
        print('Click again to change, otherwise close')

def run_punchout(qubit_dict, qchip, inst_cfg, fbw=FBW, n_freq=N_FREQ, atten_start=ATTEN_START, \
            atten_stop=ATTEN_STOP, atten_step=ATTEN_STEP, n_samples=N_SAMPLES):
    """
    Runs punchout sweep on selected qubit, plots results in clickable GUI. 
    TODO: this should also update configs.

    Parameters
    ----------
        qubitid : str
            qubit identifier
        fbw : float
            sweep bandwidth around initial res freq estimate
        n_freq : int
            number of freq points in sweep
        atten_start, atten_stop, atten_step : int
            parameters for atten sweep
        n_samples : int
            ??
    """

    punchout = c_punchout(qubitid='vna', qchip=qchip, instrument_cfg=inst_cfg)
    qubitids = np.asarray(qubitids)

    fx = np.empty((0, n_freq))
    for qubitid in qubitids:
        fcenter = punchout.opts['qchip'].getfreq(qubitid+'.readfreq')
        fstart = fcenter - fbw/2
        fstop = fcenter + fbw/2
        fx = np.vstack((fx, np.linspace(fstart, fstop, n_freq)))
    
    if fx.shape[0] > 1:
        inds = np.argsort(fx[:,0])
        fx = fx[inds]
        qubitids = qubitids[inds]
    #fx = np.squeeze(fx) #remove first axis if there's only one qubit
    punchout.run(n_samples, fx=fx, attens=np.arange(atten_start, atten_stop, atten_step), maxvatatten=0)
     
    freq = []
    atten = []
    for i in range(len(qubitids)):
        cal_gui = PunchoutGUI(punchout, i, qubitids[i])
        freq.append(cal_gui.freq)
        atten.append(cal_gui.atten)

    return freq, atten

def get_qubit_dict(qubitids, qchip):
    """
    Helper function for returning punchout
    qubit dict from qchip config object.

    Parameters
    ----------
        qubitids : list
            list of qubit IDs to punch out ('Q0', 'Q1', etc)
        qchip : c_qchip object
            object containing qubit calibration config
    Returns
    -------
        dict
            keys: qubitids (from input)
            values: readout frequencies corresponding to each 
                    qubit, as specified by qchip
    """
    qubits = qchip.paradict['Qubits']
    qubitdict = {}
    for k, v in qubits.items():
        if k in qubitids:
            qubitdict[k] = v['readfreq']

    return qubitdict

