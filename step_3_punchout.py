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
import argparse
from qubic.qcvv.punchout import c_punchout

FBW = 6e6
N_FREQ = 200
ATTEN_START = -35
ATTEN_STOP = 0.2
ATTEN_STEP = 5.0
N_SAMPLES = 103

class PunchoutGUI:
    """
    Implements clickable GUI for selecting resonator power/freq
    """

    def __init__(self, punchout, qubitid):
        """
        Parameters
        ----------
            punchout : c_punchout object
                object containing punchout data to be plotted
            qubitid : str
                qubit identifier (todo: get this from punchout)
        """
        self.fig1=plt.figure(figsize=(10, 10))
        self.sub=self.fig1.subplots(2, 2)
        self.fig1.suptitle(qubitid)
        punchout.plotamp(self.sub[0, 0])
        punchout.plotang(self.sub[0, 1])
        punchout.plotampdiff(self.sub[1, 0])
        punchout.plotangdiff(self.sub[1, 1])
        self.fig1.canvas.mpl_connect('button_press_event', self.onClick)
        print('Click any plot to select desired resonator attenuation and frequency')
        plt.show()

    def onClick(self, event):
        self.freq = event.xdata
        self.atten = event.ydata
        print('Selected resonator frequency {} and attenutation {}'.format(self.freq, self.atten))
        print('Click again to change, otherwise close')

def run_punchout(qubitid, fbw, n_freq, atten_start, atten_stop, atten_step, n_samples):
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

    punchout=c_punchout(qubitid='vna', calirepo='submodules/qchip')
    fcenter=punchout.opts['qchip'].getfreq(qubitid+'.readfreq')
    fstart=fcenter-fbw/2
    fstop=fcenter+fbw/2
    fx=np.linspace(fstart, fstop, 200)
    
    punchout.run(n_samples, fx=fx, attens=np.arange(atten_start, atten_stop, atten_step), maxvatatten=0)
    
    fig2=plt.figure()
    sub2=fig2.subplots()
    punchout.plotmaxangdiff(sub2)
    fig2.suptitle(qubitid)
    fig3=plt.figure()
    ax3=fig3.subplots(2, 2)
    ax3[0, 0].plot(punchout.s11.real.T, punchout.s11.imag.T)
    ax3[1, 0].plot(abs(punchout.s11.T))
    ax3[1, 1].plot(np.angle(punchout.s11.T))

    cal_gui = PunchoutGUI(punchout, qubitid)
    freq = cal_gui.freq
    atten = cal_gui.atten

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Run punchout freq x atten sweep')
    parser.add_argument('qubitid', help='qubit identifier; e.g. Q0, Q1, etc.')
    args = parser.parse_args()

    run_punchout(args.qubitid, FBW, N_FREQ, ATTEN_START, ATTEN_STOP, ATTEN_STEP, N_SAMPLES)
