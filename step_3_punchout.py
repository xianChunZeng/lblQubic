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

def run_punchout(qubitids, fbw, n_freq, atten_start, atten_stop, atten_step, n_samples, calirepo='submodules/qchip'):
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

    punchout = c_punchout(qubitid='vna', calirepo=calirepo)
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


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Run punchout freq x atten sweep')
    parser.add_argument('qubitids', nargs='+', help='list of qubit identifiers; e.g. Q0, Q1, etc.')
    parser.add_argument('--bandwidth', default=FBW, 
            help='frequency sweep bandwidth in Hz, default {}'.format(FBW))
    parser.add_argument('--n-freq', default=N_FREQ, help='frequency sweep bandwidth in Hz, default {}'.format(N_FREQ))
    parser.add_argument('--atten-start', default=ATTEN_START, 
            help='starting (high) attenuation, default{}. Note that atten values are negative, so lower value means lower power'.format(ATTEN_START))
    parser.add_argument('--atten-stop', default=ATTEN_STOP, 
            help='ending (low) attenuation, default{}'.format(ATTEN_STOP))
    parser.add_argument('--atten-step', default=ATTEN_STEP, 
            help='dB increment to use in atten sweep, default{}'.format(ATTEN_STEP))
    parser.add_argument('--n-samples', default=N_SAMPLES, 
            help='number of samples in readout acquisition buffer, default{}'.format(N_SAMPLES))
    args = parser.parse_args()

    run_punchout(args.qubitids, args.bandwidth, args.n_freq, args.atten_start, args.atten_stop, args.atten_step, args.n_samples)
