from matplotlib import pyplot as plt
import numpy as np
import scipy.signal as signal
from qubic.qcvv.vna import c_vna
from step_1_alignment import TLO

class VNAClickGUI:

    def __init__(self, freqs, phases, orig_qubitdict=None, snaprad=5):
        """
        Parameters
        ----------
            freqs : numpy array 
                points sweeped by VNA
            phases : numpy array
                unwrapped, detrened arg(S_11)
            orig_qubitdict : dict
                dictionary w/ entries {'Qn': freq} corresponding to 
                resonator frequencies from a previous qubit calibration
            snaprad : int
                number of points to snap peak to when clicking

        """
        self.freqs = freqs[:-1] 
        self.phasediffs = np.diff(phases)
        self.snaprad = snaprad
        self.peak_inds = find_peaks_phasediff(phases)
        self.orig_qubitdict = orig_qubitdict

        self.fig = plt.figure(figsize=(20,8))
        self.ax = self.fig.add_subplot(111)
        self._plot()
        self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        plt.show()
        print('Click peak to select, right click to delete')

    def _plot(self):
        self.ax.plot(self.freqs, self.phasediffs)
        self.ax.set_xlabel('Frequency (Hz)')
        self.ax.set_ylabel('\Delta phase')
        for ind in self.peak_inds:
            self.ax.axvline(self.freqs[ind], linestyle='-.', color='orange')
        if self.orig_qubitdict is not None:
            for qubitid, freq in self.orig_qubitdict.items():
                freqind = np.argmin(np.abs(freq - self.freqs))
                self.ax.plot(self.freqs[freqind], self.phasediffs[freqind], 'o', color='g')
                self.ax.annotate(qubitid, (self.freqs[freqind], 
                    self.phasediffs[freqind]+0.1), color='g', ha='center')


    def _on_click(self, event):
        freq = event.xdata
        freqind = np.argmin(np.abs(freq - self.freqs))

        if event.button == 1:
            peakind = np.argmax(np.abs(self.phasediffs[freqind-self.snaprad:freqind+self.snaprad+1]))\
                    + freqind - self.snaprad
            if peakind in self.peak_inds:
                print('Already found peak at {} GHz'.format(self.freqs[peakind]/1.e9))
            
            else:
                self.peak_inds = np.append(self.peak_inds, peakind)
                print('Adding peak at {} GHz'.format(self.freqs[peakind]/1.e9))
                self.ax.clear()
                self._plot()
                self.fig.canvas.draw()


        elif event.button == 3:
            min_peak_dist = np.min(np.abs(freqind - self.peak_inds))
            if min_peak_dist <= self.snaprad:
                todelete_ind = np.argmin(np.abs(freqind - self.peak_inds))
                print('Removing peak at {} GHz'.format(self.freqs[self.peak_inds[todelete_ind]]/1.e9))
                self.peak_inds = np.delete(self.peak_inds, todelete_ind)
                self.ax.clear()
                self._plot()
                self.fig.canvas.draw()

        print(self.peak_inds)
                
            
            

    

def find_peaks_phasediff(phases, sig_thresh=2):
    """
    Parameters
    ----------
        phases : numpy array
            unwrapped, detrended phases (arg(S_11))
        sig_thresh : float
            n_sigma cutoff for peak finding algorithm
    """
    phasediffs = np.diff(phases)
    #fs = np.average(np.diff(freqs))
    #filt = signal.firwin(100, (lpf_coeff*fs/2, hpf_coeff*fs/2), pass_zero=False, window=('chebwin', 100), fs=fs)
    #filt_phasediff = np.convolve(phasediffs, filt, mode='same')

    peaks = signal.find_peaks(np.abs(phasediffs), height=sig_thresh*np.std(phasediffs))
    return peaks[0]

    

def run_vna(bw, n_samples, amplitude, t0, calirepo='submodules/qchip'):
    vna=c_vna(qubitid='vna',calirepo=calirepo)
    lor=vna.opts['wiremap'].lor #where do these come from? they should either not be class attributes or stay in VNA
    bw=vna.opts['chassis'].fsample
    #fx=numpy.linspace(6.2e9,6.7e9,2000)
    fx=numpy.linspace(lor-bw/2,lor+bw/2,2000)
    vna.seqs(fx,t0=692e-9,amp=0.5)
    vna.run(100)

    #qubitdict = 
    	
    fig1=pyplot.figure(figsize=(15,8))
    sub1=fig1.subplots(3,2).reshape(6)
    vna.iqapplot(sub1)
    fig2=pyplot.figure(figsize=(15,8))
    pyplot.plot(vna.fx[0:len(vna.p)],vna.p)
    pyplot.plot(vna.fx[vna.p2],vna.p[vna.p2],'x')
    
    fig3=pyplot.figure(figsize=(8,8))
    ax3=fig3.subplots()
    vna.iqplot(ax3,hexbin=False,marker='.')
    
    #pyplot.savefig('fig1.pdf')
    pyplot.show()

if __name__=='__main__':
    pass
