"""
Script for determining delay time between readout DAC and ADC.
TODO:
    - what are the below params (N_BUF, TLO, etc), and which should 
      be hardcoded vs user specified
    - is anything to be written to the config file
"""

import matplotlib.pyplot as plt
from qubic.qcvv.alignment import c_alignment
import argparse

N_BUF = 1000
TLO = 692e-9
MON_SEL0 = 0
MON_SEL1 = 3

def run_alignment(tlo, mon_sel0, mon_sel1, nbuf, fig=None):
    """
    Parameters
    ----------
        tlo : float
            ??
        mon_sel0 : int
            ??
        mon_sel1 : int
            ??
        fig : matplotlib Figure
            (optional) figure for generating a plot
    """
    alignment=c_alignment(qubitid='alignment',calirepo='submodules/qchip',debug=False,sim=True)
    alignment.seqs(tlo=tlo,mon_sel0=0,mon_sel1=3)
    alignment.run(nbuf)
    if fig is not None:
        alignment.plot(fig)
        plt.show()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='run ADC/DAC alignment procedure')
    #parser.add_argument...
    fig1=plt.figure(1,figsize=(15,8))
    run_alignment(TLO, MON_SEL0, MON_SEL1, N_BUF, fig1)
