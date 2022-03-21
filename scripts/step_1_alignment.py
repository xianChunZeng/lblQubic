"""
Script for determining delay time between readout DAC and ADC.
TODO:
    - what are the below params (N_BUF, TLO, etc), and which should 
      be hardcoded vs user specified
    - is anything to be written to the config file
"""

import matplotlib.pyplot as plt
from chipcalibration.alignment import run_alignment
from qubic.qubic.envset import load_chip
import argparse

N_BUF = 1000
TLO = 692e-9
MON_SEL0 = 0
MON_SEL1 = 3


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='run ADC/DAC alignment procedure')
    parser.add_argument('--qchip', help='name of chip (corresponds to qchip directory; e.g. X4Y2)')
    parser.add_argument('--cfg-file', default=None,
            help='path to qubit config file (can be absolute or relative to calirepo/qchip directory. if none just use default file for qchip')
    parser.add_argument('--calirepo-dir', default='../submodules/qchip', 
            help='path to gitrepo containing chip calibrations, default {}'.format('../submodules/qchip'))
    args = parser.parse_args()

    qchip, inst_cfg = load_chip(args.calirepo_dir, args.qchip, args.cfg_file)
    fig1=plt.figure(1,figsize=(15,8))
    run_alignment(qchip, inst_cfg, TLO, MON_SEL0, MON_SEL1, N_BUF, fig1)
