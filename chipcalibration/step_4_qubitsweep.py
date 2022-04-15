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
import os
import argparse
from qubic.qubic.envset import load_chip
import chipcalibration.qubitsweep as qs


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Run qubitsweep (vnaqubit)')
    parser.add_argument('qubitids', nargs='+', help='list of qubit identifiers; e.g. Q0, Q1, etc.')
    parser.add_argument('--qchip', help='name of chip (corresponds to qchip directory; e.g. X4Y2)', required=True)
    parser.add_argument('--save', type=str, nargs='?', const='default', 
            help='OVERWRITE qchip if specified, can optionally provide another file to overwrite instead')
    parser.add_argument('--n-samples', default=po.N_SAMPLES, 
            help='number of samples in readout acquisition buffer, default{}'.format(po.N_SAMPLES))
    parser.add_argument('--cfg-file', default=None,
            help='path to qubit config file (can be absolute or relative to calirepo/qchip directory. if none just use default file for qchip')
    parser.add_argument('--calirepo-dir', default='../submodules/qchip', 
            help='path to gitrepo containing chip calibrations, default {}'.format('submodules/qchip'))
    args = parser.parse_args()

    qchip, inst_cfg, cfg_file = load_chip(args.calirepo_dir, args.qchip, args.cfg_file)

    qubit_dict = po.get_qubit_dict(args.qubitids, qchip)

    sep = qs.run_qubit_sweep(['Q3'], ['Q2', 'Q3', 'Q4'], qchip, inst_cfg, np.linspace(5.e9, 5.5e9, 750))

