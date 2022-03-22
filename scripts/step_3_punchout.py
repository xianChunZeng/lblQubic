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
import chipcalibration.punchout as po


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Run punchout freq x atten sweep')
    parser.add_argument('qubitids', nargs='+', help='list of qubit identifiers; e.g. Q0, Q1, etc.')
    parser.add_argument('--qchip', help='name of chip (corresponds to qchip directory; e.g. X4Y2)')
    parser.add_argument('--save', type=str, nargs='?', const='default', 
            help='OVERWRITE qchip if specified, can optionally provide another file to overwrite instead')
    parser.add_argument('--bandwidth', default=po.FBW, 
            help='frequency sweep bandwidth in Hz, default {}'.format(po.FBW))
    parser.add_argument('--n-freq', default=po.N_FREQ, help='frequency sweep bandwidth in Hz, default {}'.format(po.N_FREQ))
    parser.add_argument('--atten-start', default=po.ATTEN_START, 
            help='starting (high) attenuation, default{}. Note that atten values are negative, so lower value means lower power'.format(po.ATTEN_START))
    parser.add_argument('--atten-stop', default=po.ATTEN_STOP, 
            help='ending (low) attenuation, default{}'.format(po.ATTEN_STOP))
    parser.add_argument('--atten-step', default=po.ATTEN_STEP, 
            help='dB increment to use in atten sweep, default{}'.format(po.ATTEN_STEP))
    parser.add_argument('--n-samples', default=po.N_SAMPLES, 
            help='number of samples in readout acquisition buffer, default{}'.format(po.N_SAMPLES))
    parser.add_argument('--cfg-file', default=None,
            help='path to qubit config file (can be absolute or relative to calirepo/qchip directory. if none just use default file for qchip')
    parser.add_argument('--calirepo-dir', default='../submodules/qchip', 
            help='path to gitrepo containing chip calibrations, default {}'.format('submodules/qchip'))
    args = parser.parse_args()

    qchip, inst_cfg, cfg_file = load_chip(args.calirepo_dir, args.qchip, args.cfg_file)

    qubit_dict = po.get_qubit_dict(args.qubitids, qchip)

    freqs, attens, qubitids = po.run_punchout(qubit_dict, qchip, inst_cfg, args.bandwidth, args.n_freq, args.atten_start, args.atten_stop, args.atten_step, args.n_samples)

    if args.save is not None:
        globalatten = inst_cfg['wiremap'].ttydev['rdrvvat']['default']
        for i, qubitid in enumerate(qubitids):
            qchip.updatecfg({('Qubits', qubitid, 'readfreq'): freqs[i]})
            amp = 10**((attens[i] + globalatten)/20)
            qchip.updatecfg({('Gates', qubitid + 'read', 'amp'): amp})
        

        if args.save=='default':
            qchip.save(cfg_file)

        else:
            if os.path.isabs(args.save):
                qchip.save(args.save)
            else:
                qchip.save(os.path.join(args.calirepo_dir, args.qchip, args.save))
