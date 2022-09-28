import argparse
import os
from matplotlib import pyplot as plt
import numpy as np

from qubic.qubic.envset import load_chip
from chipcalibration.punchout import run_punchout
import chipcalibration.vna as vna
import chipcalibration.punchout as po


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--qchip', required=True, help='name of chip (corresponds to qchip directory; e.g. X4Y2)')
    parser.add_argument('-bw', '--bandwidth', default=vna.VNA_BANDWIDTH,\
            help='VNA sweep BW in Hz, default: {}'.format(vna.VNA_BANDWIDTH))
    parser.add_argument('-fc', '--center-freq', default=None,\
            help='VNA center freq in Hz, default: readout LO freq (from wiremap)')
    parser.add_argument('-a', '--amplitude', default=vna.AMPLITUDE,\
            help='tone amplitude, default: {}'.format(vna.AMPLITUDE))
    parser.add_argument('-n', '--n-freq', default=vna.N_FREQ_POINTS,\
            help='N points in VNA, default: {}'.format(vna.N_FREQ_POINTS))
    parser.add_argument('--n-samples', default=vna.N_SAMPLES,\
            help='number of samples in readout buffer, default: {}'.format(vna.N_SAMPLES))
    parser.add_argument('--update-cfg', action='store_true', help='overwrite config file')
    parser.add_argument('--punchout', action='store_true', help='run punchout after VNA on all identified peaks')
    parser.add_argument('--punchout-bandwidth', default=po.FBW, 
            help='frequency sweep bandwidth in Hz, default {}'.format(po.FBW))
    parser.add_argument('--n-freq-punchout', default=po.N_FREQ, 
            help='frequency sweep bandwidth in Hz, default {}'.format(po.N_FREQ))
    parser.add_argument('--atten-start', default=po.ATTEN_START, 
            help='starting (high) attenuation, default{}. Note that atten values are negative, so lower value means lower power'.format(po.ATTEN_START))
    parser.add_argument('--atten-stop', default=po.ATTEN_STOP, 
            help='ending (low) attenuation, default{}'.format(po.ATTEN_STOP))
    parser.add_argument('--atten-step', default=po.ATTEN_STEP, 
            help='dB increment to use in atten sweep, default{}'.format(po.ATTEN_STEP))
    parser.add_argument('--cfg-file', default=None,
            help='path to qubit config file (can be absolute or relative to calirepo/qchip directory. if none just use default file for qchip')
    parser.add_argument('--calirepo-dir', default='../submodules/qchip', 
            help='path to gitrepo containing chip calibrations, default {}'.format('submodules/qchip'))
    parser.add_argument('--save', type=str, nargs='?', const='default', 
            help='OVERWRITE qchip if specified, can optionally provide another file to overwrite instead')
    args = parser.parse_args()

    qchip, inst_cfg, cfg_file = load_chip(args.calirepo_dir, args.qchip, args.cfg_file)

    peak_freqs = vna.run_vna(qchip, inst_cfg, args.center_freq, args.bandwidth, args.n_freq, args.n_samples, args.amplitude)

    if args.punchout:
        qubitdict = {'peak {}'.format(i) : peak_freqs[i] for i in range(len(peak_freqs))}
        res_freqs, attens, _ = po.run_punchout(qubitdict, qchip, inst_cfg, args.punchout_bandwidth, args.n_freq_punchout,\
                args.atten_start, args.atten_stop, args.atten_step, args.n_samples)
        res_freqs = np.asarray(res_freqs)
        attens = np.asarray(attens)
        good_mask = res_freqs != None

        qubitids = ['Q{}'.format(i) for i in range(np.sum(good_mask))]
        res_freqs = res_freqs[good_mask]
        attens = attens[good_mask]

        po.update_qchip(qchip, inst_cfg, res_freqs, attens, qubitids)

    #else if args.save is not None:
    #    qubitids = ['Q{}'.format(i) for i in range(len(peak_freqs))]
    else:
        qubitids = ['Q{}'.format(i) for i in range(len(peak_freqs))]
        vna.update_qchip(qchip, peak_freqs, qubitids)

    if args.save is not None:
        if args.save=='default':
            qchip.save(cfg_file)

        else:
            if os.path.isabs(args.save):
                qchip.save(args.save)
            else:
                qchip.save(os.path.join(args.calirepo_dir, args.qchip, args.save))

