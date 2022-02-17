import argparse
from matplotlib import pyplot as plt
from chipcalibration.punchout import run_punchout
import chipcalibration.vna as vna
import chipcalibration.punchout as po


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-bw', '--bandwidth', default=vna.VNA_BANDWIDTH,\
            help='VNA sweep BW in Hz, default: {}'.format(vna.VNA_BANDWIDTH))
    parser.add_argument('-a', '--amplitude', default=vna.AMPLITDUE,\
            help='tone amplitude, default: {}'.format(vna.AMPLITUDE))
    parser.add_argument('-n', '--n-freq', default=vna.N_FREQ_POINTS,\
            help='N points in VNA, default: {}'.format(vna.N_FREQ_POINTS))
    parser.add_argument('--n-samples', default=vna.N_SAMPLES,\
            help='number of samples in readout buffer, default: {}'.format(vna.N_SAMPLES))
    parser.add_argument('--update-cfg', action='store_true', help='overwrite config file')
    parser.add_argument('--punchout', action='store_true', help='run punchout after VNA on all identified peaks')
    parser.add_argument('--punchout-bandwidth', default=vna.FBW, 
            help='frequency sweep bandwidth in Hz, default {}'.format(vna.FBW))
    parser.add_argument('--n-freq-punchout', default=po.N_FREQ, 
            help='frequency sweep bandwidth in Hz, default {}'.format(po.N_FREQ))
    parser.add_argument('--atten-start', default=po.ATTEN_START, 
            help='starting (high) attenuation, default{}. Note that atten values are negative, so lower value means lower power'.format(po.ATTEN_START))
    parser.add_argument('--atten-stop', default=po.ATTEN_STOP, 
            help='ending (low) attenuation, default{}'.format(po.ATTEN_STOP))
    parser.add_argument('--atten-step', default=po.ATTEN_STEP, 
            help='dB increment to use in atten sweep, default{}'.format(po.ATTEN_STEP))
    args = parser.parse_args()

    peak_freqs = vna.run_vna(args.bandwidth, args.n_freq, args.n_samples, args.amplitude)

    if args.punchout:
        qubitids = ['peak {}'.format(i) for i in range(len(peak_freqs))]
        res_freqs, attens = po.run_punchout(qubitids, args.punchout_bandwidth, args.n_freq_punchout,\
                args.atten_start, args.atten_stop, args.atten_step, args.n_samples)
