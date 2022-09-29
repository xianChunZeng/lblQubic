import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import ipdb
import qick
from socProxy import make_proxy

class VNAProgram(qick.AveragerProgram):
    """
    Adapted from SingleToneSpectroscopyProgram in qick_demos
    """
    def initialize(self):
        cfg=self.cfg   
        self.declare_gen(ch=cfg["dac_ch"], nqz=1) #Readout
        self.declare_readout(ch=cfg["adc_ch"], length=cfg["readout_length"],
                            freq=cfg["frequency"], gen_ch=cfg["dac_ch"])
        freq=self.freq2reg(cfg["frequency"], gen_ch=cfg["dac_ch"], ro_ch=0)  # convert frequency to dac frequency (ensuring it is an available adc frequency)
        self.set_pulse_registers(ch=cfg["dac_ch"], style="const", freq=freq, phase=0, gain=cfg["gain"],
                                length=cfg["readout_length"])

        self.synci(200)  # give processor some time to configure pulses
    
    def body(self):
        self.measure(pulse_ch=self.cfg["dac_ch"], 
             adcs=[self.cfg["adc_ch"]],
             adc_trig_offset=self.cfg["adc_trig_offset"],
             wait=True,
             syncdelay=self.us2cycles(self.cfg["relax_delay"]))

def run_vna(qchip, instrument_cfg, freqs, n_samples=1000, tof=180, amplitude=1):
    """
    freqs: list of frequencies in Hz
    tof: time of flight (or t_lo) in clock cycles (todo: change to seconds)
    """
    default_reps = 500
    default_ro_len = 600
    dac_ch = instrument_cfg['wiremap'].chanmapqubit['vna.rdrv']
    adc_ch = instrument_cfg['wiremap'].chanmapqubit['vna.read']
    cfg = {'reps': default_reps, 'soft_avgs': n_samples, 'dac_ch': dac_ch, 
            'adc_ch': adc_ch, 'readout_length': default_ro_len, 'relax_delay':500,
            'gain':amplitude*30000, 'adc_trig_offset': tof}
    soc, soccfg = make_proxy(instrument_cfg['ip'])
    soc.rfb_set_lo(instrument_cfg['wiremap'].lor/1.e6)
    #TODO: replace with wiremap params:
    soc.rfb_set_gen_rf(gen_ch=dac_ch, att1=6, att2=9)
    soc.rfb_set_ro_rf(ro_ch=adc_ch, att=30)
    #ipdb.set_trace()
    freqs = np.copy(freqs)

    freqs -= instrument_cfg['wiremap'].lor
    freqs /= 1.e6

    i = []
    q = []

    for f in tqdm(freqs):
        cfg['frequency'] = f
        rspec = VNAProgram(soccfg, cfg)
        avgi,avgq=rspec.acquire(soc, load_pulses=True)
        #ipdb.set_trace()
        i.append(avgi[0][0])
        q.append(avgq[0][0])

    i = np.asarray(i)
    q = np.asarray(q)
    phases = np.angle(i + 1j*q)
    phasefit = np.polyfit(freqs, np.unwrap(phases), deg=1)
    #plt.plot(freqs+instrument_cfg['wiremap'].lor/1.e6, np.unwrap(phases) - phasefit[0]*freqs - phasefit[1])
    #plt.show()
    return i + 1j*q

def run_punchout(qubitid, qchip, instrument_cfg, bandwidth=6.e6, n_freq=200, atten_start=5, atten_stop=25, \
        atten_step=1, amplitude=1, n_samples=1000, tof=180):
    #first run vna for global phase correction
    default_reps = 500
    default_ro_len = 600
    centerfreq = qchip.qubits[qubitid].readfreq
    vnafreqs = np.linspace(centerfreq-50e6, centerfreq+50e6, 500) #RF, in Hz
    s21 = run_vna(qchip, instrument_cfg, vnafreqs)
    phases = np.angle(s21)
    phasefit = np.polyfit(vnafreqs - instrument_cfg['wiremap'].lor, np.unwrap(phases), deg=1)

    dac_ch = instrument_cfg['wiremap'].chanmapqubit[qubitid + '.rdrv']
    adc_ch = instrument_cfg['wiremap'].chanmapqubit[qubitid + '.read']

    soc, soccfg = make_proxy(instrument_cfg['ip'])
    soc.rfb_set_lo(instrument_cfg['wiremap'].lor/1.e6)
    #TODO: replace with wiremap params:
    soc.rfb_set_gen_rf(gen_ch=dac_ch, att1=0, att2=atten_start)
    soc.rfb_set_ro_rf(ro_ch=adc_ch, att=30)
    #ipdb.set_trace()

    centerfreq_if = centerfreq - instrument_cfg['wiremap'].lor
    freqs = np.linspace(centerfreq_if - bandwidth/2, centerfreq_if + bandwidth/2, n_freq) #IF, in Hz
    attens = np.arange(atten_start, atten_stop, atten_step)

    cfg = {'reps': default_reps, 'soft_avgs': n_samples, 'dac_ch': dac_ch, 
            'adc_ch': adc_ch, 'readout_length': default_ro_len, 'relax_delay':500,
            'gain':amplitude*30000, 'adc_trig_offset': tof}

    s21 = np.zeros((len(attens), len(freqs)), dtype=complex)
    ipdb.set_trace()

    for i, atten in enumerate(tqdm(attens)):
        soc.rfb_set_gen_rf(gen_ch=dac_ch, att1=0, att2=atten)
        for j, freq in enumerate(tqdm(freqs)):
            cfg['frequency'] = freq/1.e6
            rspec = VNAProgram(soccfg, cfg)
            avgi,avgq = rspec.acquire(soc, load_pulses=True)
            #ipdb.set_trace()
            s21[i, j] = avgi[0][0] + 1j*avgq[0][0]

    phases = np.unwrap(np.angle(s21))
    for i in range(len(attens)):
        phases[i] -= (phasefit[0]*freqs + phasefit[1])

    return s21, phases, phasefit, attens, freqs


