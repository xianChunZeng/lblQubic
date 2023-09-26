import matplotlib.pyplot as plt
import numpy as np
import qubic.toolchain as tc
from qubic.state_disc import GMMManager
from scipy.optimize import curve_fit
from collections import OrderedDict

ACC_BUFSIZE = 1000

class RabiAmpSweeper:
    """
    Parameters:
        input_register: the register of target qubits to run Rabi circuits on
        amp_bounds: upper and lower bounds for the amplidue intervals
        num_partitions: number of different amplitude partitions to use
        twidth_target: the target time of the calibrated Xpi/2 gate
        qchip: calibration config 
        fpga_config
        channel_configs

    Public functions:
        run(self, num_shots): runs the rabi circuits
        remake_parititons(self, amp_bounds, num_paritions): remakes the rabi circuits with new amplitude bound
        show_rabi_oscillations(self, target_qid, sub_register=None, show_fits=False): plot the rabi oscillation vs amplitudes with options to show idling qubits and the fits



    Define circuits, take data, and fit 

    """

    def __init__(self, register, twidth_target, qchip, fpga_config, channel_configs, amp_range=(0.0, 1.0), num_partitions=100,rabigate='rabi'):
        """
        Create rabi amplitude circuits according to input parameters, then compile to asm binaries.
        """
        if type(register) is not list:
            register = [register]
        self.register = register
        self.num_partitions = num_partitions
        self.twidth = twidth_target
        self.qchip = qchip
        self.fpga_config = fpga_config
        self.channel_configs = channel_configs
        self.readout_chanmap = {qid : str(channel_configs[qid + '.rdlo'].core_ind) for qid in register}
        self._assembled_circuits = dict()
        self.fits = None
        self.gmm_manager = GMMManager(chanmap_or_chan_cfgs=channel_configs)
        
        self._set_amplitude_partitions(amp_range[0], amp_range[1], num_partitions) 
        self.circuits = OrderedDict() # a batch of circuits for each target drive qubit
        self.rabigate=rabigate
        for qid in register:
            self.circuits[qid] = self._make_rabi_circuits(qid)

    def _make_rabi_circuits(self, drvqubit):
        """
        Make list of circuits used for rabi measurement. and the list of pulse width. So there will be a total of 
        1 circuits, each of which contains len(pulse_widths) measurements. A 400 us 
        delay is inserted between each measurement.
        """
        circuits = []
        for amp in self.amplitudes:
            cur_circ = []
            cur_circ.append({'name': 'delay', 't': 400.e-6})
            cur_circ.append({'name': self.rabigate, 'qubit': [drvqubit], 'modi': {(0, 'amp'): amp, (0, 'twidth') : self.twidth}})
            cur_circ.append({'name': 'barrier', 'qubit': self.register})
            for qid in self.register:
                    cur_circ.append({'name': 'read', 'qubit': [qid]})
            circuits.append(cur_circ)
        return circuits
    
    def _set_amplitude_partitions(self, lower_bound, upper_bound, num_partitions):
        """
        set the parition of amplitudes for corresponding rabi circuits
        """
        self.amplitudes = np.linspace(lower_bound, upper_bound, num_partitions)
                
    

    def run_and_fit(self, circuit_runner, num_samples, prior_fit_params, use_fft=True):
        """
        Run Rabi amplitude pulses on all the qubits and record the Xpi/2 gate times 

        Parameters
        ----------
            circuit_runner : CircuitRunner object
            nsamples : int
            prior_fit_params : dict
                keys are qubits, values are lists of parameters:
                    [A, B, drive_period, phi] such that
                    one_state_pop = A*cos(2*pi*x/drive_period + phi) + B
        """
        self._take_all_data(circuit_runner, num_samples)
        self._fit_gmm()
        self._fit_count_data(prior_fit_params, use_fft)

        #self._fit_gmm()

    def _take_all_data(self, job_manager, num_samples):
        """
        Take all the data needed to fit the rabi oscillations on the register

        raw_shots is a dictionary of lists of dictionaries
        dataset[drive_qubit][circuit_index][read_qubit]
        """
        self.num_shots = num_samples
        self.raw_shots = dict()
        for idx, drive_qid in enumerate(self.register): 
            print(f"Taking data for qubit {drive_qid} in batch {idx + 1} of {len(self.register)}")
            self.raw_shots[drive_qid] = job_manager.build_and_run_circuits(self.circuits[drive_qid], num_samples)['s11']

    def show_count_oscillations(self, target_qid, sub_register=None, show_fits=False):
        """
        plot average count data vs drive amplitude
        """
        if sub_register is None:
            sub_register = [target_qid]
        assert self.dataset is not None, "Must collect and classify the data"

        fig, axs = plt.subplots(len(sub_register), 1)
        avg_response = {qid : np.zeros(self.num_partitions) for qid in sub_register}
        
        if len(sub_register) == 1:
            for cidx in range(self.num_partitions):
                avg_response[target_qid][cidx] = np.average(self.dataset[target_qid][target_qid][cidx].flatten())
            axs.scatter(self.amplitudes, avg_response[sub_register[0]])
            axs.set_title(f'Counts on {target_qid} with drive on {target_qid}')
        else:
            for cidx in range(self.num_partitions):
                for qid in sub_register:
                    avg_response[qid][cidx] = np.average(self.dataset[target_qid][qid][cidx].flatten())
            for idx, qid in enumerate(sub_register):
                axs[idx].set_title(f'Counts on {qid} with drive on {target_qid}')
                axs[idx].scatter(self.amplitudes, avg_response[qid])
        if show_fits is True:
            if len(sub_register) == 1:
                axs.plot(self.amplitudes, self._cos(self.amplitudes, *self.fits[target_qid][0]), c='red')
            else:
                axs[sub_register.index(target_qid)].plot(self.amplitudes, self._cos(self.amplitudes, *self.fits[target_qid][0]), c='red')
        
        plt.tight_layout()
        plt.show()
        
    

    def show_raw_rabi_oscillations(self, target_qid, sub_register=None, show_fits=False):
        """
        plot real IQ response vs drive amplitude for given target drive and sub_register
        """
        if sub_register is None:
            sub_register = [target_qid]
        assert self.dataset is not None, "Must collect the data"

        fig, axs = plt.subplots(len(sub_register), 1)
        for idx, read_qid in enumerate(sub_register): 
            average_responses = [np.average(self.dataset[target_qid][read_qid][i].real) for i in range(len(self.amplitudes))]
            if len(sub_register) == 1:
                axs.scatter(self.amplitudes, average_responses, c='black')
                axs.set_title(f'{read_qid} response with drive on {target_qid}')
            else:
                axs[idx].scatter(self.amplitudes, average_responses, c='black')
                axs[idx].set_title(f'{read_qid} response with drive on {target_qid}')
        plt.tight_layout()

        if show_fits is True:
            if len(sub_register) == 1:          
                axs.plot(self.amplitudes, self._cos(self.amplitudes, *self.fits[target_qid][0]), c='red')
            else:
                axs[sub_register.index(target_qid)].plot(self.amplitudes, self._cos_function(self.amplitudes, *self.fits[target_qid][0]), c='red')

        plt.show()

    def _cos(self, x, A, B, drive_period, phi):
        return A*np.cos(2*np.pi*x/drive_period - phi) + B

    def _fit_raw_data(self, prior_fit_params):
        """
        fit the real IQ response to a cosine
        """ 
#         self.fits = dict()
#         for qid in self.register:
#             average_response = np.array([np.average(self.dataset[qid][qid][i]) for i in range(self.num_partitions)])
#             self.fits[qid] = curve_fit(self._cos_function, self.amplitudes, average_response, prior_fit_params[qid]) 
    
    def _fit_count_data(self, prior_fit_params, use_fft):
        """
        fit the count data to a cosine
        """
        self.fits = dict()
        for qid in self.register:
            average_response = np.array([np.average(self.dataset[qid][qid][i]) for i in range(self.num_partitions)])
            try:
                if use_fft:
                    # this is "frequency" in terms of the rabi amplitude oscillation period
                    freq_ind_max = np.argmax(np.abs(np.fft.rfft(average_response)[1:])) + 1
                    freq_max = np.fft.rfftfreq(len(average_response), np.diff(self.amplitudes)[0])[freq_ind_max]
                    prior_fit_params[qid][2] = 1/freq_max
                    self.fits[qid] = curve_fit(self._cos, self.amplitudes, average_response, prior_fit_params[qid]) 
            except:
                print(f'Could not fit {qid}')

    def _fit_gmm(self):
        """
        fit the GMM manager to all the collected data and sort based 
            on response of the amplitude 0 circuit
        """
        # [0] fit the GMM's on the total raw data stream out of each readout line
        all_raw_data = {cid : np.array([]) for cid in self.readout_chanmap.values()}
        for qid in self.register:
            for chan_id in self.raw_shots[qid].keys():
                for cidx in range(self.num_partitions):
                    all_raw_data[chan_id] = np.hstack((all_raw_data[chan_id], self.raw_shots[qid][chan_id][cidx][:, 0]))
        self.all_raw_data = all_raw_data
        self.gmm_manager.fit(all_raw_data)
        # [1] relabel the GMM's based on the first (0 amplitude) circuit
        for qid in self.register:
            self.gmm_manager.gmm_dict[qid].set_labels_maxtomin(self.raw_shots[qid][self.readout_chanmap[qid]][0], [0, 1])
           
        # [2] convert the raw shot data into a dataset using the trained GMM
        self.dataset = {qid : [] for qid in self.register}
        for qid in self.register:
            self.dataset[qid] = self.gmm_manager.predict(self.raw_shots[qid])
        
        # dataset structure: 
#         # {drive qubit}{readout qubit}[circuit index, shot index]
#         self.dataset = { drive_qid:  { ro_qid : np.zeros((self.num_partitions, self.num_shots)) for ro_qid in self.register} for drive_qid in self.register}
#         for drive_qid in self.register:
#             for ro_qid in self.register:
#                 self.dataset[drive_qid][ro_qid] = self.gmm_man
        
                        
        
        
#         self.gmm_manager.set_labels_maxtomin({chan: shots.flatten() 
#                                               for chan, shots in self.raw_shots['Q0'][0].items()}, [0, 1])
#         self.state_disc_shots = self.gmm_manager.predict(self.raw_shots)
#         self.ones_frac = {qubit: np.sum(self.state_disc_shots[qubit],axis=1) for qubit in self.state_disc_shots.keys()}
#         self.zeros_frac = {qubit: np.sum(self.state_disc_shots[qubit] == 0,axis=1) for qubit in self.state_disc_shots.keys()}
    def update_qchip(self, qchip,x90gate='X90'):
        cal_drive_amps = {qid: self.fits[qid][0][2]/4 for qid in self.register} #1/4 the drive period for each qubit
        for qid in self.register:
            qchip.gates[qid + x90gate].contents[0].amp = cal_drive_amps[qid]
            qchip.gates[qid + x90gate].contents[0].twidth = self.twidth
