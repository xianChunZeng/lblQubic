import matplotlib.pyplot as plt
import numpy as np
import sys

import pygsti.data

import qubic.toolchain as tc
#import qubic.run as rc
from qubic.state_disc import GMMManager
from scipy.optimize import curve_fit
from pygsti.data import dataset
from QubicPig.chipspec import ChipSpec


def compile_circuits_to_asm(circuit_batch, qchip, fpga_config, channel_configs):
    """
    takes a dictionary with pygsti circuit keys and qubic instruction values
    outputs a dictionary with pygsti circuit keys and asm binary values
    :return:
    """
    circ_keys = list(circuit_batch.keys())
    instruction_list = [circuit_batch[c] for c in circ_keys]
    compiled_circs = tc.run_compile_stage(instruction_list,
                                          fpga_config, qchip)
    assembled_circs = tc.run_assemble_stage(compiled_circs, channel_configs)
    return {circ_keys[i]: assembled_circs[i] for i in range(len(circ_keys))}

class JobManager:
    """

    """

    def __init__(self, chipspec, fpga_config, channel_configs):
        """
        Create rabi amplitude circuits according to input parameters, then compile to asm binaries.
        """
        self.chipspec = chipspec
        self.fpga_config = fpga_config
        self.channel_configs = channel_configs
        self.readout_chanmap = {qid: channel_configs[qid + '.rdlo'].core_ind for qid in chipspec.spec.qubit_labels}


    def run_circuits(self, asm_batch, circuit_runner, num_shots):
        circ_labels = list(asm_batch.keys())
        circ_vals = [asm_batch[c] for c in circ_labels]
        return circuit_runner.run_circuit_batch(circ_vals, num_shots)

    def collect_pygsti_dataset(self, pcircs_and_qinstructions, num_shots, circuit_runner, new_chipspec=None):
        """
        assume we only read each circuit once

        :param pcircs_and_qinstructions:
        :param num_shots:
        :param circuit_runner:
        :param new_chipspec:
        :return:
        """
        # replace the chipspec if given
        if new_chipspec is not None:
            self.chipspec = new_chipspec
        # compile circuits
        asm_circs = compile_circuits_to_asm(pcircs_and_qinstructions, self.chipspec.qchip,
                                            self.fpga_config, self.channel_configs)
        # collect shots and classify
        raw_iq_shots = self.run_circuits(asm_circs, circuit_runner, num_shots)
        unsorted_counts = self.chipspec.gmm_manager.predict(raw_iq_shots)
        # sort shots into the dataset -- probably could be made more efficient
        dataset = pygsti.data.DataSet()
        for circ_idx, circ in enumerate(pcircs_and_qinstructions.keys()):
            count_dict = dict()
            for shot_idx in range(num_shots):
                shot_string = ''
                for qid in self.chipspec.spec.qubit_labels:
                    shot_string = shot_string.join(str(unsorted_counts[qid][circ_idx][shot_idx][0]))
                if shot_string in count_dict:
                    count_dict[shot_string] += 1
                else:
                    count_dict[shot_string] = 1
            dataset.add_count_dict(circ, count_dict)
        return dataset






