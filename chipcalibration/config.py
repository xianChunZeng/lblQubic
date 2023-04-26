import qubitconfig.qchip as qc
from qubic.rfsoc.hwconfig import load_channel_configs
import os

def load_configs(chipname, qchip_file='qubitcfg.json', channel_cfg_file='channel_config.json', qchip_dir='default'):
    #todo: remove hardcoded fpga config
    fpga_config = FPGAConfig(**{'fpga_clk_period': 2.e-9, 'alu_instr_clks': 5, 'jump_cond_clks': 5, 
                                'jump_fproc_clks': 5, 'pulse_regwrite_clks': 3})
    if qchip_dir=='default':
        qchip_dir = os.path.join(os.path.dirname(__file__), '../submodules/qchip')

    qchip_dir = os.path.join(qchip_dir, chipname)

    qchip = qc.QChip(os.path.join(qchip_dir, qchip_file))
    channel_configs = load_channel_configs(os.path.join(qchip_dir, channel_cfg_file))

    return fpga_config, qchip, channel_configs

