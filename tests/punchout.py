import step_3_punchout as punchout
from qubic.qubic.envset import load_chip


def getqubitdict0():
    qubitids = ['Q0', 'Q1', 'Q5']
    qchip, _ = load_chip('../submodules/qchip', 'X4Y2')
    qubitdict = punchout.get_qubit_dict(qubitids, qchip)
    for qubit in qubitids:
        assert(qubitdict[qubit] == qchip.getfreq(qubitid+'.readfreq'))
