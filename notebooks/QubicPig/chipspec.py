import pygsti.processors
from QubicPig import qupig
from pygsti.processors.processorspec import QubitProcessorSpec
from qubitconfig.qchip import QChip
from qubic.state_disc import GMMManager

class ChipSpec():
    def __init__(self, pspec: QubitProcessorSpec, qchip: QChip, gmm_manager: GMMManager):
        self.spec = pspec
        self.qchip = qchip
        self.gmm_manager = gmm_manager
        self.target_model = pygsti.models.modelconstruction.create_explicit_model(pspec)

    def update(self, keys_or_dict, value=None):
        self.qchip.update(keys_or_dict, value)
    def qubic_instructions(self, pygsti_circuit):
        qubic_circuit = list()
        qubic_circuit.append({'name': 'delay', 't': 400.e-6})
        for layer in pygsti_circuit:
            qubic_circuit.extend(qupig.parse_layer(layer))
            qubic_circuit.append({'name': 'barrier', 'qubit': list(self.spec.qubit_labels)})
        for qid in self.spec.qubit_labels:
            qubic_circuit.append({'name': 'read', 'qubit': [qid]}, )
        return qubic_circuit


