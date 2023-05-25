import networkx as nx
import logging
import matplotlib.pyplot as plt

class CalibrationGraph:

    def __init__(self):
        self.graph = nx.DiGraph()

    def add_calibration_step(self, name, cal_obj, qubits, predecessor_nodes=None, shots_per_circuit=1000):
        self.graph.add_node(name, cal_obj=cal_obj, qubits=qubits, success=False, results=None, shots_per_circuit=shots_per_circuit)
        if predecessor_nodes is not None:
            for node in predecessor_nodes:
                assert node in self.graph.nodes()
                self.graph.add_edge(node, name)

    def run_calibration(self, job_manager, initial_qchip, initial_gmm_manager=None, show_plots=False, save_plots=False):
        self.qchip = initial_qchip
        if initial_gmm_manager is not None:
            self.gmm_manager = initial_gmm_manager
        else:
            self.gmm_manager = job_manager.gmm_manager

        node_list = nx.topological_sort(self.graph)
        for node in node_list:
            exec_node = True
            for pred_node in self.graph.predecessors(node):
                if self.graph.nodes[pred_node]['success'] == False:
                    logging.getLogger(__name__).warning('{} not successful; skipping {}'.format(pred_node, node))
                    exec_node = False
                    break
            if exec_node == False:
                break

            job_manager.gmm_manager = self.gmm_manager
            
            cal_obj = self.graph.nodes()[node]['cal_obj']
            nshots = self.graph.nodes()[node]['shots_per_circuit']

            try:
                cal_obj.run_and_report(job_manager, nshots, self.qchip)

            except Exception as e:
                logging.getLogger(__name__).warning('{} not successful; error: {}'.format(node, e))
                continue

            if show_plots:
                figure = plt.figure()
                cal_obj.plot_results(figure)
                figure
                plt.show()

            cal_obj.update_qchip(self.qchip)
            cal_obj.update_gmm_manager(self.gmm_manager)
            self.graph.nodes()[node]['success'] = True
            self.graph.nodes()[node]['results'] = cal_obj.results

    def append_by_qubit(self, cal_graph):
        sum_graph = nx.union(self.graph, cal_graph.graph, rename=('0_', '1_'))
        end_nodes = ['0_' + node for node in self.graph.nodes if self.graph.succ[node] == {}]
        start_nodes = ['1_' + node for node in cal_graph.graph.nodes if cal_graph.graph.pred[node] == {}]

        for enode in end_nodes:
            for snode in start_nodes:
                if any(enode_qubit in sum_graph.nodes[snode]['qubits'] for enode_qubit in sum_graph.nodes[enode]['qubits']):
                    sum_graph.add_edge(enode, snode)

        self.graph = sum_graph




