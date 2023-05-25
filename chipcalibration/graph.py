import networkx as nx
import logging
import matplotlib.pyplot as plt

class CalibrationGraph:

    def __init__(self, job_manager, qchip, gmm_manager=None):
        self.graph = nx.DiGraph()
        self.qchip = qchip
        self.job_manager = job_manager
        self.gmm_manager = gmm_manager

    def add_calibration_step(self, name, cal_obj, predecessor_nodes=None, shots_per_circuit=1000):
        self.graph.add_node(name, cal_obj=cal_obj, success=False, results=None, shots_per_circuit=shots_per_circuit)
        if predecessor_nodes is not None:
            for node in predecessor_nodes:
                assert node in self.graph.nodes()
                self.graph.add_edge(node, name)

    def run_calibration(self, show_plots=False, save_plots=False):
        node_list = nx.topological_sort(self.graph)
        for node in node_list:
            exec_node = True
            for pred_node in self.graph.predecessors(node):
                if pred_node['success'] == False:
                    logging.getLogger(__name__).warning('{} not successful; skipping {}'.format(pred_node, node))
                    exec_node = False
                    break
            if exec_node == False:
                break

            if self.gmm_manager is not None:
                self.job_manager.gmm_manager = self.gmm_manager
            
            cal_obj = self.graph.nodes()[node]['cal_obj']
            nshots = self.graph.nodes()[node]['shots_per_circuit']

            try:
                cal_obj.run_and_report(self.job_manager, nshots, self.qchip)

            except Exception as e:
                logging.getLogger(__name__).warning('{} not successful; error: {}'.format(node, e))
                continue

            if show_plots:
                figure = plt.figure()
                cal_obj.plot_results(figure)
                plt.show()

            cal_obj.update_qchip(self.qchip)
            cal_obj.update_gmm_manager(self.gmm_manager)
            self.graph.nodes()[node]['success'] = True
            self.graph.nodes()[node]['results'] = cal_obj.results_dict



