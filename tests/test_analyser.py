import unittest
import networkx as nx
from capsidgraph.analyser.analyse import (
    get_hole_size,
    _bisection_stop_condition,
    get_framentation_probability,
    bisection,
)
from capsidgraph.analyser.fragment import (
    probability_fragment,
    energy_edges_fragment,
    energy_nodes_fragment,
)
from capsidgraph.analyser.util import _init_nodes_energy
from capsidgraph.analyser import (
    get_fragmentation_energy_threshold_edge,
    get_fragmentation_energy_threshold_node,
    get_fragmentation_probability_threshold_edge,
    get_fragmentation_probability_threshold_node,
    get_fragmentation_probability_random_edge_removal,
    get_fragmentation_probability_random_node_removal,
    get_fragment_size_distribution,
    get_hole_size_distribution,
)


class TestAnalyser(unittest.TestCase):
    def test_fragment_probability(self):
        G = nx.read_adjlist("tests/testcase1.adjlist")
        G_ = probability_fragment(
            G, {"fragmentation": 0.5, "fragmentation_type": "nodes"}
        )
        self.assertLessEqual(len(G_.nodes), len(G.nodes))

        G_ = probability_fragment(
            G, {"fragmentation": 1, "fragmentation_type": "edges"}
        )
        self.assertTrue(nx.is_empty(G_))
        self.assertEqual(len(G_.nodes), 72)
        G_ = probability_fragment(
            G, {"fragmentation": 1, "fragmentation_type": "nodes"}
        )
        self.assertEqual(len(G_.nodes), 0)

    def test_bisection_stop_condition(self):
        self.assertTrue(
            _bisection_stop_condition(100000, 0.9, {"error_probability": 0.01})
        )
        self.assertTrue(
            _bisection_stop_condition(100000000, 0.5, {"error_probability": 0.01})
        )
        self.assertFalse(
            _bisection_stop_condition(10000, 0.5, {"error_probability": 0.01})
        )

    def test_get_fragmentation_probability(self):
        G = nx.read_adjlist("tests/testcase1.adjlist")
        pfrag, n = get_framentation_probability(
            G,
            1000,
            probability_fragment,
            fragment_settings={
                "fragmentation": 0.5,
                "fragmentation_type": "nodes",
            },
        )
        self.assertLess(pfrag, 1)
        self.assertGreater(pfrag, 0)
        self.assertEqual(n, 1000)

        pfrag, n = get_framentation_probability(
            G,
            _bisection_stop_condition,
            probability_fragment,
            stop_condition_settings={"error_probability": 0.01, "min_iterations": 100},
            fragment_settings={
                "fragmentation": 0.5,
                "fragmentation_type": "nodes",
            },
        )
        self.assertLess(pfrag, 1)
        self.assertGreater(pfrag, 0)
        self.assertGreaterEqual(n, 100)

    def test_fragment_energy_nodes(self):
        G = nx.read_adjlist("tests/testcase1.adjlist")
        nx.set_edge_attributes(G, 1 / len(G.edges), "energy")
        _init_nodes_energy(G)
        for n in G.nodes:
            s = 0
            for nei in nx.neighbors(G, n):
                s += G.edges[(n, nei)]["energy"]
            self.assertEqual(G.nodes[n]["energy"], s)

    def test_fragment_energy_bonds(self):
        G = nx.read_adjlist("tests/testcase1.adjlist")
        nx.set_edge_attributes(G, 1 / len(G.edges), "energy")
        energy = 0.6
        G_ = energy_edges_fragment(G, {"fragmentation": energy})
        remaining_energy = 0
        for e in G_.edges:
            remaining_energy += G_.edges[e]["energy"]
        self.assertLessEqual(remaining_energy, 1 - energy)

    def test_fragment_energy_nodes(self):
        G = nx.read_adjlist("tests/testcase1.adjlist")
        nx.set_edge_attributes(G, 1 / len(G.edges), "energy")
        _init_nodes_energy(G)
        fragmentation_energy = 0.5
        G_ = energy_nodes_fragment(G, {"fragmentation": fragmentation_energy})
        remaining_energy = 0
        for e in G_.edges:
            remaining_energy += G_.edges[e]["energy"]
        self.assertAlmostEqual(remaining_energy, 1 - fragmentation_energy, places=2)

        for n in G.nodes:
            s = 0
            for nei in nx.neighbors(G, n):
                s += G.edges[(n, nei)]["energy"]
            self.assertEqual(G.nodes[n]["energy"], s)

    def test_bisection(self):
        steps = 3
        error_probability = 0.05
        G = nx.read_edgelist("tests/testcase2.edgelist")
        for e in G.edges:
            G.edges[e]["energy"] = 1 / len(G.edges)
        _init_nodes_energy(G)
        pf, n = bisection(
            G, steps, error_probability, energy_edges_fragment, fragment_settings={}, debug=True, debug_interval=100
        )
        self.assertAlmostEqual(pf, 0.375)

    def test_wrappers(self):
        G = nx.read_edgelist("tests/testcase2.edgelist")
        for e in G.edges:
            G.edges[e]["energy"] = 1 / len(G.edges)
        pf, n = get_fragmentation_energy_threshold_edge(G, 0.1, 3)
        self.assertAlmostEqual(pf, 0.375)
        pf, n = get_fragmentation_energy_threshold_node(G, 0.1, 3)
        self.assertAlmostEqual(pf, 0.875)
        pf, n = get_fragmentation_probability_threshold_edge(G, 0.1, 3)
        self.assertAlmostEqual(pf, 0.375)
        pf, n = get_fragmentation_probability_threshold_node(G, 0.1, 3)
        self.assertAlmostEqual(pf, 0.375)

        G = nx.read_adjlist("tests/AaLS_24.adjlist")
        p = get_fragmentation_probability_random_edge_removal(G, 0.4, 10000)
        self.assertAlmostEqual(p, 0.564, places=1)
        p = get_fragmentation_probability_random_node_removal(G, 0.4, 10000)
        self.assertAlmostEqual(p, 0.614, places=1)

    def test_fragment_size_distribution(self):
        iterations = 100
        G = nx.read_adjlist("tests/AaLS_24.adjlist")
        res = get_fragment_size_distribution(
            G,
            iterations,
            probability_fragment,
            {"fragmentation": 0.4, "fragmentation_type": "nodes"},
        )
        self.assertLessEqual(len(res), len(G.nodes))
        for i in res:
            self.assertGreaterEqual(i, 0)
            self.assertLessEqual(i, len(G.nodes))

    def test_hole_size(self):
        G = nx.from_edgelist(
            [
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (1, 5),
                (3, 6),
                (4, 7),
                (2, 8),
                (5, 8),
                (8, 7),
                (7, 6),
                (6, 5),
            ]
        )
        G_frag = G.copy()
        G_frag.remove_nodes_from([1, 2, 3, 4])
        self.assertEqual(get_hole_size(G_frag, G), 5)
        self.assertEqual(get_hole_size(G, G), 0)
        self.assertEqual(get_hole_size(nx.empty_graph(0), G), len(G.nodes))

    def test_hole_distribution(self):
        iterations = 100
        G = nx.read_adjlist("tests/AaLS_24.adjlist")
        res = get_hole_size_distribution(
            G,
            iterations,
            probability_fragment,
            {"fragmentation": 0.4, "fragmentation_type": "nodes"},
        )
        self.assertEquals(len(res), len(G.nodes) + 1)
        for i in res:
            self.assertGreaterEqual(i, 0)


if __name__ == "__main__":
    unittest.main()
