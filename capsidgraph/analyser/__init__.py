import networkx as nx
from .analyse import _get_framentation_probability, _bisection, get_fragment_size_distribution, get_hole_size_distribution
from .fragment import (
    _probability_fragment,
    _energy_edges_fragment,
    _energy_nodes_fragment,
)
from .util import _init_nodes_energy
from typing import Tuple


def get_fragmentation_probability_random_node_removal(
    G: nx.Graph, removal_probability: float, iterations: int
) -> float:
    pfrag, n = _get_framentation_probability(
        G,
        iterations,
        _probability_fragment,
        fragment_settings={
            "fragmentation": removal_probability,
            "fragmentation_type": "nodes",
        },
    )
    return pfrag


def get_fragmentation_probability_random_edge_removal(
    G: nx.Graph, removal_probability: float, iterations: int
) -> float:
    pfrag, n = _get_framentation_probability(
        G,
        iterations,
        _probability_fragment,
        fragment_settings={
            "fragmentation": removal_probability,
            "fragmentation_type": "edges",
        },
    )
    return pfrag

def get_fragmentation_probability_energy_node_removal(
    G: nx.Graph, removal_energy: float, iterations: int
) -> float:
    G_ = G.copy()
    _init_nodes_energy(G_)
    pfrag, n = _get_framentation_probability(
        G_,
        iterations,
        _energy_nodes_fragment,
        fragment_settings={
            "fragmentation": removal_energy,
        },
    )
    return pfrag


def get_fragmentation_probability_energy_edge_removal(
    G: nx.Graph, removal_energy: float, iterations: int
) -> float:
    pfrag, n = _get_framentation_probability(
        G,
        iterations,
        _energy_edges_fragment,
        fragment_settings={
            "fragmentation": removal_energy,
        },
    )
    return pfrag


# TODO: add tests for the wrappers
def get_fragmentation_probability_threshold_node(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    debug: bool = False,
) -> Tuple[float, int]:
    pf, n = _bisection(
        G,
        steps,
        error_probability,
        _probability_fragment,
        fragment_settings={"fragmentation_type": "nodes"},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
    )
    return pf, n


def get_fragmentation_probability_threshold_edge(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    debug: bool = False,
) -> Tuple[float, int]:
    pf, n = _bisection(
        G,
        steps,
        error_probability,
        _probability_fragment,
        fragment_settings={"fragmentation_type": "edges"},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
    )
    return pf, n


def get_fragmentation_energy_threshold_edge(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    debug: bool = False,
) -> Tuple[float, int]:
    pf, n = _bisection(
        G,
        steps,
        error_probability,
        _energy_edges_fragment,
        fragment_settings={},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
    )
    return pf, n


def get_fragmentation_energy_threshold_node(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    debug: bool = False,
) -> Tuple[float, int]:
    G_ = G.copy()
    _init_nodes_energy(G_)
    pf, n = _bisection(
        G_,
        steps,
        error_probability,
        _energy_nodes_fragment,
        fragment_settings={},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
    )
    return pf, n
