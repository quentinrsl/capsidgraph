import networkx as nx
from .analyse import (
    get_fragmentation_probability,
    bisection,
    get_fragment_size_distribution,
    get_hole_size_distribution,
    get_hole_size
)
from .fragment import (
    probability_fragment,
    strength_edges_fragment,
    strength_nodes_fragment,
)
from .util import _init_nodes_strength as init_nodes_strength
from typing import Tuple


def get_fragmentation_probability_random_node_removal(
    G: nx.Graph, removal_probability: float, iterations: int, process_number: int = 1 ,debug: bool = False, debug_interval: int=100000
) -> float:
    """
    Compute the probability of a graph fragmenting when randomly removing every node with a given probability.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    removal_probability: float
        The probability to to remove every node, a float between 0 and 1.
    iterations: int
        The number of iterations used for the estiamtion
    process_number: int, optional
        The number of processes to use for the estimation
    debug: bool, optional
        Whether to print debug information
    debug_interval: int, optional
        The interval at which to print debug information

    Returns
    -------
    float
        The estimated probability of the graph fragmenting.
    """
    pfrag, n = get_fragmentation_probability(
        G,
        iterations,
        probability_fragment,
        fragment_settings={
            "fragmentation": removal_probability,
            "fragmentation_type": "nodes",
        },
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pfrag


def get_fragmentation_probability_random_edge_removal(
    G: nx.Graph, removal_probability: float, iterations: int, process_number : int=1,debug: bool = False, debug_interval: int=100000
) -> float:
    """
    Compute the probability of a graph fragmenting when randomly removing every edge with a given probability.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    removal_probability : float
        The probability to to remove every edge, a float between 0 and 1.
    iterations : int
        The number of iterations used for the estiamtion
    process_number: int, optional
        The number of processes to use for the estimation
    debug: bool, optional
        Whether to print debug information
    debug_interval: int, optional
        The interval at which to print debug information

    Returns
    -------
    float
        The estimated probability of the graph fragmenting.
    """
    pfrag, n = get_fragmentation_probability(
        G,
        iterations,
        probability_fragment,
        fragment_settings={
            "fragmentation": removal_probability,
            "fragmentation_type": "edges",
        },
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pfrag


def get_fragmentation_probability_strength_node_removal(
    G: nx.Graph, removed_strength: float, iterations: int, process_number: int = 1, debug: bool = False, debug_interval: int=100000
) -> float:
    """
    Compute the probability of a graph fragmenting when removing random nodes util a fraction of the graph "strength" has been removed.
    Each node has a probability weight proportional to the inverse of the sum of the `strength` attributes of the edges connected to it.
    The strength of the graph is defined as the sum of the strength of all edges.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    removed_strength : float
        The fraction of the graph stength to remove, a float between 0 and 1.
    iterations : int
        The number of iterations used for the estiamtion
    process_number: int, optional
        The number of processes to use for the estimation
    debug: bool, optional
        Whether to print debug information
    debug_interval: int, optional
        The interval at which to print debug information

    Returns
    -------
    float
        The estimated probability of the graph fragmenting.

    Notes
    -----
    The strength of neighbouring nodes is not updated when a node is removed.
    """
    G_ = G.copy()
    init_nodes_strength(G_)
    pfrag, n = get_fragmentation_probability(
        G_,
        iterations,
        strength_nodes_fragment,
        fragment_settings={
            "fragmentation": removed_strength,
        },
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pfrag


def get_fragmentation_probability_strength_edge_removal(
    G: nx.Graph, removed_strength: float, iterations: int, process_number: int = 1, debug: bool = False, debug_interval: int=100000
) -> float:
    """
    Compute the probability of a graph fragmenting when removing random edges util a fraction of the graph "strength" has been removed.
    Each edge has a probability weight proportional to the inverse of its `strength` attribute.
    The strength of the graph is defined as the sum of the strength of all edges.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    removed_strength : float
        The fraction of the graph strength to remove, a float between 0 and 1.
    iterations : int
        The number of iterations used for the estiamtion
    process_number: int, optional
        The number of processes to use for the estimation
    debug: bool, optional
        Whether to print debug information
    debug_interval: int, optional
        The interval at which to print debug information

    Returns
    -------
    float
        The estimated probability of the graph fragmenting.
    """
    pfrag, n = get_fragmentation_probability(
        G,
        iterations,
        strength_edges_fragment,
        fragment_settings={
            "fragmentation": removed_strength,
        },
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pfrag


def get_fragmentation_probability_threshold_node(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    process_number: int = 1,
    debug: bool = False,
    debug_interval: int = 100000,
) -> Tuple[float, int]:
    """
    Estimate the probability of node removal that will fragment the graph with a probability of 1/2.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    error_probability : float
        An upper bound for the probability that the returned value is incorrect, a float between 0 and 1.
    steps : int
        The maximum number of steps for the bisection process. If correct the returned value will be within 2**(-steps) of the true value.
    min_iterations : int, optional
        The minimum number of iterations for each bisection step
    max_iterations : int, optional
        The maximum number of iterations for each bisection step
    process_number: int, optional
        The number of processes to use for the estimation
    debug : bool, optional
        If True, print debug information
    debug_interval : int, optional
        The interval at which to print debug information

    Returns
    -------
    Tuple[float, int]
        The estimated probability of node removal and the number of step reached.
    """
    pf, n = bisection(
        G,
        steps,
        error_probability,
        probability_fragment,
        fragment_settings={"fragmentation_type": "nodes"},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pf, n


def get_fragmentation_probability_threshold_edge(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    process_number: int = 1,
    debug: bool = False,
    debug_interval: int = 100000,
) -> Tuple[float, int]:
    """
    Estimate the probability of edge removal that will fragment the graph with a probability of 1/2.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    error_probability : float
        An upper bound for the probability that the returned value is incorrect, a float between 0 and 1.
    steps : int
        The maximum number of steps for the bisection process. If correct the returned value will be within 2**(-steps) of the true value.
    min_iterations : int, optional
        The minimum number of iterations for each bisection step
    max_iterations : int, optional
        The maximum number of iterations for each bisection step
    process_number: int, optional
        The number of processes to use for the estimation
    debug: bool, optional
        Whether to print debug information
    debug_interval: int, optional
        The interval at which to print debug information

    Returns
    -------
    Tuple[float, int]
        The estimated probability of edge removal and the number of step reached.
    """
    pf, n = bisection(
        G,
        steps,
        error_probability,
        probability_fragment,
        fragment_settings={"fragmentation_type": "edges"},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pf, n


def get_fragmentation_strength_threshold_edge(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    process_number: int = 1,
    debug: bool = False,
    debug_interval: int = 100000,
) -> Tuple[float, int]:
    """
    Estimate the fraction of the graph strength that needs to be removed (by randomly removing edges) to fragment the graph with a probability of 1/2.
    The probability weight of edges is inversly proportional to its `strength` attribute.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    error_probability : float
        An upper bound for the probability that the returned value is incorrect, a float between 0 and 1.
    steps : int
        The maximum number of steps for the bisection process. If correct the returned value will be within 2**(-steps) of the true value.
    min_iterations : int, optional
        The minimum number of iterations for each bisection step
    max_iterations : int, optional
        The maximum number of iterations for each bisection step
    process_number: int, optional
        The number of processes to use for the estimation
    debug: bool, optional
        Whether to print debug information
    debug_interval: int, optional
        The interval at which to print debug information

    Returns
    -------
    Tuple[float, int]
        The estimated strength to remove and the number of step reached.
    """
    pf, n = bisection(
        G,
        steps,
        error_probability,
        strength_edges_fragment,
        fragment_settings={},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pf, n


def get_fragmentation_strength_threshold_node(
    G: nx.Graph,
    error_probability: float,
    steps: int,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    process_number: int = 1,
    debug: bool = False,
    debug_interval: int = 100000,
) -> Tuple[float, int]:
    """
    Estimate the fraction of the graph strength that needs to be removed (by randomly removing nodes) to fragment the graph with a probability of 1/2.
    The probability weight of nodes is inversly proportional to its strength which is defined as the sum of the strength of all edges connected to the node.

    Parameters
    ----------
    G : nx.Graph
        The graph to analyse
    error_probability : float
        An upper bound for the probability that the returned value is incorrect, a float between 0 and 1.
    steps : int
        The maximum number of steps for the bisection process. If correct the returned value will be within 2**(-steps) of the true value.
    min_iterations : int, optional
        The minimum number of iterations for each bisection step
    max_iterations : int, optional
        The maximum number of iterations for each bisection step
    process_number: int, optional
        The number of processes to use for the estimation
    debug: bool, optional
        Whether to print debug information
    debug_interval: int, optional
        The interval at which to print debug information

    Returns
    -------
    Tuple[float,int]
        The estimated strength to remove and the number of step reached.
    """
    G_ = G.copy()
    init_nodes_strength(G_)
    pf, n = bisection(
        G_,
        steps,
        error_probability,
        strength_nodes_fragment,
        fragment_settings={},
        min_iterations=min_iterations,
        max_iterations=max_iterations,
        debug=debug,
        debug_interval=debug_interval,
        process_number=process_number
    )
    return pf, n
