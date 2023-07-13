import networkx as nx
from .analyse import get_framentation_probability, bisection
from .fragment import probability_fragment, energy_edges_fragment, energy_nodes_fragment
from typing import Tuple

def get_fragmentation_probability_random_node_removal(G:nx.Graph,removal_probability:float,iterations:int)->float:
    pfrag, n = get_framentation_probability(
        G,
        iterations,
        probability_fragment,
        fragment_settings={
            "fragmentation": removal_probability,
            "fragmentation_type": "nodes",
        },
    )
    return pfrag

def get_fragmentation_probability_random_edge_removal(G:nx.Graph,removal_probability:float,iterations:int)->float:
    pfrag, n = get_framentation_probability(
        G,
        iterations,
        probability_fragment,
        fragment_settings={
            "fragmentation": removal_probability,
            "fragmentation_type": "edges",
        },
    )
    return pfrag

#TODO: add tests for the wrappers
def get_fragmentation_probability_threshold_node(G:nx.Graph,error_probability:float, steps:int, min_iterations:int=1000, max_iterations:int=1000000, debug:bool=False) -> Tuple[float,int]:
    pf, n = bisection(G, steps, error_probability, probability_fragment, fragment_settings={"fragmentation_type": "nodes"}, min_iterations=min_iterations, max_iterations=max_iterations, debug=debug)
    return pf,n

def get_fragmentation_probability_threshold_edge(G:nx.Graph,error_probability:float, steps:int, min_iterations:int=1000, max_iterations:int=1000000, debug:bool=False) -> Tuple[float,int]:
    pf, n = bisection(G, steps, error_probability, probability_fragment, fragment_settings={"fragmentation_type": "edges"}, min_iterations=min_iterations, max_iterations=max_iterations, debug=debug)
    return pf,n

def get_fragmentation_energy_threshold_edge(G:nx.Graph, error_probability:float, steps:int, min_iterations:int=1000, max_iterations:int=1000000, debug:bool=False) -> Tuple[float,int]:
    pf, n = bisection(G, steps, error_probability, energy_edges_fragment, fragment_settings={}, min_iterations=min_iterations, max_iterations=max_iterations, debug=debug)
    return pf,n

def get_fragmentation_energy_threshold_node(G:nx.Graph, error_probability:float, steps:int, min_iterations:int=1000, max_iterations:int=1000000, debug:bool=False) -> Tuple[float,int]:
    pf, n = bisection(G, steps, error_probability, energy_nodes_fragment, fragment_settings={}, min_iterations=min_iterations, max_iterations=max_iterations, debug=debug)
    return pf,n