import time
import networkx as nx
from typing import List, Tuple, Dict, Callable
from inspect import signature
import random


def get_framentation_probability(
    G: nx.Graph,
    stop_condition: int | Callable[[int, float, Dict], bool],
    fragment: Callable[[nx.Graph, Dict], nx.Graph] | Callable[[nx.Graph], nx.Graph],
    stop_condition_settings: Dict | None = None,
    fragment_settings: Dict | None = None,
    debug: bool = False,
) -> Tuple[float, int]:
    """
    Compute the fragmentation probability of a graph G using a given fragmentation method

    Parameters
    ----------
    G : nx.Graph
        The graph to fragment
    stop_condition : int | Callable[[int, float, Dict], bool]
        The stop condition for the fragmentation process. If an int is given, the process will stop after this number of iterations. If a callable is given, the process will stop when the callable returns True. The callable takes as parameters the number of iterations, the current estimated fragmentation probability and the stop_condition_settings
    fragment : Callable[[nx.Graph, Dict], None] | Callable[[nx.Graph], None]
        The fragmentation method to use. It must take as parameter a graph and a dict of settings and return the fragmented graph. If the fragmentation method does not take settings, it can be given without the settings parameter
    stop_condition_settings : Dict | None
        The settings to pass to the stop condition callable
    fragment_settings : Dict | None
        The settings to pass to the fragmentation method
    debug : bool
        If True, print debug information

    Returns
    -------
    Tuple[float, bool]
        The estimated fragmentation probability and an int representing the number of iterations used to compute it
    """
    start = time.time()
    fragmentation_count = 0
    pfrag = 0
    n = 0
    while (
        type(stop_condition) == int
        and n < stop_condition
        or (
            callable(stop_condition)
            and not stop_condition(n, pfrag, stop_condition_settings)
        )
    ):
        if len(signature(fragment).parameters) == 2:
            G_ = fragment(G, fragment_settings)
        else:
            G_ = fragment(G)

        n += 1
        if len(G_.nodes) > 0 and not nx.is_connected(G_):
            fragmentation_count += 1
        pfrag = fragmentation_count / n
    if debug:
        print(
            "fragmentation setttings=",
            fragment_settings,
            "stop_condition_settings=",
            stop_condition_settings,
            "with n=",
            n,
            "got p(frag)=",
            pfrag,
            000 * (time.time() - start) / n,
            "ms/sim",
        )
    return pfrag, n


def _bisection_stop_condition(n: int, pfrag: float, settings: Dict) -> bool:
    """
    Stop condition for the bisection method
    Returns True if 4*n*(pfrag-0.5)^2 >= 1/error_probability or n >= max_iterations
    This ensures that the probability of making an error in ont of the bissection step (estimating pfrag>1/2  when pfrag<1/2 or pfrag<1/2 when pfrag>1/2) is less than error_probability

    Parameters
    ----------
    n : int
        The number of iterations
    pfrag : float
        The current estimated fragmentation probability
    settings : Dict
        The settings of the bisection method\n
        The entry "error_probability" is the probability of making an error in one bisection step\n
        The entry "max_iterations" is the maximum number of iterations to perform before stopping the bisection method\n
        The entry "min_iterations" is the minimum number of iterations to perform

    Returns
    -------
    bool
        True if the stop condition is met, False otherwise
    """
    INV_PROBA = 1 / settings["error_probability"]
    max_iterations = settings.get("max_iterations", 1000000)
    min_iterations = settings.get("min_iterations", 1000)
    return n >= min_iterations and (
        4 * n * ((pfrag - 0.5) ** 2) >= INV_PROBA or n >= max_iterations
    )


def bisection(
    G: nx.Graph,
    steps: int,
    error_probability: float,
    fragment: Callable[[nx.Graph, Dict], nx.Graph],
    fragment_settings: Dict | None = None,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    debug: bool = False,
) -> Tuple[float, int]:
    """
    Compute the fragmentation threshold of a graph G using a given fragmentation method, ie the "fragmentation" parameter of the fragmentation method for which the graph is fragmented with probability 1/2

    Parameters
    ----------
    G : nx.Graph
        The graph to fragment
    steps : int
        The number of bisection steps to perform
    error_probability : float
        The probability of making an error during the entire bisection process
    fragment : Callable[[nx.Graph, Dict], None]
        The fragmentation method to use. It must take as parameter a graph and a dict of settings and return the fragmented graph.
        The fragmentation method must take a settings parameter which is a dict containing the "fragmentation" parameter, a value between 0 and 1
    fragment_settings : Dict | None
        The settings to pass to the fragmentation method, except for the "fragmentation", "min_iterations", "max_iterations" parameter
    min_iterations : int
        The minimum number of iterations to perform for each step
    max_iterations : int
        The maximum number of iterations to perform for each step

    Returns
    -------
    Tuple[float, int]
        The estimated fragmentation threshold and the number of steps reached

    """
    # Compute the upper bond of error for one bisection step from the upper bond of making a mistake in the entire process
    eps = 1 - (1 - error_probability) ** (1 / steps)
    lower_bound = 0
    upper_bound = 1
    step_count = 0
    while step_count < steps:
        middle = (lower_bound + upper_bound) / 2
        step_count += 1
        fragment_settings["fragmentation"] = middle
        pfrag, iteration_count = get_framentation_probability(
            G,
            _bisection_stop_condition,
            fragment,
            stop_condition_settings={
                "error_probability": eps,
                "min_iterations": min_iterations,
                "max_iterations": max_iterations,
            },
            fragment_settings=fragment_settings,
            debug=debug,
        )
        if iteration_count == max_iterations:
            return middle, step_count
        elif pfrag > 0.5:
            upper_bound = middle
        else:
            lower_bound = middle
    return middle, step_count


def get_fragment_size_distribution(
    G: nx.Graph,
    iterations: int,
    fragment: Callable[[nx.Graph, Dict], nx.Graph] | Callable[[nx.Graph], nx.Graph],
    fragment_settings: Dict | None = None,
) -> List[float]:
    """
    Compute the distribution of the size of the fragments obtained by fragmenting a graph G with a given fragmentation method

    Parameters
    ----------
    G : nx.Graph
        The graph to fragment
    iterations : int
        The number of iterations to perform
    fragment : Callable[[nx.Graph, Dict], None] | Callable[[nx.Graph], None]
        The fragmentation method to use. It must take as parameter a graph and a dict of settings and return the fragmented graph. The function may not take the settings parameter.
    fragment_settings : Dict | None
        The settings to pass to the fragmentation method.

    Returns
    -------
    List[float]
        The distribution of the size of the fragments obtained by fragmenting a graph G with a given fragmentation method.
        The entry i of the list is the estimated expected number of fragments of size i
    """
    fragments_size = {}
    max_size = 0
    for i in range(iterations):
        if len(signature(fragment).parameters) == 2:
            G_ = fragment(G, fragment_settings)
        else:
            G_ = fragment(G)
        for component in nx.connected_components(G_):
            size = len(component)
            max_size = max(size, max_size)
            fragments_size[size] = fragments_size.get(size, 0) + 1
    return [
        fragments_size[i] / iterations if i in fragments_size else 0
        for i in range(max_size + 1)
    ]


# ===================
# Hole Size Detection
# ===================

def _get_hole_size(fragmented_graph:nx.Graph, original_graph:nx.Graph)->int:
    """
    Compute the size of the hole in a fragmented graph.

    Parameters
    ----------
    fragmented_graph : nx.Graph
        The fragmented graph.
    original_graph : nx.Graph
        The graph that has been fragmented to give `fragmented_graph`
    
    Returns
    -------
    int
        the size of the hole size of the fragmented graph
    
    Notes
    -----
    To define a hole we consider the largest connected component of `fragmented_graph`. We then take the dual of that graph (ie. The graph contructed by taking out those nodes from `original_graph`).
    The largest connected component of that graph is the hole.
    In case of multiple largest conencted components for `gragmented_graph` we repeat the process for each of them and return the smallest value.
    """
    connected_components = list(nx.connected_components(fragmented_graph))
    if len(connected_components) == 0:
        return len(original_graph.nodes)
    largest_component_size = len(max(connected_components, key=len))
    largest_components = [
        c for c in connected_components if len(c) == largest_component_size
    ]
    min_hole_size = None
    for largest_component in largest_components:
        hole = original_graph.subgraph(original_graph.nodes - largest_component)
        if len(hole.nodes) == 0:
            return 0
        hole_size = len(max(nx.connected_components(hole), key=len))
        if min_hole_size == None or hole_size < min_hole_size:
            min_hole_size = hole_size
    return min_hole_size


def get_hole_size_distribution(
    G: nx.Graph,
    iterations: int,
    fragment: Callable[[nx.Graph, Dict], nx.Graph] | Callable[[nx.Graph], nx.Graph],
    fragment_settings: Dict | None = None,
) -> List[float]:
    """
    Compute the distribution of the size of the holes obtained by fragmenting a graph G with a given fragmentation method.

    Parameters
    ----------
    G : nx.Graph
        The graph to fragment
    iterations : int
        The number of iterations to perform
    fragment : Callable[[nx.Graph, Dict], None] | Callable[[nx.Graph], None]
        The fragmentation method to use, it must take as a parameter the graph to fragment, and may take a second Dict paramter containing settings.
    fragment_settings : Dict
        The settings to pass to the fragment method

    Returns
    -------
    List[float]:
        The probability distribution of hole sizes. The i-th entry of the list is the probability of obtaining a hole of size i.
    """
    # Initialize the list of hole sizes
    holes_size = {}
    # For each iteration
    for i in range(iterations):
        if len(signature(fragment).parameters) == 2:
            G_ = fragment(G, fragment_settings)
        else:
            G_ = fragment(G)
        m = _get_hole_size(G_, G)
        holes_size[m] = holes_size.get(m, 0) + 1
    return [
        holes_size[i] / iterations if i in holes_size else 0
        for i in range(len(G.nodes) + 1)
    ]
