import time
import networkx as nx
from typing import List, Tuple, Dict, Callable
from inspect import signature
import random


def _get_framentation_probability(
    G: nx.Graph,
    stop_condition: int | Callable[[int, float, Dict], bool],
    fragment: Callable[[nx.Graph, Dict], None] | Callable[[nx.Graph], None],
    stop_condition_settings: Dict | None = None,
    fragment_settings: Dict | None = None,
    debug: bool = False,
) -> Tuple[float, bool]:
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
    INV_PROBA = 1 / settings["error_probability"]
    max_iterations = settings.get("max_iterations", 1000000)
    min_iterations = settings.get("min_iterations", 1000)
    return n >= min_iterations and (
        4 * n * ((pfrag - 0.5) ** 2) >= INV_PROBA or n >= max_iterations
    )


def _bisection(
    G: nx.Graph,
    steps: int,
    error_probability: float,
    fragment: Callable[[nx.Graph, Dict], None] | Callable[[nx.Graph], None],
    fragment_settings: Dict | None = None,
    min_iterations: int = 1000,
    max_iterations: int = 1000000,
    debug: bool = False,
) -> Tuple[float, int]:
    # Compute the upper bond of error for one bisection step from the upper bond of making a mistake in the entire process
    eps = 1 - (1 - error_probability) ** (1 / steps)
    lower_bound = 0
    upper_bound = 1
    step_count = 0
    while step_count < steps:
        middle = (lower_bound + upper_bound) / 2
        step_count += 1
        fragment_settings["fragmentation"] = middle
        pfrag, iteration_count = _get_framentation_probability(
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


# G is the graph to percolate on
# p is the probability of removal of each nodes/edges
# n is the number of iterations
# fragtype is either "nodes" or "edges" and determine if edges or nodes are to be removed
# Returns an array A where A[i] is the number of fragments with size i encountered
def get_fragment_size_distribution(
    G: nx.Graph,
    iterations: int,
    fragment: Callable[[nx.Graph, Dict], None] | Callable[[nx.Graph], None],
    fragment_settings: Dict | None = None,
) -> List[float]:
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

## Unweighted graphs


def _get_hole_size(fragmneted_graph, original_graph):
    connected_components = list(nx.connected_components(fragmneted_graph))
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
    fragment: Callable[[nx.Graph, Dict], None] | Callable[[nx.Graph], None],
    fragment_settings: Dict | None = None,
) -> List[float]:
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
        for i in range(max(holes_size.keys()) + 1)
    ]
