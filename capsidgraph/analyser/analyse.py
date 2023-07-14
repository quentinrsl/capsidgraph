import time
import networkx as nx
from typing import List, Tuple, Dict, Callable
from inspect import signature


def _bisection_stop_condition(n: int, pfrag: float, settings: Dict) -> bool:
    INV_PROBA = 1 / settings["error_probability"]
    max_iterations = settings.get("max_iterations", 1000000)
    min_iterations = settings.get("min_iterations", 1000)
    return n >= min_iterations and (
        4 * n * ((pfrag - 0.5) ** 2) >= INV_PROBA or n >= max_iterations
    )


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
            1000 * (time.time() - start) / n,
            "ms/sim",
        )
    return pfrag, n


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
