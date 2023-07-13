import time
import networkx as nx
from typing import List, Tuple, Dict, Callable
from inspect import signature


def bisection_stop_condition(n: int, pfrag: float, settings: Dict) -> bool:
	INV_PROBA = 1 / settings.get("error_probability", 0.05)
	max_iterations = settings.get("max_iterations", 1000000)
	min_iterations = settings.get("min_iterations", 1000)
	return n >= min_iterations and (
		4 * n * ((pfrag - 0.5) ** 2) >= INV_PROBA or n >= max_iterations
	)


def get_framentation_probability(
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
	while (type(stop_condition) == int and n < stop_condition) or (
		callable(stop_condition) 
		and (not stop_condition(n, pfrag, stop_condition_settings))
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
