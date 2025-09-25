"""
Ecological graph theory methods, containing the functions to simulate species competition on a graph.
"""
import numpy as np
from copy import deepcopy
import random as rnd
import networkx as nx
import pandas as pd
from itertools import product, combinations


def get_prop_neighbor(graph):
    """Get the proportion of conspecific and allospecific neighbors.  

    Args:
        graph (NetworkX Graph)

    Returns:
        dictionary
    """
    neighbor_proportions = {}
    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))
        neighbors_with_attribute = sum([1 for neighbor in neighbors if graph.nodes[neighbor]['Phenotype'] != graph.nodes[node]['Phenotype']])
        total = len(neighbors) 
        proportion = (neighbors_with_attribute, total - neighbors_with_attribute, total)
        neighbor_proportions[node] = proportion
    return neighbor_proportions

def uniform_rewiring(G, m):
    """Rewire graph G according to uniform rewiring algorithm.

    Args:
        G (NetworkX Graph)
        m (int): The number of rewirings.

    Returns:
        NetworkX Graph
    """
    edges = list(G.edges)
    for _ in range(min(m, len(edges))):
        edge_to_remove = rnd.choice(edges)
        edges.remove(edge_to_remove)
        node_to_rewire = rnd.choice(edge_to_remove)
        current_neighbors = set(G.neighbors(node_to_rewire))
        current_neighbors.add(node_to_rewire) # to avoid self-loop 
        potential_new_nodes = set(G.nodes) - current_neighbors
        G.remove_edge(*edge_to_remove)
        if potential_new_nodes:
            new_node = rnd.choice(list(potential_new_nodes))   
            G.add_edge(node_to_rewire, new_node)
    return G

def RR_rewiring(G, m):
    """Rewire graph G according to rich-gets-richer rewiring algorithm.

    Args:
        G (NetworkX Graph)
        m (int): The number of rewirings.

    Returns:
        NetworkX Graph
    """
    edges = list(G.edges)
    for _ in range(min(m, len(edges))):
        edge_to_remove = rnd.choice(edges)
        edges.remove(edge_to_remove)
        node_to_rewire = rnd.choice(edge_to_remove)
        current_neighbors = set(G.neighbors(node_to_rewire))
        current_neighbors.add(node_to_rewire)  # to avoid self-loop
        potential_new_nodes = set(G.nodes) - current_neighbors
        G.remove_edge(*edge_to_remove)
        if potential_new_nodes:
            degrees = np.array([G.degree(node) for node in potential_new_nodes])
            total_degree = degrees.sum()
            if total_degree > 0:
                probabilities = degrees / total_degree
                new_node = rnd.choices(list(potential_new_nodes), weights=probabilities, k=1)[0]
            else:
                new_node = rnd.choice(list(potential_new_nodes))
            G.add_edge(node_to_rewire, new_node)
    return G

def PR_rewiring(G, m):
    """Rewire graph G according to poor-gets-rich rewiring algorithm.

    Args:
        G (NetworkX Graph)
        m (int): The number of rewirings.

    Returns:
        NetworkX Graph
    """
    edges = list(G.edges)
    for _ in range(min(m, len(edges))):
        edge_to_remove = rnd.choice(edges)
        edges.remove(edge_to_remove)
        node_to_rewire = rnd.choice(edge_to_remove)
        current_neighbors = set(G.neighbors(node_to_rewire))
        current_neighbors.add(node_to_rewire)  # to avoid self-loop
        potential_new_nodes = set(G.nodes) - current_neighbors
        G.remove_edge(*edge_to_remove)
        if potential_new_nodes:
            degrees = np.array([G.degree(node) for node in potential_new_nodes])
            inv_degress = np.max(degrees) - degrees
            total_inv_degress = inv_degress.sum()
            if total_inv_degress > 0:
                probabilities = inv_degress / total_inv_degress
                new_node = rnd.choices(list(potential_new_nodes), weights=probabilities, k=1)[0]
            else:
                new_node = rnd.choice(list(potential_new_nodes))
            G.add_edge(node_to_rewire, new_node)
    
    return G

def rnd_rewiring_agg(G, m, pref={'A':{'A': 4, 'B': 1}, 'B':{'A': 1, 'B': 4}}):
    """Rewire graph G according to aggregator or repulsive rewiring algorithms.

    Args:
        G (NetworkX Graph)
        m (int): The number of rewirings.
        pref (dict, optional): The preference to attach to conspecific and allospecific neighbors.

    Returns:
        NetworkX Graph
    """
    edges = list(G.edges)
    for _ in range(min(m, len(edges))):
        edge_to_remove = rnd.choice(edges)
        edges.remove(edge_to_remove)
        node_to_rewire = rnd.choice(edge_to_remove)
        current_neighbors = set(G.neighbors(node_to_rewire))
        current_neighbors.add(node_to_rewire)  # to avoid self-loop
        potential_new_nodes = set(G.nodes) - current_neighbors
        P_node_to_rewire = G.nodes[node_to_rewire]['Phenotype']
        G.remove_edge(*edge_to_remove)
        if potential_new_nodes:
            prefs = np.array([pref[P_node_to_rewire][G.nodes[j]['Phenotype']] for j in potential_new_nodes])            
            total_prefs = prefs.sum()
            if total_prefs > 0:
                probabilities = prefs / total_prefs
                new_node = rnd.choices(list(potential_new_nodes), weights=probabilities, k=1)[0]
            else:
                new_node = rnd.choice(list(potential_new_nodes))
            G.add_edge(node_to_rewire, new_node)    
    return G

def SC_on_dodec_graph(t, alpha = {"A":0.0, "B":0.0}):
    """Species competition on a dodecahedron graph.

    Args:
        t (int): Number of simulation steps.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
    """
    dodec_pop = nx.dodecahedral_graph()
    pop_size = len(dodec_pop.nodes)
    nodes = list(dodec_pop.nodes)
    species_A = np.random.choice(range(pop_size), size=pop_size//2, replace=False)
    attributes = {nodes[i]:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in range(pop_size)}
    nx.set_node_attributes(dodec_pop, attributes)
    init_pop = deepcopy(dodec_pop)
    freq_B = []
    freq_B.append([dodec_pop.nodes[j]['Phenotype'] for j in dodec_pop.nodes].count('B'))
    for i in range(t):
        # density-dependent death
        prop_neigh = get_prop_neighbor(dodec_pop)
        death_rate = [np.divide(v[1] + alpha[dodec_pop.nodes[k]['Phenotype']]*v[0], v[2]) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(dodec_pop.nodes, p=normalized_d)
        # random birth
        neighbors = [n for n in dodec_pop.neighbors(rand_ind)]
        if len(neighbors):
            selected_neigh = np.random.choice(neighbors)
            dodec_pop.nodes[rand_ind]['Phenotype'] = dodec_pop.nodes[selected_neigh]['Phenotype']
            dodec_pop.nodes[rand_ind]['W'] = dodec_pop.nodes[selected_neigh]['W']
        freq_B.append([dodec_pop.nodes[j]['Phenotype'] for j in dodec_pop.nodes].count('B'))
    return init_pop, dodec_pop, freq_B

def create_kings_graph_with_probability(n, p):
    """Creates an n x n King's graph with edges added with probability p.

    Args:
        n (int): Size of the chessboard (n x n).
        p (float): Probability of adding an edge between two nodes.

    Returns:
        NetworkX graph 
    """
    G = nx.Graph()
    for i in range(n):
        for j in range(n):
            current = (i, j)
            G.add_node(current)
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    if dx == 0 and dy == 0:
                        continue
                    neighbor = (i + dx, j + dy)
                    if 0 <= neighbor[0] < n and 0 <= neighbor[1] < n:
                        if rnd.random() < p:
                            G.add_edge(current, neighbor)
    return G

def SC_on_kings_graph(N, p, t, alpha = {"A":0.0, "B":0.0}, t_series=False): 
    """Species competition on a King's graph.

    Args:
        N (int): Population dimension, resulting in n x n nodes.
        p (float): Probability of adding an edge between two nodes.
        t (int): Number of simulation steps.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
        t_series (bool, optional): If True, return the graphs during the simulation. Defaults to False.
    """
    pop = create_kings_graph_with_probability(N, p)
    nodes = list(pop.nodes)
    pop_size = len(nodes)
    species_A = np.random.choice(range(pop_size), size=pop_size//2, replace=False)
    species_A = sorted(species_A)
    attributes = {nodes[i]:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in range(pop_size)}
    nx.set_node_attributes(pop, attributes)
    init_pop = deepcopy(pop)
    freq_B = []
    freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    pop_t_series = []
    pop_t_series.append(init_pop)
    for _ in range(t):
        # density-dependent death
        prop_neigh = get_prop_neighbor(pop)
        death_rate = [np.divide(v[1] + alpha[pop.nodes[k]['Phenotype']]*v[0], v[2]) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(range(pop_size), p=normalized_d)
        # random birth
        neighbours = [n for n in pop.neighbors(nodes[rand_ind])]
        selected_neigh = np.random.choice(range(len(neighbours)))
        selected_neigh = neighbours[selected_neigh]
        pop.nodes[nodes[rand_ind]]['Phenotype'] = pop.nodes[selected_neigh]['Phenotype']
        pop.nodes[nodes[rand_ind]]['W'] = pop.nodes[selected_neigh]['W']
        freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
        current_pop = deepcopy(pop)
        if t_series:
            pop_t_series.append(current_pop)
    if t_series:
        return pop_t_series, freq_B
    else:
        return init_pop, current_pop, freq_B

def SC_on_comp_graph(N, t, alpha = {"A":0.0, "B":0.0}, t_series=False):
    """Species competition on a complete graph.

    Args:
        N (int): Population size.
        t (int): Number of simulation steps.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
        t_series (bool, optional): If True, return the graphs during the simulation. Defaults to False.
    """
    pop = nx.complete_graph(N)
    nodes = list(pop.nodes)
    species_A = np.random.choice(range(N), size=N//2, replace=False)
    species_A = sorted(species_A)
    attributes = {nodes[i]:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in range(N)}
    nx.set_node_attributes(pop, attributes)
    init_pop = deepcopy(pop)
    freq_B = []
    freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    pop_t_series = []
    pop_t_series.append(init_pop)
    for _ in range(t):
        # density-dependent death
        prop_neigh = get_prop_neighbor(pop)
        death_rate = [np.divide(v[1] + alpha[pop.nodes[k]['Phenotype']]*v[0], v[2]) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(range(N), p=normalized_d)
        # random birth
        neighbours = [n for n in pop.neighbors(nodes[rand_ind])]
        selected_neigh = np.random.choice(range(len(neighbours)))
        selected_neigh = neighbours[selected_neigh]
        pop.nodes[nodes[rand_ind]]['Phenotype'] = pop.nodes[selected_neigh]['Phenotype']
        pop.nodes[nodes[rand_ind]]['W'] = pop.nodes[selected_neigh]['W']
        freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
        current_pop = deepcopy(pop)
        if t_series:
            pop_t_series.append(current_pop)
    if t_series:
        return pop_t_series, freq_B
    else:
        return init_pop, current_pop, freq_B
    
def SC_on_NWSgraph(N, p, k, t, alpha = {"A":0.0, "B":0.0}):
    """Species competition on a Newman-Watts-Strogatz graph.

    Args:
        N (int): Population size.
        p (float): Probability of adding new edges.
        k (int): The initial degree of each node.
        t (int): Number of simulation steps.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
    """
    pop = nx.newman_watts_strogatz_graph(N, k, p)
    species_A = np.random.choice(pop.nodes, size=N//2, replace=False)
    attributes = {i:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in pop.nodes }
    nx.set_node_attributes(pop, attributes)
    init_pop = deepcopy(pop)
    freq_B = []
    freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    for i in range(t):
        # density-dependent death
        prop_neigh = get_prop_neighbor(pop)

        death_rate = [np.divide(v[1] + alpha[pop.nodes[k]['Phenotype']]*v[0], v[2]+ 1e-5) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        if total != 0:
            normalized_d = [p / total for p in death_rate]
            rand_ind = np.random.choice(pop.nodes, p=normalized_d)
        else:
            rand_ind = np.random.choice(pop.nodes)
        # random birth
        neighbours = [n for n in pop.neighbors(rand_ind)]
        if len(neighbours):
            selected_neigh = np.random.choice(neighbours)
            pop.nodes[rand_ind]['Phenotype'] = pop.nodes[selected_neigh]['Phenotype']
            pop.nodes[rand_ind]['W'] = pop.nodes[selected_neigh]['W']
        freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    return init_pop, pop, freq_B

def SC_on_ER_graph(N, p, t, alpha = {"A":0.0, "B":0.0}, t_series=False):
    """Species competition on an Erdős–Rényi graph.

    Args:
        N (int): Population size.
        p (float): Probability for edge creation.
        t (int): Number of simulation steps.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
        t_series (bool, optional): If True, return the graphs during the simulation. Defaults to False.
    """
    pop = nx.erdos_renyi_graph(N, p)
    nodes = list(pop.nodes)
    species_A = np.random.choice(range(N), size=N//2, replace=False)
    species_A = sorted(species_A)
    attributes = {nodes[i]:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in range(N)}
    nx.set_node_attributes(pop, attributes)
    init_pop = deepcopy(pop)
    freq_B = []
    freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    pop_t_series = []
    pop_t_series.append(init_pop)
    for _ in range(t):
        # density-dependent death
        prop_neigh = get_prop_neighbor(pop)
        death_rate = [np.divide(v[1] + alpha[pop.nodes[k]['Phenotype']]*v[0], v[2]+ 1e-5) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(range(N), p=normalized_d)
        # random birth
        neighbours = [n for n in pop.neighbors(nodes[rand_ind])]
        if len(neighbours):
            selected_neigh = np.random.choice(range(len(neighbours)))
            selected_neigh = neighbours[selected_neigh]
            pop.nodes[nodes[rand_ind]]['Phenotype'] = pop.nodes[selected_neigh]['Phenotype']
            pop.nodes[nodes[rand_ind]]['W'] = pop.nodes[selected_neigh]['W']
        freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
        current_pop = deepcopy(pop)
        if t_series:
            pop_t_series.append(current_pop)
    if t_series:
        return pop_t_series, freq_B
    else:
        return init_pop, current_pop, freq_B
    
def SC_on_ER_graph_with_UR(N, p, t, m, alpha = {"A":0.0, "B":0.0}, t_series=False):
    """Species competition on an Erdős–Rényi graph with uniform rewiring algorithm.

    Args:
        N (int): Population size.
        p (float): Probability for edge creation.
        t (int): Number of simulation steps.
        m (int): Number of rewirings per step.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
        t_series (bool, optional): If True, return the graphs during the simulation. Defaults to False.
    """
    pop = nx.erdos_renyi_graph(N, p)
    nodes = list(pop.nodes)
    species_A = np.random.choice(range(N), size=N//2, replace=False)
    species_A = sorted(species_A)
    attributes = {nodes[i]:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in range(N)}
    nx.set_node_attributes(pop, attributes)
    init_pop = deepcopy(pop)
    freq_B = []
    freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    pop_t_series = []
    pop_t_series.append(init_pop)
    for _ in range(t):
        # Random rewiring
        pop = uniform_rewiring(pop, m)
        # density-dependent death
        prop_neigh = get_prop_neighbor(pop)
        death_rate = [np.divide(v[1] + alpha[pop.nodes[k]['Phenotype']]*v[0], v[2]+ 1e-5) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(range(N), p=normalized_d)
        # random birth
        neighbours = [n for n in pop.neighbors(nodes[rand_ind])]
        if len(neighbours):
            selected_neigh = np.random.choice(range(len(neighbours)))
            selected_neigh = neighbours[selected_neigh]
            pop.nodes[nodes[rand_ind]]['Phenotype'] = pop.nodes[selected_neigh]['Phenotype']
            pop.nodes[nodes[rand_ind]]['W'] = pop.nodes[selected_neigh]['W']
        freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
        current_pop = deepcopy(pop)
        if t_series:
            pop_t_series.append(current_pop)
    if t_series:
        return pop_t_series, freq_B
    else:
        return init_pop, current_pop, freq_B

def SC_on_ER_graph_with_PA(N, p, t, m, alpha = {"A":0.0, "B":0.0}, RR_rewire=True, t_series=False):
    """Species competition on an Erdős–Rényi graph with preferential-attachment rewiring algorithms.

    Args:
        N (int): Population size.
        p (float): Probability for edge creation.
        t (int): Number of simulation steps.
        m (int): Number of rewirings per step.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
        RR_rewire (bool, optional): If True, use rich-gets-richer algorithm, otherwise poor-gets-rich algorithm. Defaults to True.
        t_series (bool, optional): If True, return the graphs during the simulation. Defaults to False.
    """
    pop = nx.erdos_renyi_graph(N, p)
    nodes = list(pop.nodes)
    species_A = np.random.choice(range(N), size=N//2, replace=False)
    species_A = sorted(species_A)
    attributes = {nodes[i]:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in range(N)}
    nx.set_node_attributes(pop, attributes)
    init_pop = deepcopy(pop)
    freq_B = []
    freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    pop_t_series = []
    pop_t_series.append(init_pop)
    for _ in range(t):
        # Random rewiring
        if RR_rewire:
            pop = RR_rewiring(pop, m)
        else:
            pop = PR_rewiring(pop, m)
        # density-dependent death
        prop_neigh = get_prop_neighbor(pop)
        death_rate = [np.divide(v[1] + alpha[pop.nodes[k]['Phenotype']]*v[0], v[2]+ 1e-5) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(range(N), p=normalized_d)
        # random birth
        neighbours = [n for n in pop.neighbors(nodes[rand_ind])]
        if len(neighbours):
            selected_neigh = np.random.choice(range(len(neighbours)))
            selected_neigh = neighbours[selected_neigh]
            pop.nodes[nodes[rand_ind]]['Phenotype'] = pop.nodes[selected_neigh]['Phenotype']
            pop.nodes[nodes[rand_ind]]['W'] = pop.nodes[selected_neigh]['W']
        freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
        current_pop = deepcopy(pop)
        if t_series:
            pop_t_series.append(current_pop)
    if t_series:
        return pop_t_series, freq_B
    else:
        return init_pop, current_pop, freq_B
    
def SC_on_ER_graph_agg(N, p, t, m, alpha = {"A":0.0, "B":0.0}, pref={'A':{'A': 1, 'B': 1}, 'B':{'A': 1, 'B': 1}}, t_series=False):
    """Species competition on an Erdős–Rényi graph with aggregator or repulsive rewiring algorithm.

    Args:
        N (int): Population size.
        p (float): Probability for edge creation.
        t (int): Number of simulation steps.
        m (int): Number of rewirings per step.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
        pref (dict, optional): Rewiring parameters. For aggregator, {'A':{'A': 4, 'B': 1}, 'B':{'A': 1, 'B': 4}}.
                               For repulsive, {'A':{'A': 1, 'B': 4}, 'B':{'A': 4, 'B': 1}.
        t_series (bool, optional): If True, return the graphs during the simulation. Defaults to False.
    """
    pop = nx.erdos_renyi_graph(N, p)
    nodes = list(pop.nodes)
    species_A = np.random.choice(range(N), size=N//2, replace=False)
    species_A = sorted(species_A)
    attributes = {nodes[i]:{'Phenotype': 'A', 'W' : 1.} if i in species_A else {'Phenotype': 'B', 'W' : 1.} for i in range(N)}
    nx.set_node_attributes(pop, attributes)
    init_pop = deepcopy(pop)
    freq_B = []
    freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
    pop_t_series = []
    pop_t_series.append(init_pop)
    for _ in range(t):
        # Random rewiring
        pop = rnd_rewiring_agg(pop, m, pref)
        # density-dependent death
        prop_neigh = get_prop_neighbor(pop)
        death_rate = [np.divide(v[1] + alpha[pop.nodes[k]['Phenotype']]*v[0], v[2]+ 1e-5) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(range(N), p=normalized_d)
        # random birth
        neighbours = [n for n in pop.neighbors(nodes[rand_ind])]
        if len(neighbours):
            selected_neigh = np.random.choice(range(len(neighbours)))
            selected_neigh = neighbours[selected_neigh]
            pop.nodes[nodes[rand_ind]]['Phenotype'] = pop.nodes[selected_neigh]['Phenotype']
            pop.nodes[nodes[rand_ind]]['W'] = pop.nodes[selected_neigh]['W']
        freq_B.append([pop.nodes[j]['Phenotype'] for j in pop.nodes].count('B'))
        current_pop = deepcopy(pop)
        if t_series:
            pop_t_series.append(current_pop)
    if t_series:
        return pop_t_series, freq_B
    else:
        return init_pop, current_pop, freq_B
    
### metapop methods

def creat_edges_btwn_pops(pop1, pop2, β):
    """Create edges between populations in a metapopulation (or metacommunity).

    Args:
        pop1 (int): pop1 index.
        pop2 (int): pop2 index.

    Returns:
        list: list of edges.
    """
    pairs = list(product(pop1, pop2))
    total_pairs = len(pairs)
    num_pairs = np.random.poisson(β)
    num_pairs = min(num_pairs, total_pairs)
    sampled_pairs = np.random.choice(range(total_pairs), num_pairs, replace=False)
    return [pairs[i] for i in sampled_pairs]

def splice_list(lst, n):
    """Splice a list in intervals of length n.
    """
    return [lst[i:i + n] for i in range(0, len(lst), n)]
    
def SC_on_static_metapop(n_comm, N_pop, t, omega, p_ER, alpha = {"A":0.0, "B":0.0}):
    """Species competition on a metapopulation (or metacommunity).

    Args:
        n_comm (int): Number of communities.
        N_pop (int): Population size in each subpopulation.
        t (int): Number of simulation steps.
        omega (float): Probability of edge creation between communities.
        p_ER (float): Probability for edge creation in Erdős–Rényi model.
        alpha (dict, optional): Species-interaction parameters. Defaults to {"A":0.0, "B":0.0}.
    """
    # Specify nodes in each subpopulation
    metapop = nx.Graph()
    attributes = []
    for pop in range(n_comm):
        species_A = np.random.choice(range(N_pop), size=N_pop//2, replace=False)
        attributes_temp = [{'Phenotype': 'A', 'W' : 1., 'subpop': pop} if i in species_A else {'Phenotype': 'B', 'W' : 1., 'subpop': pop}  for i in range(N_pop)]
        attributes.extend(attributes_temp)
    # Create random edges in each subpop baesd on the ER method
    indices = np.arange(0, int(n_comm*N_pop))
    np.random.shuffle(indices)
    supops_ind = splice_list(indices, N_pop)
    poss_edges = [list(combinations(i, 2)) for i in supops_ind]
    poss_edges = [item for sublist in poss_edges for item in sublist]
    L = int(np.divide(N_pop*(N_pop - 1), 1))
    real_edges_prob = [list(np.random.binomial(1, p_ER, L)) for i in range(n_comm)]
    real_edges_prob = [item for sublist in real_edges_prob for item in sublist]
    real_edges = [i for i,j in zip(poss_edges, real_edges_prob) if j]
    # Create a metapop graph contaning all the communities
    metapop = nx.Graph()
    metapop.add_nodes_from(indices)
    attributes = {i:j for i,j in zip(indices, attributes)}
    nx.set_node_attributes(metapop, attributes)
    metapop.add_edges_from(real_edges)
    # create random edges between each pair of communities according to parameter β
    pop_combs = list(combinations(np.arange(0,n_comm,1), 2))
    new_edges = []
    for i in pop_combs:
        popA = supops_ind[i[0]]
        popB = supops_ind[i[1]]
        AB_new_edges = creat_edges_btwn_pops(popA, popB, omega)
        new_edges.extend(AB_new_edges)
    metapop.add_edges_from(new_edges)
    f0 = deepcopy(metapop)
    metapop_size = len(metapop.nodes)
    nodes = list(metapop.nodes)
    freq_B = pd.DataFrame()
    freq_per_pop = [[metapop.nodes[i]['Phenotype'] for i in j].count('B') for j in supops_ind]
    freq_total = [metapop.nodes[j]['Phenotype'] for j in nodes].count('B')
    temp = {'sub_' + str(i):[np.round(j/N_pop, decimals=3)] for i,j in zip(range(N_pop), freq_per_pop)}
    temp['tot'] = [np.round(freq_total/metapop_size, decimals=3)]
    temp_DF = pd.DataFrame.from_dict(temp)
    freq_B = pd.concat([freq_B, temp_DF])
    for _ in range(t):
        prop_neigh = get_prop_neighbor(metapop)
        death_rate = [np.divide(v[1] + alpha[metapop.nodes[k]['Phenotype']]*v[0], v[2]+ 1e-5) for i, (k,v) in enumerate(prop_neigh.items())]
        total = sum(death_rate)
        normalized_d = [p / total for p in death_rate]
        rand_ind = np.random.choice(range(metapop_size), p=normalized_d)
        # random birth
        neighbours = [n for n in metapop.neighbors(nodes[rand_ind])]
        if len(neighbours):
            selected_neigh = np.random.choice(range(len(neighbours)))
            selected_neigh = neighbours[selected_neigh]
            metapop.nodes[nodes[rand_ind]]['Phenotype'] = metapop.nodes[selected_neigh]['Phenotype']
            metapop.nodes[nodes[rand_ind]]['W'] = metapop.nodes[selected_neigh]['W']
        freq_per_pop = [[metapop.nodes[i]['Phenotype'] for i in j].count('B') for j in supops_ind]
        freq_total = [metapop.nodes[j]['Phenotype'] for j in nodes].count('B')
        temp = {'sub_' + str(i):[np.round(j/N_pop, decimals=3)] for i,j in zip(range(N_pop), freq_per_pop)}
        temp['tot'] = [np.round(freq_total/metapop_size, decimals=3)]
        temp_DF = pd.DataFrame.from_dict(temp)
        freq_B = pd.concat([freq_B, temp_DF])
    freq_B = freq_B.reset_index()
    return f0, metapop, freq_B