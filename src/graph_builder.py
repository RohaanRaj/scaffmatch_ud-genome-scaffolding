import networkx as nx
from tqdm import tqdm

def build_graph(kmer_counts, k):

    G = nx.DiGraph()

    mask = (1 << (2*(k-1))) - 1

    for kmer, count in tqdm(kmer_counts.items(), desc="Building graph"):

        if count < 2:
            continue

        prefix = kmer >> 2
        suffix = kmer & mask

        G.add_edge(prefix, suffix, weight=count)

    return G
