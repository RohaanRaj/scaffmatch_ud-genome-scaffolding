from read_loader import load_reads
from kmer_counter import count_kmers
from graph_builder import build_graph
from graph_cleaner import remove_tips, remove_bubbles
from contig_extractor import extract_contigs
from metrics import compute_metrics
from plotter import plot_lengths

ENABLE_GRAPH_VISUALIZATION = False
MAX_VISUALIZATION_READS = 200

file_path = "dataset/sample.fastq"
k = 21

print("\n--- Genome Assembly Pipeline ---")

reads = load_reads(file_path)

kmer_counts, read_count = count_kmers(reads, k)

print("Reads processed:", read_count)

G = build_graph(kmer_counts, k)

print("Graph nodes:", G.number_of_nodes())
print("Graph edges:", G.number_of_edges())

print("\n--- Graph Cleaning ---")

tips = remove_tips(G)
print("Tips removed:", tips)

bubbles = remove_bubbles(G)
print("Bubbles removed:", bubbles)

print("\n--- Contig Extraction ---")

contigs = extract_contigs(G, k)

print("Contigs generated:", len(contigs))

print("\n--- Assembly Metrics ---")

total, avg, max_len, n50, lengths = compute_metrics(contigs)

print("Total contigs:", total)
print("Average length:", round(avg,2))
print("Max length:", max_len)
print("N50:", n50)

plot_lengths(lengths)

print("\n--- Pipeline Complete ---")


# -----------------------------
# Optional Graph Visualization
# -----------------------------

if ENABLE_GRAPH_VISUALIZATION:

    import itertools
    import networkx as nx
    import matplotlib.pyplot as plt

    print("\n--- Safe Graph Visualization ---")

    small_reads = list(itertools.islice(load_reads(file_path), MAX_VISUALIZATION_READS))

    small_kmers, _ = count_kmers(small_reads, k)

    small_graph = build_graph(small_kmers, k)

    print("Visualization graph nodes:", small_graph.number_of_nodes())
    print("Visualization graph edges:", small_graph.number_of_edges())

    plt.figure(figsize=(10,8))

    pos = nx.spring_layout(small_graph)

    nx.draw(
        small_graph,
        pos,
        node_size=30,
        node_color="blue",
        edge_color="gray",
        with_labels=False
    )

    plt.title("Sample de Bruijn Graph")

    plt.savefig("output/debruijn_graph_sample.png")

    print("Graph visualization saved.")
