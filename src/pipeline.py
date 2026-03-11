from read_loader import load_reads
from kmer_counter import count_kmers
from graph_builder import build_graph
from graph_cleaner import remove_tips, remove_bubbles
from contig_extractor import extract_contigs
from metrics import compute_metrics
from plotter import plot_lengths

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
