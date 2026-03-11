import matplotlib.pyplot as plt
from collections import Counter
import itertools

from read_loader import load_reads
from kmer_counter import count_kmers


# Experimental results from your runs
reads = [10000, 30000, 50000, 75000, 90000]

n50 = [183, 188, 206, 239, 250]

avg_contig = [111, 116, 124, 139, 149]

contigs = [3007, 13272, 24191, 34248, 38125]

graph_nodes = [279232, 1287341, 2555906, 4126412, 4977674]


# -----------------------------
# Figure 4.2 - Reads vs N50
# -----------------------------

plt.figure()
plt.plot(reads, n50, marker='o')
plt.title("Reads vs N50")
plt.xlabel("Number of Reads")
plt.ylabel("N50 (bp)")
plt.grid(True)

plt.savefig("output/fig_4_1_reads_vs_n50.png", dpi=300)


# -----------------------------
# Figure 4.3 - Reads vs Average Contig Length
# -----------------------------

plt.figure()
plt.plot(reads, avg_contig, marker='o')
plt.title("Reads vs Average Contig Length")
plt.xlabel("Number of Reads")
plt.ylabel("Average Contig Length (bp)")
plt.grid(True)

plt.savefig("output/fig_4_2_reads_vs_avg_contig.png", dpi=300)


# -----------------------------
# Figure 4.4 - Reads vs Contig Count
# -----------------------------

plt.figure()
plt.plot(reads, contigs, marker='o')
plt.title("Reads vs Number of Contigs")
plt.xlabel("Number of Reads")
plt.ylabel("Contig Count")
plt.grid(True)

plt.savefig("output/fig_4_3_reads_vs_contigs.png", dpi=300)


# -----------------------------
# Figure 4.5 - Reads vs Graph Size
# -----------------------------

plt.figure()
plt.plot(reads, graph_nodes, marker='o')
plt.title("Reads vs Graph Size")
plt.xlabel("Number of Reads")
plt.ylabel("Graph Nodes")
plt.grid(True)

plt.savefig("output/fig_4_4_reads_vs_graph_nodes.png", dpi=300)


# -----------------------------
# Figure 4.6 - Assembly Quality vs Read Coverage
# -----------------------------

plt.figure()

fig, ax1 = plt.subplots()

ax1.plot(reads, n50, marker='o', label="N50", color="blue")
ax1.plot(reads, avg_contig, marker='o', label="Avg Contig Length", color="orange")

ax1.set_xlabel("Number of Reads")
ax1.set_ylabel("Length (bp)")

ax2 = ax1.twinx()

ax2.plot(reads, contigs, marker='o', label="Contig Count", color="green")

ax2.set_ylabel("Contig Count")

plt.title("Assembly Quality Summary")

ax1.legend(loc="upper left")
ax2.legend(loc="upper right")

plt.grid(True)

plt.savefig("output/fig_4_5_assembly_summary.png", dpi=300)


# -----------------------------
# Figure 4.7 - k-mer Spectrum
# -----------------------------

print("\nGenerating k-mer spectrum plot...")

file_path = "dataset/sample.fastq"
k = 21
MAX_READS = 2000

reads_subset = list(itertools.islice(load_reads(file_path), MAX_READS))

kmer_counts, _ = count_kmers(reads_subset, k)

frequency_counts = Counter(kmer_counts.values())

x = list(frequency_counts.keys())
y = list(frequency_counts.values())

plt.figure()

plt.scatter(x, y)

plt.title("k-mer Frequency Spectrum")
plt.xlabel("k-mer Frequency")
plt.ylabel("Number of k-mers")

plt.savefig("output/fig_4_6_kmer_spectrum.png", dpi=300)

print("k-mer spectrum saved.")


print("\nAll experiment plots saved in output/")
