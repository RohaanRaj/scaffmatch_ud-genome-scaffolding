from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt

file_path = "dataset/sample.fastq"
k = 21

reads = [str(record.seq) for record in SeqIO.parse(file_path, "fastq")]


def generate_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]


# Build graph
G = nx.DiGraph()

for read in reads:
    for kmer in generate_kmers(read, k):

        prefix = kmer[:-1]
        suffix = kmer[1:]

        if G.has_edge(prefix, suffix):
            G[prefix][suffix]["weight"] += 1
        else:
            G.add_edge(prefix, suffix, weight=1)


# Extract contigs
def extract_contigs(graph):

    contigs = []

    for node in graph.nodes():

        if graph.in_degree(node) != 1 or graph.out_degree(node) != 1:

            if graph.out_degree(node) > 0:

                for next_node in graph.successors(node):

                    contig = node
                    current = next_node

                    while graph.in_degree(current) == 1 and graph.out_degree(current) == 1:

                        contig += current[-1]
                        current = list(graph.successors(current))[0]

                    contig += current[-1]

                    contigs.append(contig)

    return contigs


contigs = extract_contigs(G)

lengths = [len(c) for c in contigs]

# Plot histogram
plt.hist(lengths, bins=30)

plt.title("Contig Length Distribution")
plt.xlabel("Contig Length (bp)")
plt.ylabel("Number of Contigs")

plt.savefig("output/contig_length_distribution.png")
