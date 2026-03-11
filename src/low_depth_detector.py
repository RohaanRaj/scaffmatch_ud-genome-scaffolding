from Bio import SeqIO
import networkx as nx

file_path = "dataset/sample.fastq"
k = 21

reads = []

for record in SeqIO.parse(file_path, "fastq"):
    reads.append(str(record.seq))


def generate_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]


G = nx.DiGraph()

for read in reads:
    for kmer in generate_kmers(read, k):

        prefix = kmer[:-1]
        suffix = kmer[1:]

        if G.has_edge(prefix, suffix):
            G[prefix][suffix]["weight"] += 1
        else:
            G.add_edge(prefix, suffix, weight=1)


low_depth_edges = []

for u, v, data in G.edges(data=True):

    if data["weight"] <= 2:
        low_depth_edges.append((u, v))


print("Total edges:", G.number_of_edges())
print("Low-depth edges:", len(low_depth_edges))

print("Example low-depth edge:", low_depth_edges[0])
