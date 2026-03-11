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

        G.add_edge(prefix, suffix)


print("Graph nodes:", G.number_of_nodes())
print("Graph edges:", G.number_of_edges())

example_edge = list(G.edges())[0]
print("Example edge:", example_edge)
