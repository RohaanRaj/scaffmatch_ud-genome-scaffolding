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


def extract_contigs(graph):

    contigs = []
    visited = set()

    for node in graph.nodes():

        if graph.in_degree(node) != 1 or graph.out_degree(node) != 1:

            if graph.out_degree(node) > 0:

                for next_node in graph.successors(node):

                    contig = node
                    current = next_node

                    while graph.in_degree(current) == 1 and graph.out_degree(current) == 1:

                        contig += current[-1]

                        visited.add(current)

                        current = list(graph.successors(current))[0]

                    contig += current[-1]

                    contigs.append(contig)

    return contigs


contigs = extract_contigs(G)

print("Total contigs generated:", len(contigs))
print("Example contig:", contigs[0])
print("Example contig length:", len(contigs[0]))
