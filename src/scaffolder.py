from Bio import SeqIO

file_path = "dataset/sample.fastq"
k = 21

reads = []

for record in SeqIO.parse(file_path, "fastq"):
    reads.append(str(record.seq))


def generate_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i+k]


# Build graph
import networkx as nx

G = nx.DiGraph()

for read in reads:
    for kmer in generate_kmers(read, k):

        prefix = kmer[:-1]
        suffix = kmer[1:]

        if G.has_edge(prefix, suffix):
            G[prefix][suffix]["weight"] += 1
        else:
            G.add_edge(prefix, suffix, weight=1)


# Contig extraction
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

print("Contigs:", len(contigs))


# Build scaffold graph
scaffolds = []

for contig in contigs:
    scaffolds.append(contig)

print("Scaffolds generated:", len(scaffolds))

print("Example scaffold:", scaffolds[0])
