from tqdm import tqdm

def decode_kmer(val, length):

    bases = []

    for _ in range(length):
        bases.append("ACGT"[val & 3])
        val >>= 2

    return "".join(reversed(bases))


def extract_contigs(G, k):

    contigs = []

    for node in tqdm(G.nodes(), desc="Traversing graph"):

        if G.in_degree(node) != 1 or G.out_degree(node) != 1:

            if G.out_degree(node) > 0:

                for nxt in G.successors(node):

                    contig = decode_kmer(node, k-1)
                    current = nxt

                    while G.in_degree(current) == 1 and G.out_degree(current) == 1:

                        contig += decode_kmer(current,1)
                        current = list(G.successors(current))[0]

                    contig += decode_kmer(current,1)

                    contigs.append(contig)

    return contigs
