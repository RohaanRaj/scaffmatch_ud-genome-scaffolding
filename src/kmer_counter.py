from collections import defaultdict
from tqdm import tqdm

DNA_MAP = {'A':0,'C':1,'G':2,'T':3}

def encode_kmer(seq):
    val = 0
    for base in seq:
        val = (val << 2) | DNA_MAP[base]
    return val


def count_kmers(reads, k):

    kmer_counts = defaultdict(int)
    read_count = 0

    for seq in tqdm(reads, desc="Processing reads"):

        read_count += 1

        for i in range(len(seq) - k + 1):
            kmer = encode_kmer(seq[i:i+k])
            kmer_counts[kmer] += 1

    return kmer_counts, read_count
