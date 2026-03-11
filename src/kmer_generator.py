from Bio import SeqIO

file_path = "dataset/sample.fastq"

k = 21  # typical k-mer size

reads = []

for record in SeqIO.parse(file_path, "fastq"):
    reads.append(str(record.seq))


def generate_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers.append(kmer)
    return kmers


all_kmers = []

for read in reads:
    kmers = generate_kmers(read, k)
    all_kmers.extend(kmers)


print("Total reads:", len(reads))
print("Total k-mers generated:", len(all_kmers))
print("Example k-mer:", all_kmers[0])
