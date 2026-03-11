from Bio import SeqIO

def load_reads(file_path):
    for record in SeqIO.parse(file_path, "fastq"):
        seq = str(record.seq)

        if "N" in seq:
            continue

        yield seq
