from Bio import SeqIO

file_path = "dataset/sample.fastq"

reads = []

for record in SeqIO.parse(file_path, "fastq"):
    reads.append(str(record.seq))

print("Total reads loaded:", len(reads))
print("Example read:", reads[0])
print("Read length:", len(reads[0]))
