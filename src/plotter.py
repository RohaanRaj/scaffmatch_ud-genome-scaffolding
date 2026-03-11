import matplotlib.pyplot as plt

def plot_lengths(lengths):

    plt.hist(lengths, bins=30)

    plt.title("Contig Length Distribution")
    plt.xlabel("Contig Length (bp)")
    plt.ylabel("Number of Contigs")

    plt.savefig("output/contig_length_distribution.png")
