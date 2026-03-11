def compute_metrics(contigs):

    lengths = [len(c) for c in contigs]

    total_contigs = len(lengths)
    avg_length = sum(lengths) / total_contigs
    max_length = max(lengths)

    sorted_lengths = sorted(lengths, reverse=True)

    total_length = sum(sorted_lengths)
    half_length = total_length / 2

    running_sum = 0

    for length in sorted_lengths:

        running_sum += length

        if running_sum >= half_length:
            n50 = length
            break

    return total_contigs, avg_length, max_length, n50, lengths
