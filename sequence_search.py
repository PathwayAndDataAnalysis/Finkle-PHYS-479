def get_sequence(gene_name):
    """Return the index of the named gene in the database."""
    with open('data/gene_names.txt') as f:
        gene_names = [line[:-1] for line in f.readlines()]
        sequence_index = gene_names.index(gene_name)
    with open ('data/raw_sequences.txt') as sequences:
        for i, sequence in enumerate(sequences):
            if i == sequence_index:
                return sequence


def windowed_sequence(gene_name, site_number, window_width):
    """
    Return a view of sequences aligned by phosphorylation site and limited by
    a window width to either site of it.
    """
    site_index, site_count = 0, 0
    sequence = get_sequence(gene_name)
    windowed_sequence = []
    for i, letter in enumerate(sequence):
        if letter.upper() == "S": site_count += 1
        if site_count == site_number: site_index = i; break
    for i in range(site_index - window_width, site_index + window_width + 1):
        windowed_sequence.append(sequence[i])
    return windowed_sequence
