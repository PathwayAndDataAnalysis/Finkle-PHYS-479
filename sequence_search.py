def get_sequences(gene_names):
    """Return the sequences of the named genes in the database."""
    with open('data/gene_names.txt') as f:
        lines = f.readlines()
    gene_index_dictionary = {line[:-1] : i for i, line in enumerate(lines)}
    with open ('data/raw_sequences.txt') as f:
        sequences = f.readlines()
    results = []
    for i, name in enumerate(gene_names):
        try:
            results.append(sequences[gene_index_dictionary[name]])
        except:
            continue
    return results


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
