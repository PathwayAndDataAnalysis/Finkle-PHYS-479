def get_sequences(gene_names):
    """Return the sequences of the named genes in the database."""
    with open('data/gene_names.txt') as f:
        lines = f.readlines()
    gene_index_dictionary = {line[:-1] : i for i, line in enumerate(lines)}
    with open ('data/raw_sequences.txt') as f:
        sequences = f.readlines()
    results = []
    for name in gene_names:
        try:
            results.append(sequences[gene_index_dictionary[name]])
        except:
            continue
    return results


def windowed_sequence(sequence, site_index, window_width):
    """
    Return a view of sequences aligned by phosphorylation site and limited by
    a window width to either site of it.
    """
    left = site_index - window_width; right = site_index + window_width + 1
    if left < 0 or right > len(sequence): print("ding")
    return sequence[left:right]
