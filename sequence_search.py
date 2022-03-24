def windowed_sequence(gene_name, site_number, window_width):
    """
    Return a view of sequences aligned by phosphorylation site and limited by
    a window width to either site of it.
    """
    site_index, site_count = 0, 0
    sequence = sequence_search.get_sequence(gene_name)
    windowed_sequence = []
    for i, letter in enumerate(sequence):
        if letter.upper() == "S": site_count += 1
        if site_count == site_number: site_index = i; break
    for i in range(site_index - window_width, site + window_width):
        windowed_sequence.append(sequence[i])
    return windowed_sequence
