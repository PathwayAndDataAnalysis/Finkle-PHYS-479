import sequence_search


def windowed_ranked_sequences(names, sites, p_values, window):
    original_order = {p_value : i for i, p_value in enumerate(p_values)}
    p_values.sort()
    negative_p_values, positive_p_values = [], []
    for p_value in p_values:
        if p_value < 0: negative_p_values.append(p_value)
        else: positive_p_values.append(p_value)
    sorted_p_values = positive_p_values + negative_p_values
    sequences = sequence_search.get_sequences(names)
    return [sequence_search.windowed_sequence(sequence, site, window)
            for sequence, site in zip(sequences, sites)]
    
