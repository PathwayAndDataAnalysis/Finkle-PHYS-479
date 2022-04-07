import sequence_search


def windowed_ranked_sequences(path, window_width):
    with open(path) as f:
        rows = [line.split("\t") for line in f.readlines()[1:]]

    p_values = [float(row[-1]) for row in rows]
    site_numbers = [int(row[2][1:]) for row in rows]

    original_order = {p_value : i for i, p_value in enumerate(p_values)}

    p_values.sort()
    negative_p_values, positive_p_values = [], []
    for p_value in p_values:
        if p_value < 0: negative_p_values.append(p_value)
        else: positive_p_values.append(p_value)

    sorted_p_values = positive_p_values + negative_p_values
    names = [rows[original_order[p_value]][1] for p_value in sorted_p_values]

    sequences = sequence_search.get_sequences(names)
    return [sequence_search.windowed_sequence(sequence, site_number, window_width)
            for sequence, site_number in zip(sequences, site_numbers)]
    


    
