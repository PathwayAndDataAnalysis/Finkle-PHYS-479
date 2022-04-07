import sequence_search
import time


def ranked_sequences(path):
    with open(path) as f:
        lines = f.readlines()[1:]
        
    p_values = [float("".join(line.split("\t")[-1])) for line in lines]

    original_order = {p_value : i for i, p_value in enumerate(p_values)}

    p_values.sort()
    negative_p_values, positive_p_values = [], []
    for p_value in p_values:
        if p_value < 0: negative_p_values.append(p_value)
        else: positive_p_values.append(p_value)

    sorted_p_values = positive_p_values + negative_p_values
    ranked_list = [lines[original_order[p_value]] for p_value in sorted_p_values]

    names = [line.split("\t")[1] for line in ranked_list]

    return sequence_search.get_sequences(names)
