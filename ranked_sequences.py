import sequence_search


def ranked_sequences()
    with open("sample-data.txt") as f:
        lines = f.readlines()[1:]

    p_values = []
    for line in lines:
        p_value = []
        for character in reversed(line):
            if character == " ": break
            p_value.append(character)
        p_value = float("".join(reversed(p_value)))
        p_values.append(p_value)

    original_order = {original_order[p_value] : i for i, p_value in enumerate(p_values)}

    p_values.sort()
    negative_p_values, positive_p_values = [], []
    for p_value in p_values:
        if p_values < 0: negative_p_values.append(p_value)
        else: positive_p_values.append(p_value)

    negative_p_values = reversed(negative_p_values)

    ranked_list = []
    for p_value in positive_p_values: ranked_list.append(lines[original_order[p_value]])
    for p_value in negative_p_values: ranked_list.append(lines[original_order[p_value]])

    names = []
    for line in ranked_list:
        white_space_begun = False
        white_space_ended = False
        name = []
        for character in line:
            if not white_space_begun:
                if character == " ": white_space_begun = True
            else:
                if character != " ":
                    white_space_ended = True
                    name.append(character)
                else: break
        names.append(name)

    return [sequence_search.get_sequence(name) for name in names]        
