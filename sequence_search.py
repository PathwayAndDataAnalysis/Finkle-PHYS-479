import random



def get_sequences(gene_names):
    """
    Return the sequences of the named genes in the database and the
    index of any sequence that could not be returned.
    
    :param gene_names: identifiers of the genes
    :type gene_names: list[str]
    :return: sequences of the named genes
    :rtype: list[str]
    """
    with open('data/gene_names.txt') as f:
        lines = f.readlines()
    gene_index_dictionary = {line[:-1] : i for i, line in enumerate(lines)}
    with open ('data/raw_sequences.txt') as f:
        sequences = f.readlines()
    results, failures = [], []
    for n, name in enumerate(gene_names):
        try:
            result = sequences[gene_index_dictionary[name]]
            results.append(result if result[-1] != "\n" else result[:-1])
        except:
            failures.append(n)
    return results, failures


def windowed_sequence(sequence, index, window):
    """
    Return a view of the letters one window width left and right of an 
    index, including the letter there.
    
    :param sequence: sequence of which a windowed view is sought
    :type sequence: str
    :param index: index to either side of which the window spans
    :type site_index: int
    :param window: width of the view to either side of the index
    :type window: int
    :return: windowed view of the letters in the sequence
    :rtype str
    """
    return sequence[index - window - 2: index + window - 1]


def ranked_windowed_sequences(sequences, indices, data_p_values, window):
    """
    Return a list of windowed sequences ranked by p-value.
    
    Acceptable p-values range from [-1, 1], with negative p-values
    indicating deficiency and positive ones enrichment of an associated
    protein in a cell.
    
    The ranking is 0 -> 1, -1 -> 0.
    
    :param names: amino acid sequences
    :type names: list[str]
    :param indices: positions at which to view each sequence
    :type indices: list[int]
    :param data_p_values: p values of the sequences
    :type data_p_values: list[float]
    :param window: width of the view to either side of the index
    :type window: int
    :return: list of windowed sequences ranked by p value
    :rtype: list[str]
    
    """
    original_order = {p_value : i for i, p_value in enumerate(data_p_values)}
    data_p_values.sort()
    deficient_sequences, enriched_sequences = [], []
    for p_value in data_p_values:
        if p_value < 0:
            deficient_sequences.append(sequences[original_order[p_value]])
        else:
            enriched_sequences.append(sequences[original_order[p_value]])
    return enriched_sequences + deficient_sequences
    
