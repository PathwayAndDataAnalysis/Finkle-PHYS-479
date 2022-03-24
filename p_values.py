import analysis
import sequence_search


def p_values(selected_sequences, remaining_sequences):
    """
    Return the p-values of the count of every amino acid letter in each columnar
    position of aligned rows of amino acid sequences selected from a corpus.
    """
    letters = ["A", "G", "I", "L", "P", "V", "F", "W", "Y", "D", 
               "E", "R", "H", "K", "S", "T", "C", "M", "N", "Q"]
    columns = len(selected_sequences)
    p_values = [[] for _ in columns]
    selected_letter_counts = analysis.letter_counts(selected_sequences)
    remaining_letter_counts = analysis.letter_counts(remaining_sequences)
    for i in range(columns)
        for letter in letters:
            p_values[i].append((selected_sequences[i],
                                remaining_sequences[i],
                                selected_column_counts[i],
                                remaining_column_count[i]))
    return p_values
