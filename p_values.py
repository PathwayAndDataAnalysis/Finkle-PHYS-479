import analysis
import sequence_search
import scipy


def sequence_p_values(selected_sequences, remaining_sequences):
    """
    Return the p-values of the count of every amino acid letter in each columnar
    position of aligned rows of amino acid sequences selected from a corpus.
    """
    letters = ["A", "G", "I", "L", "P", "V", "F", "W", "Y", "D", 
               "E", "R", "H", "K", "S", "T", "C", "M", "N", "Q"]
    enrichment_p_values = [[] for _ in columns]
    deficiency_p_values = [[] for _ in columns]
    remaining_rows = len(remaining_sequences[0])
    selected_rows = len(selected_sequences[0])
    selected_letter_counts = analysis.letter_counts(selected_sequences)
    remaining_letter_counts = analysis.letter_counts(remaining_sequences)
    for column in selected_sequences:
        for letter in letters:
            table = [[remaining_rows - remaining_letter_counts[column],
                     remaining_letter_counts[column]],
                     [selected_rows - selected_letter_counts[column],
                     selected_letter_counts[column]]
            enrichment_p_values[column].append(scipy.stats.fisher_exact(table, 'less')[1])
            deficiency_p_values[column].append(scipy.stats.fisher_exact(table, 'greater')[1])
    return enrichment_p_values, deficiency_p_values
