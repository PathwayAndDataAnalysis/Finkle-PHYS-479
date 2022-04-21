import analysis
import sequence_search
from scipy.stats import fisher_exact


def sequence_p_values(selected_sequences, remaining_sequences):
    """
    Return the p-values of the count of every amino acid letter in each columnar
    position of aligned rows of amino acid sequences selected from a corpus.
    """
    letters = ["A", "G", "I", "L", "P", "V", "F", "W", "Y", "D", 
               "E", "R", "H", "K", "S", "T", "C", "M", "N", "Q"]
    columns = len(selected_sequences)
    enrichment_p_values = [[] for _ in range(columns)]
    deficiency_p_values = [[] for _ in range(columns)]
    remaining_rows = len(remaining_sequences[0])
    selected_rows = len(selected_sequences[0])
    selected_letter_counts = analysis.letter_counts(selected_sequences)
    remaining_letter_counts = analysis.letter_counts(remaining_sequences)
    for column in selected_sequences:
        for letter in letters:
            table = [[remaining_rows - remaining_letter_counts[column],
                     remaining_letter_counts[column]],
                     [selected_rows - selected_letter_counts[column],
                     selected_letter_counts[column]]]
            enrichment_p_values[column].append(fisher_exact(table, 'less')[1])
            deficiency_p_values[column].append(fisher_exact(table, 'greater')[1])
    return enrichment_p_values, deficiency_p_values


def most_significant_p_values(sequences, column_index, letter, favorable):
    """
    Return the p-values of least absolute magnitude calculated by Fisher's
    exact test on a list of windowed sequences, ranked by phenotypic p-value from
    least to greatest positive and then least to greatest negative, by sliding a
    threshold up from greatest negative and then  down from least positive.
    """
    least_deficiency_p_value, least_enrichment_p_value = 1, 1
    
    selected_favorable, selected_unfavorable = 0, 0
    remaining_favorable = favorable
    remaining_unfavorable = len(sequences) - favorable
    for sequence in reversed(sequences):
        if sequence[column_index] == letter:
            remaining_favorable -= 1; selected_favorable += 1;
        else:
            remaining_unfavorable -= 1; selected_unfavorable += 1
        table = [[selected_unfavorable, selected_favorable],
                 [remaining_unfavorable, remaining_favorable]]
        deficiency_p_value = fisher_exact(table, alternative = 'less')[1]
        least_deficiency_p_value = min(least_deficiency_p_value,
                                       deficiency_p_value)
        
    selected_favorable, selected_unfavorable = 0, 0
    remaining_favorable = favorable
    remaining_unfavorable = len(sequences) - favorable
    for sequence in sequences:
        if sequence[column_index] == letter:
            remaining_favorable -= 1; selected_favorable += 1;
        else:
            remaining_unfavorable -= 1; selected_unfavorable += 1
        table = [[selected_unfavorable, selected_favorable],
                 [remaining_unfavorable, remaining_favorable]]
        enrichment_p_value = fisher_exact(table, alternative = 'greater')[1]
        least_enrichment_p_value = min(least_enrichment_p_value,
                                       enrichment_p_value)
        
    return least_deficiency_p_value, least_enrichment_p_value


def all_most_significant_p_values(sequences, letter_counts):
    """
    Return the the enrichments and deficiencies of each letter at each position
    in the ranked sequences.
    """
    columns = len(sequences[0])
    all_most_significant_p_values_list = [[] for _ in range(columns)]
    for column in range(columns):
        for i, count in enumerate(letter_counts[column]["A":]):
            all_most_significant_p_values_list[column].append(
                [most_significant_p_values(sequences, column, chr(i), count)]
            )
    return all_most_significant_p_values_list
    
        
   
