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
            enrichment_p_values[column].append(scipy.stats.fisher_exact(table, 'less')[1])
            deficiency_p_values[column].append(scipy.stats.fisher_exact(table, 'greater')[1])
    return enrichment_p_values, deficiency_p_values


def most_significant_p_values(sequences, column, letter, favorable):
    """
    Return the p-values of least absolute magnitude calculated by Fisher's
    exact test on a list of windowed sequences ranked by phenotypic p-value from
    most least positive to greatest positive to least negative to greatest
    negative, first by sliding a threshold down from greatest positive and then
    by sliding up from least negative.
    """
    greatest_negative_p_value, least_positive_p_value = -1, 0
    selected_favorable, selected_unfavorable = 0, 0
    remaining_favorable = favorable
    remaining_unfavorable = len(sequences) - favorable
    for selected, sequence in enumerate(reversed(sequences)):
        is_favorable = sequence[column] == letter
        is_unfavorable = sequence[column] != letter
        remaining_unfavorable -= is_unfavorable
        remaining_favorable -= is_favorable
        selected_unfavorable += is_unfavorable
        selected_favorable += is_favorable
        table = [[selected_unfavorable, selected_favorable],
                 [remaining_unfavorable, remaining_favorable]]
        deficiency_p_value = scipy.fisher_exact(table, alternative = 'less')
        greatest_negative_p_value = max(greatest_negative_p_value,
                                        deficiency_p_value)
        
    selected_favorable, selected_unfavorable = 0, 0
    remaining_favorable = favorable
    remaining_unfavorable = len(sequences) - favorable
    for selected, sequence in enumerate(sequences):
        is_favorable = sequence[column] == letter
        is_unfavorable = sequence[column] != letter
        remaining_unfavorable -= is_unfavorable
        remaining_favorable -= is_favorable
        selected_unfavorable += is_unfavorable
        selected_favorable += is_favorable
        table = [[selected_unfavorable, selected_favorable],
                 [remaining_unfavorable, remaining_favorable]]
        enrichment_p_value = scipy.fisher_exact(table, alternative = 'greater')
        least_positive_p_value = min(least_positive_p_value,
                                     enrichment_p_value)
        
    return greatest_negative_p_value, least_positive_p_value


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
    
        
   
