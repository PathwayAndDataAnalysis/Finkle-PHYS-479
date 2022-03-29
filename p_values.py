import analysis
import sequence_search
from scipy.stats import hypergeom
    

def deficiency_p_value(remaining_unfavorable,
                       selected_unfavorable,
                       remaining_favorable,
                       selected_favorable):
    """
    a b
    c d

      3     10
    a - 1 b + 1
    c + 1 d - 1

    15 8
    20 

    
    """
    unfavorable = remaining_unfavorable + selected_unfavorable
    favorable = remaining_favorable + selected_favorable
    selected = selected_unfavorable + selected_favorable
    trials = unfavorable + favorable
    least_favorable = max(0, selected + favorable - trials)
    less_favorable = [i for i in range(least_favorable, selected_favorable + 1)]
    return sum(hypergeom.pmf(less_favorable, trials, selected_favorable, selected))


def enrichment_p_value(remaining_unfavorable,
                       selected_unfavorable,
                       remaining_favorable,
                       selected_favorable):
    unfavorable = remaining_unfavorable + selected_unfavorable
    favorable = remaining_favorable + selected_favorable
    selected = selected_unfavorable + selected_favorable
    trials = unfavorable + favorable
    most_favorable = min(favorable, selected)
    more_favorable = [i for i in range(selected_favorable, most_favorable + 1)]
    return sum(hypergeom.pmf(more_favorable, trials, favorable, selected))


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
            table = (remaining_rows - remaining_letter_counts[column],
                     remaining_letter_counts[column],
                     selected_rows - selected_letter_counts[column],
                     selected_letter_counts[column])
            enrichment_p_values[column].append(enrichment_p_value(*table))
            deficiency_p_values[column].append(deficiency_p_value(*table))
    return enrichment_p_values, deficiency_p_values
