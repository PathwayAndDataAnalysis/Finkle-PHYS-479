#import analysis
#import sequence_search
from scipy.stats import hypergeom
    

def p_value(remaining_unfavorable,
            selected_unfavorable,
            remaining_favorable,
            selected_favorable):
    
    unfavorable = remaining_unfavorable + selected_unfavorable
    favorable = remaining_favorable + selected_favorable
    selected = selected_unfavorable + selected_favorable
    trials = unfavorable + favorable
    limit = favorable + max(remaining_unfavorable - favorable, 0)
    possibilities = [i for i in range(favorable, limit)]
    return sum(hypergeom.pmf(possibilities, trials, favorable, selected))


def p_values(selected_sequences, remaining_sequences):
    """
    Return the p-values of the count of every amino acid letter in each columnar
    position of aligned rows of amino acid sequences selected from a corpus.
    """
    letters = ["A", "G", "I", "L", "P", "V", "F", "W", "Y", "D", 
               "E", "R", "H", "K", "S", "T", "C", "M", "N", "Q"]
    p_values = [[] for _ in columns]
    remaining_rows = len(remaining_sequences[0])
    selected_rows = len(selected_sequences[0])
    selected_letter_counts = analysis.letter_counts(selected_sequences)
    remaining_letter_counts = analysis.letter_counts(remaining_sequences)
    for column in selected_sequences:
        for letter in letters:
            p_values[column].append(p_value(
                remaining_rows - remaining_letter_counts[column],
                remaining_letter_counts[column],
                selected_rows - selected_letter_counts[column],
                selected_letter_counts[column])
                )
    return p_values
