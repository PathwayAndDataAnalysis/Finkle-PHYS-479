import analysis
import sequence_search
from scipy.stats import fisher_exact


def least_p_value(sequences, index, letter, total_favorable, alternative):
    least_p_value, favorable, unfavorable, length = 1, 0, 0, len(sequences)
    total_unfavorable = length - total_favorable
    if alternative == "less": chosen_sequences = reversed(sequences[length//2:])
    elif alternative == "greater": chosen_sequences = sequences[:length//2]
    for i, sequence in enumerate(chosen_sequences):
        if sequence[index] == letter: favorable += 1
        else: unfavorable += 1
        if i % 4 == 0:
            table = [[total_unfavorable - unfavorable, unfavorable],
                     [total_favorable - favorable, favorable]]
            least_p_value = min(least_p_value, fisher_exact(table, alternative)[1])
    return least_p_value


def most_significant_p_values(sequences, index, letter, total_favorable):
    """
    Return the p-values of least absolute magnitude calculated by Fisher's
    exact test on a list of windowed sequences, ranked by phenotypic p-value from
    least to greatest positive and then least to greatest negative, by sliding a
    threshold up from greatest negative and then  down from least positive.
    """
    return (least_p_value(sequences, index, letter, total_favorable, "less"),
            least_p_value(sequences, index, letter, total_favorable, "greater"))


def all_most_significant_p_values(sequences, letter_counts):
    """
    Return the the enrichments and deficiencies of each letter at each position
    in the ranked sequences.
    """
    columns = len(sequences[0]); middle = columns // 2
    all_most_significant_p_values_list = [[] for _ in range(columns)]
    for column in range(columns):
        if column == middle: continue
        print(-(middle - column))
        for i, count in enumerate(letter_counts[column]):
            if count == 0:
                all_most_significant_p_values_list[column].append((0,0)); continue
            values = most_significant_p_values(sequences, column, chr(i), count)
            print(chr(i), f'{values[0]:.3f}', f'{values[1]:.3f}')
            all_most_significant_p_values_list[column].append(values)
        print()
    return all_most_significant_p_values_list
    
        
   
