import analysis
import sequence_search
from scipy.stats import fisher_exact



def new_least_p_value(least_p_value, table, alternative):
    return min(least_p_value, )


def least_p_value(sequences, index, letter, total_favorable, alternative):
    least_p_value, favorable, unfavorable, length = 1, 0, 0, len(sequences)
    total_unfavorable = length - total_favorable
    best_threshold = 0; best_table = [[0,0],[0,0]]
    if alternative == "less": chosen_sequences = reversed(sequences[length//8:])
    elif alternative == "greater": chosen_sequences = sequences[:length//8]
    for threshold, sequence in enumerate(chosen_sequences):
        if sequence[index] == letter: favorable += 1
        else: unfavorable += 1
        if threshold % 4 == 0:
            table = [[total_unfavorable - unfavorable, unfavorable],
                     [total_favorable - favorable, favorable]]
            new_p_value = fisher_exact(table, alternative)[1]
            if new_p_value < least_p_value:
                least_p_value = new_p_value
                best_table = table
                best_threshold = threshold
    print(best_threshold, best_table, least_p_value,
          fisher_exact(best_table, alternative)[1])
    return least_p_value


def most_significant_p_values(sequences, index, letter, total_favorable):
    """
    Return the p-values of least absolute magnitude calculated by Fisher's
    exact test on a list of windowed sequences, ranked by phenotypic p-value from
    least to greatest positive and then least to greatest negative, by sliding a
    threshold up from greatest negative and then  down from least positive.
    """
    print(letter)
    return (least_p_value(sequences, index, letter, total_favorable, "less"),
            least_p_value(sequences, index, letter, total_favorable, "greater"))


def all_most_significant_p_values(sequences, letter_counts):
    """
    Return the the enrichments and deficiencies of each letter at each position
    in the ranked sequences.
    """
    print(len(sequences), "sequences\n")
    indices = len(sequences[0]); middle = indices // 2
    all_most_significant_p_values_list = [[] for _ in range(indices)]
    for index in range(indices):
        if index == middle: continue
        print(-(middle - index))
        for i, count in enumerate(letter_counts[index]):
            if count == 0:
                all_most_significant_p_values_list[index].append((0,0)); continue
            values = most_significant_p_values(sequences, index, chr(i), count)
            all_most_significant_p_values_list[index].append(values)
        print()
    return all_most_significant_p_values_list
    
        
   
