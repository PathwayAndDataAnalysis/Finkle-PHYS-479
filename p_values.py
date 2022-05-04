import analysis
import sequence_search
from scipy.stats import fisher_exact



def least_p_value(sequences, index, letter, total_favorable, alternative, step):
    least_p_value, favorable, unfavorable, length = 1, 0, 0, len(sequences)
    total_unfavorable = length - total_favorable
    for threshold, sequence in enumerate(sequences):
        if sequence[index] == letter: favorable += 1
        else: unfavorable += 1
        if threshold % step == 0:
            table = [[total_unfavorable - unfavorable, unfavorable],
                     [total_favorable - favorable, favorable]]
            least_p_value = min(least_p_value, fisher_exact(table, alternative)[1])
    return least_p_value


def most_significant_p_values(sequences, index, letter, total_favorable, step):
    """
    Return the p-values of least absolute magnitude calculated by Fisher's
    exact test on a list of windowed sequences, ranked by phenotypic p-value from
    least to greatest positive and then least to greatest negative, by sliding a
    threshold up from greatest negative and then  down from least positive.
    """
    return (least_p_value(sequences, index, letter, total_favorable, "less", step),
            least_p_value(sequences, index, letter, total_favorable, "greater", step))


def all_most_significant_p_values(sequences, letter_counts, step):
    """
    Return the the enrichments and deficiencies of each letter at each position
    in the ranked sequences.
    """
    indices = len(sequences[0]); middle = indices // 2
    all_most_significant_p_values_list = [[] for _ in range(indices)]
    for index in range(indices):
        if index == middle: continue
        print(-(middle - index))
        for i, count in enumerate(letter_counts[index]):
            if count == 0:
                all_most_significant_p_values_list[index].append((0,0)); continue
            values = most_significant_p_values(sequences, index, chr(i), count, step)
            all_most_significant_p_values_list[index].append(values)
            print(chr(i), "|", values[0], "|", values[1])
    return all_most_significant_p_values_list
    
        
   
