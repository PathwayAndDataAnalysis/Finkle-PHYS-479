import analysis
from scipy.stats import fisher_exact



def least_p_value(sequences, index, letter, total_favorable, step, alternative):
    """
    Return the least p-value calculated by Fisher's Exact Test on a 
    list of windowed sequences by sliding a threshold down it.

    The sequences are ranked by phenotypic p-value from most significant
    to most insignificant enrichment most significant to most
    insignificant deficiency.  The appearances of the letter at the index
    in the sequences are counted down to the threshold, whereupon
    Fisher's Exact test is applied to calculate the p-value of the count,
    with this p-value is compared to the least one found so far and saved
    if lesser.  Then, the threshold is advanced, and counting resumes,
    with this process repeating down the list.
        
    :param sequences: sequences of which the p value of a letter at an
                      index is desired
    :type sequences: list[str]
    :param index: position in each sequence at which the p-value of a 
                  letter is to be calculated
    :param letter: amino acid of which the p-value is to be calculated
    :type letter: str
    :param total_favorable: sum of the appearances of the letter in the
                            column
    :param step: how many more letters should be considered in each 
                 iteration of the expensive p-value calculation
    :param alternative: whether the deficiency ("lesser") or enrichment
                        ("greater") p-value should be calculated
    :return: least p-value of the letter
    :rtype: float
    """
    least_p_value, favorable, unfavorable, length = 1, 0, 0, len(sequences)
    total_unfavorable = length - total_favorable
    for threshold, sequence in enumerate(sequences):
        if sequence[index] == letter: favorable += 1
        else: unfavorable += 1
        if threshold % step == 0:
            table = [[total_unfavorable - unfavorable, unfavorable],
                     [total_favorable - favorable, favorable]]
            least_p_value = min(least_p_value,
                                fisher_exact(table, alternative)[1])
    return least_p_value


def most_significant_p_values(sequences, index, letter, total_favorable, step):
    """
    Return the least p-values calculated by Fisher's Exact Test on a 
    list of windowed sequences, ranked by phenotypic p-value from
    most significant to most insignificant enrichment most significant
    to most insignificant deficiency, by sliding a threshold down the list.
    
    :param sequences: sequences of which the p value of a letter at an index is
                      desired
    :type sequences: list[str]
    :param index: position in each sequence at which the p-value of a letter is
                  to be calculated
    :param letter: amino acid of which the p-value is to be calculated
    :type letter: str
    :param total_favorable: sum of the appearances of the letter in the column
    :param step: how many more letters should be considered in each iteration of
                 the expensive p-value calculation
    :return: least deficiency and enrichment p-values of the letter
    :rtype: tuple(float, float)
    """
    return (
        least_p_value(sequences, index, letter, total_favorable, step, "less"),
        least_p_value(sequences, index, letter, total_favorable, step, "greater")
    )


def all_most_significant_p_values(sequences, letter_counts, step):
    """
    Return the the enrichments and deficiencies of each letter at each
    position in the ranked windowed sequences.
    
    :param sequences: ranked windowed sequences of which the 
                      enrichments and deficiencies of each letter at 
                      each position are to be calculated
    :type sequences: list[str]
    :param step: how many more letters should be considered in each
                 iteration of the expensive p-value calculation
    :return: enrichments and deficiencies of each letter at each
             position in the ranked windowed sequences
    :rtype: tuple(tuple(tuple(float, float)))
    """
    indices = len(sequences[0]); middle = indices // 2
    all_most_significant_p_values_list = [[] for _ in range(indices)]
    for index in range(indices):
        if index == middle: continue # the middle is not to be counted
        for c, count in enumerate(letter_counts[index]):
            if count == 0: # ignore letters that don't appear at all
                all_most_significant_p_values_list[index].append((None, None))
                continue
            values = most_significant_p_values(
                sequences, index, chr(c), count, step
            )
            all_most_significant_p_values_list[index].append(values)
        all_most_significant_p_values_list[index] = tuple(
            all_most_significant_p_values_list[index]
        )
    del all_most_significant_p_values_list[middle]
    return tuple(all_most_significant_p_values_list)
    
        
   
