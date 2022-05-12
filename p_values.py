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
    :type index: int
    :param letter: amino acid of which the p-value is to be calculated
    :type letter: str
    :param total_favorable: sum of the appearances of the letter in the
                            column
    :param step: how many more letters should be considered in each 
                 iteration of the expensive p-value calculation
    :type step: int
    :param alternative: whether the deficiency ("lesser") or enrichment
                        ("greater") p-value should be calculated
    :type alternative: str
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

    Return (0, 0) if the total_favorable is zero.
    
    :param sequences: sequences of which the p value of a letter at an
                      index is desired
    :type sequences: list[str]
    :param index: position in each sequence at which the p-value of a
                  letter is to be calculated
    :param index: int
    :param letter: amino acid of which the p-value is to be calculated
    :type letter: str
    :param total_favorable: sum of the appearances of the letter in the
                            column
    :type total_favorable: int
    :param step: how many more letters should be considered in each
                 iteration of the expensive p-value calculation
    :type step: int
    :return: least deficiency and enrichment p-values of the letter
    :rtype: tuple(float, float)
    """
    if total_favorable == 0: return 0, 0
    return (
        least_p_value(sequences, index, letter, total_favorable, step, "less"),
        least_p_value(sequences, index, letter, total_favorable, step, "greater")
    )


def index_p_values(sequences, letter_counts, index, step):
    """
    Return the deficiency and enrichment p values of all letters at an
    index in tuple of aligned, windowed sequences.

    :param letter_counts: count of every letter at every index of the
                          sequences
    :type letter_counts: tuple(tuple(int))
    :index: index in question
    :type index: int
    :param step: how many more letters should be considered in each
                 iteration of the expensive p-value calculation
    :type step: int
    :return: the deficiency and enrichment p values of all letters at an
             index in tuple of aligned, windowed sequences
    :rtype: tuple(tuple(int, int))
    """
    return tuple(most_significant_p_values(sequences, index, chr(c), count, step)
                 for c, count in enumerate(letter_counts[index]))


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
    return tuple(
        index_p_values(sequences, letter_counts, index, step)
        for index in range(indices) if index != middle
    )
    
        
   
