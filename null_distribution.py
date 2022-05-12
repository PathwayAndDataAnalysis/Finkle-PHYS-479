import analysis
import p_values

import random

def null_distribution(sequences, window, repetitions, step, seed):            
    """
    Return the null distribution of letters in the columnns of given
    windowed sequences.

    Deficiency and enrichment p-values are not uniformly distributed
    because of the threshold searching step and therefore not p-values
    but scores. Estimating the nominal p-values (statistical
    significance) of these scores requires their null distribution.
    
    This method generates a null distribution of the deficiency and
    enrichment scores of the given sequences by randomizing their
    ranking of the, determining the enrichment and deficiency scores
    for each amino acid at each index, collecting the scores, and
    repeating this process to collect enough random scores for their
    distribution to be statistically significant.

    :param sequences: aligned, windowed amino acid sequences, the
                      letters of which are to be filtered
    :type sequences: list[str]
    :param window: width of the view to either side of the index
    :type window: int
    :param repetitions: how many times to shuffle the sequences and
                        record their new scores.  A number between
                        1,000 and 1,000,000 is ideal.
    :param step: how many more letters should be considered in each 
                 iteration of the expensive p-value calculation
    :type step: int
    :param seed: seed for the random shuffling of sequences
    :type seed: float
    :return: the deficiency and enrichment scores and for each column,
             for each letter
    """
    length = 2 * window + 1
    columns = [[sequence[i] for sequence in sequences] for i in range(length)]
    letter_counts = analysis.letter_counts(columns)
    random.seed(seed)
    results = []
    for _ in range(repetitions):
        random.shuffle(sequences)
        results.append(p_values.all_most_significant_p_values(
            sequences, letter_counts, step
            )
        )
    return [[result[i] for result in results] for i in range(length - 1)]
            
