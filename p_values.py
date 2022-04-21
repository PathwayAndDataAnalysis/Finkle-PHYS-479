import analysis
import sequence_search
from scipy.stats import fisher_exact


def most_significant_p_values(sequences, index, letter, total_favorable):
    """
    Return the p-values of least absolute magnitude calculated by Fisher's
    exact test on a list of windowed sequences, ranked by phenotypic p-value from
    least to greatest positive and then least to greatest negative, by sliding a
    threshold up from greatest negative and then  down from least positive.
    """
    least_deficiency_p_value, least_enrichment_p_value = 1, 1
    length = len(sequences)
    total_unfavorable = length - total_favorable
    
    favorable, unfavorable = 0, 0
    for i, sequence in enumerate(reversed(sequences[length//2:])):
        if sequence[index] == letter: favorable += 1
        else: unfavorable += 1
        if not i % 4 == 0: continue
        table = [[total_unfavorable - unfavorable, unfavorable],
                 [total_favorable - favorable, favorable]]
        least_deficiency_p_value = min(least_deficiency_p_value,
                                       fisher_exact(table, 'less')[1])
        
    favorable, unfavorable = 0, 0
    for i, sequence in enumerate(sequences[:length//2]):
        if sequence[index] == letter: favorable += 1
        else: unfavorable += 1
        if not i % 4 == 0: continue
        table = [[total_unfavorable - unfavorable, unfavorable],
                 [total_favorable - favorable, favorable]]
        least_enrichment_p_value = min(least_enrichment_p_value,
                                       fisher_exact(table, 'greater')[1])
        
    return least_deficiency_p_value, least_enrichment_p_value


def all_most_significant_p_values(sequences, letter_counts):
    """
    Return the the enrichments and deficiencies of each letter at each position
    in the ranked sequences.
    """
    columns = len(sequences[0])
    all_most_significant_p_values_list = [[] for _ in range(columns)]
    for column in range(columns):
        for i, count in enumerate(letter_counts[column]):
            if count == 0:
                all_most_significant_p_values_list[column].append((0,0)); continue
            values = most_significant_p_values(sequences, column, chr(i), count)
            print(chr(i), values)
            all_most_significant_p_values_list[column].append(values)
    return all_most_significant_p_values_list
    
        
   
