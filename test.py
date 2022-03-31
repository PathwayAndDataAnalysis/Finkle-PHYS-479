#import p_values
import sequence_search
import scipy


def hypergeometric_test():
    """
    Assume the letters a, b, c, and d represent the counts as displayed below.

      - +
    - a b
    + c d

    Then if a = 4, b = 6, c = 6, d = 4,
    The enrichment p-value should be 0.9105522960012121
    The deficiency p-value should be 0.3281408993483296

    If a = 20, b = 30, c = 25, d = 30,
    The enrichment p-value should be 0.7766662300662146
    The deficiency p-value should be 0.357198690677301

    If a = 15, b = 8, c = 20, d = 42
    The enrichment p-value should be 0.006445865568610187
    The deficiency p-value should be 0.9985899806396821

    If a = 21, b = 20, c = 34, d = 13
    The enrichment p-value should be 0.9883420938210076
    The deficiency p-value should be 0.0341612031176084
    """
    tables = [[4, 6, 6, 4],
              [20, 30, 25, 20],
              [15, 8, 20, 42],
              [21, 20, 34, 13]]
    accepted_p_values_list = [[0.9105522960012121, 0.3281408993483296],
                              [0.7766662300662146, 0.357198690677301],
                              [0.006445865568610187, 0.9985899806396821],
                              [0.9883420938210076, 0.0341612031176084]]
    failures = 0
    print()
    for table, accepted_p_values in zip(tables, accepted_p_values_list):
        deficiency_p_value = p_values.deficiency_p_value(*table)
        enrichment_p_value = p_values.enrichment_p_value(*table)
        print(*table[:2])
        print(*table[2:])
        print("           Enrichment      |     Deficiency")
        print("Accepted", *accepted_p_values)
        print("Found   ", enrichment_p_value, "", deficiency_p_value)
        fisher_table = [[table[0], table[1]],
                        [table[2], table[3]]]
        right = scipy.stats.fisher_exact(fisher_table, 'less')[1]
        left = scipy.stats.fisher_exact(fisher_table, 'greater')[1]
        print(left, right)

def windowed_sequence_test():
    gene_names = ["UL38", "CVC1", "U69", "nef", "RTN4"]
    site_numbers = [1, 2, 3, 4, 107]
    window_sizes = [3, 4, 5, 6, 5]
    windowed_sequences = ["T T H S T A A",
                          "Q E V L S N E E A",
                          "L K K Q I S A C S D M",
                          "D G V G A A S R D L E K H",
                          "P E R Q P S W D P S P"]
    print()
    for name, number, size, sequence in zip(gene_names, site_numbers,
                                            window_sizes, windowed_sequences):
        print("Gene:    ", name)
        print("Accepted:", sequence)
        print("Found:   ", *sequence_search.windowed_sequence(name, number, size))
        print()

print("p-value")        
#hypergeometric_test()
print()
print("Windowed sequence")
windowed_sequence_test()


        
                           
