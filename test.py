import p_values
import sequence_search
import windowed_ranked_sequences
import processing
import analysis

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

    If a = 15, b = 8, c = 20, d = 42
    The enrichment p-value should be 0.006445865568610187
    The deficiency p-value should be 0.9985899806396821

    If a = 21, b = 20, c = 34, d = 13
    The enrichment p-value should be 0.9883420938210076
    The deficiency p-value should be 0.0341612031176084
    """
    tables = [[[4, 6], [6, 4]],
              [[15, 8], [20, 42]],
              [[21, 20], [34, 13]]]
    accepted_p_values_list = [[0.9105522960012121, 0.3281408993483296],
                              [0.006445865568610187, 0.9985899806396821],
                              [0.9883420938210076, 0.0341612031176084]]
    failures = 0
    for table, accepted_p_values in zip(tables, accepted_p_values_list):
        deficiency_p_value = scipy.stats.fisher_exact(table, 'less')[1]
        enrichment_p_value = scipy.stats.fisher_exact(table, 'greater')[1]
        print()
        print(table[0])
        print(table[1])
        print("           Enrichment      |     Deficiency")
        print("Accepted", *accepted_p_values)
        print("Found   ", enrichment_p_value, deficiency_p_value)


def windowed_sequence_test():
    gene_names = ["ADAM2", "MCU", "RYR3", "TG", "KCP", "RTN4"]
    site_indexes = [10, 20, 30, 40, 50, 107]
    window_sizes = [3, 4, 5, 6, 5, 5]
    accepted_windowed_sequences = ["L L S G L G G",
                                   "R G G G G G G A G",
                                   "C I A T I H K E Q R K",
                                   "C E L Q R E T A F L K Q A",
                                   "L A G N S Q E Q W H P",
                                   "P E R Q P S W D P S P"]

    whole_sequences = sequence_search.get_sequences(gene_names)
    for name, index, size, whole_sequence, accepted_windowed_sequence in zip(
        gene_names, site_indexes, window_sizes, whole_sequences,
        accepted_windowed_sequences):

        print()
        print("Gene:    ", name)
        print("Accepted:", accepted_windowed_sequence)
        print("Found:   ", *sequence_search.windowed_sequence(whole_sequence, index, size))


def windowed_ranked_sequence_test():
    expected = ["EPSEVPTPKRP","LSLVAASPTLS","RRADNCSPVAE","PPYPQSRKLSY"]
    found = windowed_ranked_sequences.windowed_ranked_sequences(
        "test_data/ranked_sequences_test_data.txt", 5)
    for expected_sequence, found_sequence in zip(expected, found):
        print("Expected:", expected_sequence)
        print("Found:   ", found_sequence)
        print()


def amino_acid_substitution_test():
    sequence = "ARNDCEQGHILKMFPSTWYV"
    print("Sequence:             ", sequence)
    print("Expected substitution: AKQDCDQAHIIKIFPSSWFI")
    print("Found substitution:   ",
          "".join(processing.substitute_amino_acids(sequence)))


def letter_counts_test():
    columns = ["HHK", "DDY", "JYH", "SSD", "YTH"]
    expected_count_list = [["H", 2, "K", 1],
                           ["D", 2, "Y", 1],
                           ["J", 1, "Y", 1, "H", 1],
                           ["S", 2, "D", 1],
                           ["Y", 1, "T", 1, "H", 1]]
    expected_counts = [[0 for _ in range(ord("Z"))] for _ in range(5)]
    for i, row in enumerate(expected_count_list):
        for letter, count in zip(row[::2], row[1::2]):
            expected_counts[i][ord(letter)] = count
    found_counts = analysis.letter_counts(columns)
    print(found_counts)
    for i, row in enumerate(found_counts):
        for j, count in enumerate(row):
            if expected_counts[i][j] > 0 or count > 0:
                print("Expected:", chr(j), expected_counts[i][j])
                print("Found:", chr(j), count)
                print()
            if count != expected_counts[i][j]:
                print("TEST FAILED!")
                print()


def filtered_sequences_test():
    sequences = ["DHFJD", "DHFGT", "FGHTU", "KJFHY", "DHFHH", "DJFSP"]
    requirements = ((-2, [], ["D"]),
                    (1, ["J", "S"], []))
    expected_filtered_sequences = ["DHFGT", "DHFHH"]
    found_filtered_sequences = processing.filter_sequences(requirements,
                                                           sequences)
    print("Expected:", expected_filtered_sequences)
    print("Found:   ", found_filtered_sequences)


def most_significant_p_values_test():
    path = "test_data/sample-phosphoproteomic-data.txt"
    letter, window, index = "L", 4, 2
    
    names, sites, data_p_values = [], [], []
    with open(path) as f:
        for line in f.readlines()[1:]:
            row = line.split("\t")
            if row[3] == "P":
                names.append(row[0].split("-")[0])
                sites.append(int(row[2].split("|")[0][1:]))
                data_p_values.append(float(row[5]))
    
    sequences = windowed_ranked_sequences.windowed_ranked_sequences(
        names, sites, data_p_values, window
    )
    sequences = [sequence for sequence in sequences if len(sequence) >= 2*window]
    column = [sequence[index] for sequence in sequences]
    favorable = analysis.column_letter_counts(column)[ord(letter)]
    import time
    start = time.time()
    print(
        "Most significiant p-values:",
        p_values.most_significant_p_values(sequences, index, letter, favorable)
    )
    print(round(time.time() - start, 2))


def all_most_significant_p_values_test():
    path = "test_data/sample-phosphoproteomic-data.txt"
    window = 4; length = 2 * window
    offset = 2
    
    names, sites, data_p_values = [], [], []
    with open(path) as f:
        for line in f.readlines()[1:]:
            row = line.split("\t")
            if row[3] == "P":
                names.append(row[0].split("-")[0])
                sites.append(int(row[2].split("|")[0][1:]))
                data_p_values.append(float(row[5]))

    sequences = windowed_ranked_sequences.windowed_ranked_sequences(
        names, sites, data_p_values, window
    )
    sequences = [sequence for sequence in sequences if len(sequence) >= length]
    columns = [[sequence[i] for sequence in sequences] for i in range(length)]
    letter_counts = analysis.letter_counts(columns)
    result = p_values.all_most_significant_p_values(sequences, letter_counts)
    for i, row in enumerate(result): print(chr(i), row)
        
    


def main():
    #print("\n\nP-Value"); hypergeometric_test()
    #print("\n\nWindowed Sequence"); windowed_sequence_test()
    #print("\n\nWindowed Ranked Sequence\n"); windowed_ranked_sequence_test()
    #print("\nSubstitution test\n"); amino_acid_substitution_test()
    #print("\nLetter Counts Test\n"); letter_counts_test()
    #print("\nFiltered Sequences Test\n"); filtered_sequences_test()
    most_significant_p_values_test()
    #all_most_significant_p_values_test()
if __name__ == "__main__": main()
        
                           
