import p_values
import sequence_search
import processing
import analysis
import motif_search
import null_distribution

import scipy


def test_fisher_exact():
    """
    Assume the letters a, b, c, and d represent the counts 
    as displayed below.

      - +
    - a b
    + c d

    Then if a = 4, b = 6, c = 6, d = 4,
    the deficiency p-value should be 0.3281408993483299
    the enrichment p-value should be 0.9105522960012129

    If a = 15, b = 8, c = 20, d = 42
    the deficiency p-value should be 0.9985899806396819
    the enrichment p-value should be 0.006445865568610196

    If a = 21, b = 20, c = 34, d = 13
    the deficiency p-value should be 0.03416120311760853
    the enrichment p-value should be 0.9883420938210062
    """
    tables = [[[4, 6], [6, 4]],
              [[15, 8], [20, 42]],
              [[21, 20], [34, 13]]]
                                 
    accepted_p_values_list = [[0.3281408993483299, 0.9105522960012129],
                              [0.9985899806396819, 0.006445865568610196],
                              [0.03416120311760853, 0.9883420938210062]]
    for table, accepted_p_values in zip(tables, accepted_p_values_list):
        deficiency_p_value = scipy.stats.fisher_exact(table, 'less')[1]
        enrichment_p_value = scipy.stats.fisher_exact(table, 'greater')[1]
        assert (accepted_p_values[0] == deficiency_p_value
                and accepted_p_values[1] == enrichment_p_value)


def test_windowed_sequence():
    gene_names = ["ADAM2", "MCU", "RYR3", "TG", "KCP", "RTN4"]
    site_indexes = [10, 20, 30, 40, 50, 107]
    window_sizes = [3, 4, 5, 6, 5, 5]
    accepted_windowed_sequences = ["FLLSGLG",
                                   "SRGGGGGGA",
                                   "QCIATIHKEQR",
                                   "PCELQRETAFLKQ",
                                   "VLAGNSQEQWH",
                                   "APERQPSWDPS"]

    whole_sequences = sequence_search.get_sequences(gene_names)
    for name, index, size, whole_sequence, accepted_windowed_sequence in zip(
        gene_names, site_indexes, window_sizes, whole_sequences,
        accepted_windowed_sequences):
        assert (sequence_search.windowed_sequence(whole_sequence, index, size)
                == accepted_windowed_sequence)


def test_ranked_windowed_sequences():
    #return # this test needs rewriting
    data_path = "test_data/ranked_sequences_test_data.txt"
    expected = ["EPSEVPTPKRP","LSLVAASPTLS","RRADNCSPVAE","PPYPQSRKLSY"]
    window = 5

    names, indices, data_p_values = [], [], []
    with open(data_path) as f:
        for line in f.readlines()[1:]:
            row = line.split("\t")
            if row[3] == "P":
                names.append(row[0].split("-")[0])
                indices.append(int(row[2].split("|")[0][1:]))
                data_p_values.append(float(row[5]))
    found = sequence_search.ranked_windowed_sequences(
        names, indices, data_p_values, window
    )
    for expected_sequence, found_sequence in zip(expected, found):
        assert found_sequence == expected_sequence


def test_amino_acid_substitution():
    sequence = "ARNDCEQGHILKMFPSTWYV"
    expected_substitution = "AKQDCDQAHIIKIFPSSWFI"
    found_substitution = "".join(processing.substitute_amino_acids(sequence))
    assert expected_substitution == found_substitution


def test_letter_counts():
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
                assert expected_counts[i][j] == count


def test_filtered_sequences():
    sequences = ["DHFJD", "DHFGT", "FGHTU", "KJFHY", "DHFHH", "DJFSP"]
    requirements = ((-2, [], ["D"]),
                    (1, ["J", "S"], []))
    expected_filtered_sequences = ["DHFGT", "DHFHH"]
    found_filtered_sequences = processing.filter_sequences(requirements,
                                                           sequences)
    assert found_filtered_sequences == expected_filtered_sequences


def test_most_significant_p_values():
    return # this test needs rewriting
    # needs expected values
    data_path = "test_data/sample-phosphoproteomic-data.txt"
    letter, window, index = "A", 4, 2
    length = 2 * window
    step = 1024
    
    names, indices, data_p_values = [], [], []
    with open(data_path) as f:
        for line in f.readlines()[1:]:
            row = line.split("\t")
            if row[3] == "P":
                names.append(row[0].split("-")[0])
                indices.append(int(row[2].split("|")[0][1:]))
                data_p_values.append(float(row[5]))
                
    sequences = sequence_search.ranked_windowed_sequences(
        names, indices, data_p_values, window
    )
    sequences = [sequence for sequence in sequences if len(sequence) >= length]
    sequences = processing.substitute_amino_acids(sequences)
    column = [sequence[index] for sequence in sequences]
    favorable = analysis.column_letter_counts(column)[ord(letter)]
    print(
        "Most significant p-values:",
        p_values.most_significant_p_values(
            sequences, index, letter, favorable, step
        )
    )


def all_most_significant_p_values_test():
    # needs expected values
    path = "test_data/simulated-phosphoproteomic-data.txt"
    window = 2; length = 2 * window + 1
    step = 1024
    
    names, indices, data_p_values = [], [], []
    with open(path) as f:
        for line in f.readlines()[1:]:
            row = line.split("\t")
            if row[3] == "P":
                names.append(row[0].split("-")[0])
                indices.append(int(row[2].split("|")[0][1:]))
                data_p_values.append(float(row[5]))

    sequences = sequence_search.ranked_windowed_sequences(
        names, indices, data_p_values, window
    )
    sequences = [sequence for sequence in sequences if len(sequence) >= length]
    sequences = processing.substitute_amino_acids(sequences)
    columns = [[sequence[i] for sequence in sequences] for i in range(length)]
    letter_counts = analysis.letter_counts(columns)
    results = p_values.all_most_significant_p_values(
        sequences, letter_counts, step
    )
    
    for c, column in enumerate(results):
        print(c - window if c < window else "+" + str(c - window + 1))
        for p, pair in enumerate(column):
            if pair != (0,0): 
                print(chr(p), "|", f'{pair[0]:.6f}', "|", f'{pair[1]:.6f}')
        print()
    
    
def null_distribution_test():
    # needs expected values
    path = "test_data/simulated-phosphoproteomic-data.txt"
    window = 1; length = 2 * window + 1
    step = 1024
    repetitions = 10
    seed = 0

    if repetitions > 10: return
    
    names, indices, = [], []
    with open(path) as f:
        for line in f.readlines()[1:]:
            row = line.split("\t")
            if row[3] == "P":
                names.append(row[0].split("-")[0])
                indices.append(int(row[2].split("|")[0][1:]))
                
    sequences = sequence_search.get_sequences(names)
    sequences = [sequence_search.windowed_sequence(sequence, index, window)
                for (sequence, index) in zip(sequences, indices)]
    sequences = [sequence for sequence in sequences if len(sequence) >= length]
    sequences = processing.substitute_amino_acids(sequences)
    results = null_distribution.null_distribution(
        sequences, window, repetitions, step, seed
    )
        

def iterative_motif_search_test():
    # needs expected values
    path = "test_data/simulated-phosphoproteomic-data.txt"
    window = 5; length = 2 * window + 1
    step = 1024
    threshold = 0.0005
    
    names, indices, = [], []
    with open(path) as f:
        for line in f.readlines()[1:]:
            row = line.split("\t")
            if row[3] == "P":
                names.append(row[0].split("-")[0])
                indices.append(int(row[2].split("|")[0][1:]))
                
    sequences = sequence_search.get_sequences(names)
    sequences = [sequence_search.windowed_sequence(sequence, index, window)
                for (sequence, index) in zip(sequences, indices)]
    sequences = [sequence for sequence in sequences if len(sequence) >= length]
    sequences = processing.substitute_amino_acids(sequences)
    
    print("Key:")
    print("motif => [newfound_motifs]")
    print("(index, letter, presence)")
    print("Index is relative to 0 at the left, letter is the amino acid,")
    print("and presence is whether the acid must appear (True) or absent (False)")
    print("Final Graph", motif_search.motif_search(sequences, step, threshold))
    
                        


def main():
    #print("\n\nWindowed Sequence"); windowed_sequence_test()
    #print("\n\Ranked Windowed Sequences Test\n"); ranked_windowed_sequences_test()
    #print("\nSubstitution Test\n"); amino_acid_substitution_test()
    #print("\nLetter Counts Test\n"); letter_counts_test()
    #print("\nFiltered Sequences Test\n"); filtered_sequences_test()
    #print("\nMost Significant P Values Test\n"); most_significant_p_values_test()
    #print("\nAll Most Significant P Values Test\n"); all_most_significant_p_values_test()
    #print("\nNull Distribution Test\n"); null_distribution_test()
    print("\nIterative Motif Search:"); iterative_motif_search_test()
if __name__ == "__main__": main()
        
                           
