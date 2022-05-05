import analysis
import p_values
import processing

import copy


def match(sequence, index, letter, presence):
    return (sequence[index] == letter) == presence


def motif_search(passed_sequences, step, threshold, motif = None):
    """
    Return a directed acyclic graph of any motifs in the sequences.
    
    :param passed_sequences: the sequences of which motifs are to be
                             graphed
    :type sequences: list[str]
    :param step: the step size of letters by which to increment a
                 p-value sweep of a column
    :type step: int
    :param threshold: the maximum p-value of the presence or absence of
                      a letter to be considered a motif
    :type threshold: float
    :param motif: an index, letter, and the presence or absence of the 
                  letter at the column
    :return: a directed acyclic graph of any motifs in the sequences
    """

    # create a deep copy of sequences
    sequences = copy.deepcopy(passed_sequences)

     
    # If a motif is given, exclude any sequence without it.
    if motif:
        sequences = [sequence for sequence in sequences
                     if match(sequence, *motif)]
    
    # Calculate the p values of the letter counts of the remaining sequences
    columns = tuple(tuple(sequence[i] for sequence in sequences) 
                    for i in range(len(sequences[0])))
    letter_counts = analysis.letter_counts(columns)
    results = p_values.all_most_significant_p_values(
        sequences, letter_counts, step
    )
    
    # Find any new motifs
    newfound_motifs = []
    for index, column in enumerate(results):
        for letter, p_value_pair in enumerate(column):
            if p_value_pair[0] > threshold or p_value_pair[1] > threshold:
                if p_value_pair[0] < threshold:
                    newfound_motifs.append(tuple((index, chr(letter), False)))
                if p_value_pair[1] < threshold:
                    newfound_motifs.append(tuple((index, chr(letter), True)))
    print()
    print(len(sequences))
    print(motif, "=>", newfound_motifs)
    # If no new motif has been found, return the given motif,
    # Else, return the next branch in the graph for each new motif
    return {motif : motif_search(sequences, step, threshold, newfound_motif)
            for newfound_motif in newfound_motifs} if newfound_motifs else motif
        
