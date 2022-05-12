import analysis
import p_values
import processing

import copy


def matches(sequence, index, letter, presence):
    """
    Return whether the character of the sequence at the index satisfies
    the requirement of the presence or absence of a letter.

    :param sequence: an amino acid sequence
    :type sequence: str
    :param index: a position in the sequence, where 0 is left-most
    :type index: int
    :param letter: the amino acid character to be matched
    :type letter: str
    :param presence: whether the character must be present (True) or
                     must be absent (False)
    """
    return (sequence[index] == letter) == presence


def filtered_sequences(sequences, motif):
    """
    Return the sequences that have a motif.

    :param sequences: amino acid sequences of equal length and aligned
                      on a site
    :type sequences: list[str]
    :param motif: an index, letter, and the presence or absence of the 
                  letter at the index in each sequence
    :type motif: tuple(int, str, bool)
    :return: the sequences that have the motif
    :rtype: list[str]
    
    """
    return [sequence for sequence in sequences if matches(sequence, *motif)]


def enrichment_and_deficiency_p_values(sequences, step):
    """
    Return the enrichment and deficiency p-values of the letter counts
    of the sequences.

    Reorganize the sequences from rows to columns, count how many
    times each letter appears in each column, and return the
    enrichment and deficiency p-value of each count.

    :param sequences: amino acid sequences of equal length and aligned
                      on a site
    :type sequences: list[str]
    
    """
    columns = tuple(tuple(sequence[i] for sequence in sequences) 
                    for i in range(len(sequences[0])))
    letter_counts = analysis.letter_counts(columns)
    return p_values.all_most_significant_p_values(
        sequences, letter_counts, step
    )


def newfound_motifs(sequences, enrichments_and_deficiencies, threshold):
    """
    Return any new motifs in the sequences.

    :param sequences: windowed amino acid sequences of equal length and
                      aligned on a site
    :type sequences: list[str]
    :param enrichments_and_deficiencies: the enrichment and deficiency
                                         p-value of every amino acid
                                         at every index of the sequences
    :type enrichments_and_deficiencies: tuple(tuple(float, float))
    """
    newfound_motifs = []
    for index, column in enumerate(enrichments_and_deficiencies):
        for letter, p_value_pair in enumerate(column):
            if p_value_pair[0] > threshold or p_value_pair[1] > threshold:
                if p_value_pair[0] < threshold:
                    newfound_motifs.append(tuple((index, chr(letter), False)))
                if p_value_pair[1] < threshold:
                    newfound_motifs.append(tuple((index, chr(letter), True)))
    return newfound_motifs


def motif_search(sequences, step, threshold, motif = None):
    """
    Return a directed acyclic graph of any motifs in the sequences.

    Find any motifs in the sequences and then, recursively, find any
    motifs in the sequences that match that motif, tracking each 'path'
    of this search on a directed acyclic graph, wherein each node is
    a filtered list of sequences and each edge is a motif found therein.
    Return this graph in a convenient data structure.

    A motif is a letter, the presence or absence of which at an index
    in these aligned sequences has been found to be a pattern, and the
    sequences are to be filtered down to the ones exhibiting it.
    
    :param sequences: windowed amino acid sequences of equal length and
                      aligned on a site
    :type sequences: list[str]
    :param step: the step size of letters by which to increment a
                 p-value sweep of a column
    :type step: int
    :param threshold: the maximum p-value of the presence or absence of
                      a letter to be considered a motif
    :type threshold: float
    :param motif: an index, letter, and the presence or absence of the 
                  letter at the index in each sequence
    :type motif: tuple(int, string, bool)
    :return: a directed acyclic graph of any motifs in the sequences
    """

    if motif: sequences = filtered_sequences(sequences, motif)

    enrichments_and_deficiencies = enrichment_and_deficiency_p_values(
        sequences, step
    )
    
    motifs = newfound_motifs(sequences, enrichments_and_deficiencies, threshold)    
    print()
    print(len(sequences))
    print(motif, "=>", motifs)
    # If no new motif has been found, return the given motif,
    # Else, return the next branch in the graph for each new motif
    return {motif : motif_search(sequences, step, threshold, newfound_motif)
            for newfound_motif in motifs} if motifs else motif
        
