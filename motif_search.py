import analysis
import p_values
import processing

def motif_search(sequences, step, threshold, motif = None):
    """
    Return a directed acyclic graph of any motifs in the sequences.
    
    :param sequences: the sequences of which motifs are to be graphed
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
    
    # If a motif is given, exclude any sequence without it.
    if motif:
        index, letter, requirement = motif
        old_length = len(sequences)
        for s, sequence in enumerate(sequences):
            if sequence[index] == letter and requirement == False: 
                del sequences[s]
            elif sequence[index] != letter and requirement == True: 
                del sequences[s]
        # If no sequence has been excluded, instead return the motif
        if len(sequences) == old_length: return motif
    
    # Calculate the p values of the letter counts of the remaining sequences
    columns = tuple(tuple(sequence[i] for sequence in sequences) 
                    for i in range(len(sequences[0])))
    letter_counts = analysis.letter_counts(columns)
    results = p_values.all_most_significant_p_values(
        sequences, letter_counts, step
    )
    
    # Find any new motifs
    motifs = []
    for c, column in enumerate(results):
        for p, pair in enumerate(column):
            if pair[0] < threshold and pair[1] < threshold: continue
            if pair[0] < threshold: motifs.append(tuple((c, p, False)))
            elif pair[1] < threshold: motifs.append(tuple((c, p, True)))
    
    # If no new motif has been found, return the given motif
    if motifs == []: return motif
    # Return the next branch in the graph for each new motif
    return {motif : motif_search(sequences, step, threshold, motif) 
                 for motif in motifs}
        