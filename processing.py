def substitute_amino_acids(sequences):
    """
    Return the sequences with their amino acids replaced as specified 
    in amino-acid-mapping.tsv
    """
    with open("amino-acid-mapping.tsv") as f:
        substitutions = {line.split("\t")[0] : line.split("\t")[1][:-1]
                         for line in f.readlines()[1:]}

    new_sequences = []
    for sequence in sequences:
        try:
            new_sequences.append("".join(
                substitutions[amino_acid] for amino_acid in sequence
            ))
        except:
            pass
            
    return new_sequences


def filter_sequences(requirements, sequences):
    """
    Return the sequences that pass a filter of requirements.

    Requirements must be an list or tuple of the form:
    (
        (index, letter, presence_or_absence, letter, presence_or_absence ...),
        (index, letter, presence_or_absence, letter, presence_or_absence ...),
        ...
    )
    
    :param requirements: letters required to be present or absent at 
                         each index of a sequence
    :type requirements: tuple(tuple(int, int, boolean, int, boolean, ...))
    :param sequences: amino acid sequences, the letters of which are to be
                      filtered
    :type sequences: list[str]
    :return: the sequences that pass the filters
    :rtype: list[str]
    """
    filtered_sequences = []
    middle = len(sequences[0]) // 2
    for sequence in sequences:        
        passes_filter = True
        for row in requirements:
            letter = sequence[middle + row[0]]            
            for prohibited_character in row[1]:
                if letter == prohibited_character:
                    passes_filter = False; break
            if not passes_filter: break
            for required_character in row[2]:               
                if not letter == required_character:
                    passes_filter = False; break            
        if passes_filter: filtered_sequences.append(sequence)
    return filtered_sequences
