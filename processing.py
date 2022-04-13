def substitute_amino_acids(sequences):
    """
    Return the sequences with their amino acids replaced as specified in
    amino-acid-mapping.tsv
    """
    with open("amino-acid-mapping.tsv") as f:
        substitutions = {line.split("\t")[0] : line.split("\t")[1][:-1]
                         for line in f.readlines()}
    return ["".join([substitutions[amino_acid] for amino_acid in sequence])
             for sequence in sequences]


def filter_sequences(requirements, sequences):
    """
    Return the sequences that pass a filter of requirements.

    Requirements must be an list or tuple of the form:
    (
        (column, letter, presence_or_absence, letter, presence_or_absence ...),
        (column, letter, presence_or_absence, letter, presence_or_absence ...),
        ...
    )
    """
    filtered_sequences = []
    for sequence in sequences:
        for row in requirements:
            passes_filter = True
            for i in range(1, len(row) - 1, 2):
                if sequence[row[0]] == row[i] != row[i+1]:
                    passes_filter = False; break
            if passes_filter: filtered_sequences.append(sequence)
