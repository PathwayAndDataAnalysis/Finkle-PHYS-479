def substitute_amino_acids(sequences):
    with open("amino-acid-mapping.tsv") as f:
        substitutions = {line.split("\t")[0] : line.split("\t")[1][:-1]
                         for line in f.readlines()}
    return ["".join([substitutions[amino_acid] for amino_acid in sequence])
             for sequence in sequences]
