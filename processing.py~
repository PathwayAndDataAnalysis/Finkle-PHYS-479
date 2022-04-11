def substitute_amino_acids(sequences):
    with open("amino-acid-mapping.tsv") as f:
        substitutions = {line.split("\t")[0] : line.split("\t")[1] for line in f.readlines()}
    return [[substitutions[character] for character in lines] for sequence in sequences]
