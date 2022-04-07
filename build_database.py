def read_gene_list(filename = "data/uniprot_sprot.fasta"):
    """
    Return a list of the names and sequences of the genes in a fasta file.

    Defaults to the uniprot_sprot database.
    """
    with open(filename) as f:
        return f.read().split(">")

    
def extract_names_and_sequences(gene_list):
    """
    Returns one list of the names and another of the sequences of a unified list
    of gene names and sequences.
    """
    names, sequences = [], []
    unique_names = set()
    for entry in gene_list[1:]:
        lines = entry.split("\n")
        annotation = lines[0]
        if not ("OS=Homo sapiens" in annotation and "GN=" in annotation): continue
        print(annotation)
        name = annotation.split("GN=")[1].split(" ")[0]
        if name in unique_names: continue
        unique_names.add(name)
        fasta_sequence = "".join(lines[1:])
        names.append(name)
        sequences.append(fasta_sequence)
    return names, sequences
        
    
def save_names_and_sequences(names, sequences):
    """
    
    """
    names_and_indices = [[name, n] for n, name in enumerate(names)]
    names_and_indices.sort()
    with open('data/raw_sequences.txt', 'w') as sequences_file, open('data/gene_names.txt', 'w') as names_file:
        for name, index in names_and_indices:
            names_file.write(name + "\n")
            sequences_file.write(sequences[index] + "\n")


def main():
    gene_list = read_gene_list()
    names, sequences = extract_names_and_sequences(gene_list)
    save_names_and_sequences(names, sequences)
if __name__ == "__main__": main()
