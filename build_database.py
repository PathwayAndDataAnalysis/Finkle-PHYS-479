def read_gene_list(filename = "data/uniprot_sprot.fasta"):
    """
    Return a list of the names and sequences of the genes in a .fasta
    file.

    Defaults to the uniprot_sprot database.

    :param filename: The name of the .fasta file containing the amino
                     acid sequences of the genes to be studied
    :type filename: str
    :returns A list of the names and sequences in the .fasta file.
    :rtype list[str, str]:
                     
    """
    with open(filename) as f:
        return f.read().split(">")

    
def extract_names_and_sequences(gene_list, organism = "Homo sapiens"):
    """
    Returns one list of the names and another of the sequences of a
    unified list of gene names and sequences for an organism.

    :param gene_list: a list of genes and their names
    :type gene_list: list[list[str], list[str]]
    :param organism: organism of which genes are to be selected,
                     defaulting to human
    :type organism: str
    :returns: a list of gene names and associated sequences
    :rtype: list[list[str], list[str]]
    """
    names, sequences = [], []
    unique_names = set()
    for entry in gene_list[1:]:
        lines = entry.split("\n")
        annotation = lines[0]
        if not ("OS=" + organism in annotation and "GN=" in annotation): continue
        name = annotation.split("GN=")[1].split(" ")[0]
        if name in unique_names: continue
        unique_names.add(name)
        fasta_sequence = "".join(lines[1:])
        names.append(name)
        sequences.append(fasta_sequence)
    return names, sequences
        
    
def save_names_and_sequences(names, sequences):
    """
    Save the names and sequences of a list of genes to a text file.

    :param names: The names of the genes
    :type names: list[str]
    :param sequences: The amino acid sequences of the genes
    :type sequences: list[str]
    """
    names_and_indices = [[name, n] for n, name in enumerate(names)]
    with (open('data/raw_sequences.txt', 'w') as sequences_file,
          open('data/gene_names.txt', 'w') as names_file):
        for name, index in names_and_indices:
            names_file.write(name + "\n")
            sequences_file.write(sequences[index] + "\n")


def main():
    """
    Build a concise text file database of the names and sequences of
    the genes of an organism from a comprehensive .fasta file
    """
    gene_list = read_gene_list()
    names, sequences = extract_names_and_sequences(gene_list)
    save_names_and_sequences(names, sequences)
if __name__ == "__main__": main()
