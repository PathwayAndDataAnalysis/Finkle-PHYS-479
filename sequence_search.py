import mmap


def read_gene_list(filename = "uniprot_sprot.fasta"):
    with open(filename) as f:
        return f.read().split(">")

    
def extract_names_and_sequences(gene_list):
    names, sequences = [], []
    greatest_name_length = 0
    greatest_sequence_length = 0
    for gene_entry in gene_list[1:]:
        lines = gene_entry.split("\n")
        gene_annotation = lines[0]
        if "OS=Human" not in gene_annotation: continue
        if "GN=" not in gene_annotation: continue
        else: gene_name = gene_annotation.split("GN=")[1].split(" ")[0]
        fasta_sequence = "".join(lines[1:])
        names.append(gene_name)
        sequences.append(fasta_sequence)
    return names, sequences
        
    
def save_names_and_sequences(names, sequences):
    names_and_indices = [[name, n] for n, name in enumerate(names)]
    names_and_indices.sort()
    with open('raw_sequences.txt', 'w') as sequences_file, open('gene_names.txt', 'w') as names_file:
        for name, index in names_and_indices:
            names_file.write(name + "\n")
            sequences_file.write(sequences[index] + "\n")
        
            
def get_sequence(gene_name):
    with open('gene_names.txt') as gene_names:
        for i, name in gene_names:
            if name == gene_name:
                sequence_number = i
                break
    with open ('raw_sequences.txt') as sequences:
        for i, sequence in sequences:
            if i == sequence_number:
                return sequence
        
