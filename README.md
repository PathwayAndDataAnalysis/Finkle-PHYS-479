# Finkle-PHYS-479

## Work in progress: not feature complete or integrated

## Required Software: Command Line, Python, SciPy

This project has two features:

- Ordering a list of amino sequences by p-value, from smallest positive to largest positive followed by negatives likewise
- Calculating, for one selected and another remaining set of windowed amino acid sequences centered and aligned on phosphorylation site, the p-value of the count of each amino acid at each position left or right of the site in the selected set.

The features are useable in a python script.

To use the first feature import ranked_sequences into your python file and call ranked_sequences.ranked_sequences(your/sequence/file/path) where your file is a tab-separated (.TSV) file

To use the second feature, import ranked_sequences into your python file and call p_values.sequence_p_values(selected_sequences, remaining_sequences) where the sequences are lists or tuples