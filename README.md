# Finkle-PHYS-479

## Work in progress: not feature complete or integrated

## Setup:
1 Download Python
2 Open your command line and keep it open
3 Change directories to your Desktop
4 Type: python3 -m venv your_project_name
5 Type: source bin/activate
6 Type: pip3 install scipy
7 Open to the folder you have created
8 Create a folder called Data
9 Download the UniProt gene sequence database as a .fasta file
10 Move the .fasta file to the Data folder you just created
11 Return to your command line
12 Type: python3 -m build_database
13 Type: python3 -m test
14 If all the expected and found results align, you are ready to use the module in your own Python code!


This project has two features:

- Ordering a list of amino sequences by p-value, from smallest positive to largest positive followed by negatives likewise
- Calculating, for one selected and another remaining set of windowed amino acid sequences centered and aligned on phosphorylation site, the p-value of the count of each amino acid at each position left or right of the site in the selected set.

The features are useable in a python script.

To use the first feature import windowed_ranked_sequences into your python file and call ranked_sequences.ranked_sequences(your/sequence/file/path) where your file is a tab-separated (.TSV) file

To use the second feature, import p_values into your python file and call p_values.sequence_p_values(selected_sequences, remaining_sequences) where the sequences are lists or tuples


## Test Examples

To run the test, change directories to the file containing the code and type python3 -m test

p_values.hypergeometric_test(selected_sequences, remaining_sequences)

Assume the letters a, b, c, and d represent the counts as displayed below.

  - +
- a b
+ c d

Then if a = 4, b = 6, c = 6, d = 4,
The enrichment p-value should be 0.9105522960012121
The deficiency p-value should be 0.3281408993483296

If a = 15, b = 8, c = 20, d = 42
The enrichment p-value should be 0.006445865568610187
The deficiency p-value should be 0.9985899806396821

If a = 21, b = 20, c = 34, d = 13
The enrichment p-value should be 0.9883420938210076
The deficiency p-value should be 0.0341612031176084


sequence_search.windowed_sequence(sequence, site_number, window_width)

Gene:     ADAM2
Accepted: F L L S G L G
Found:    F L L S G L G

Gene:     MCU
Accepted: L L L L S S R G G
Found:    L L L L S S R G G

Gene:     RYR3
Accepted: L E Q S L S V R A L Q
Found:    L E Q S L S V R A L Q

Gene:     TG
Accepted: C Q N D G R S C W C V G A
Found:    C Q N D G R S C W C V G A

Gene:     KCP
Accepted: V R Q L E S C E C H P
Found:    V R Q L E S C E C H P


windowed_ranked_sequences.windowed_ranked_sequences(sequences, site_number, window_width)

Expected: SQKEPSEVPTP
Found:    SQKEPSEVPTP

Expected: NDPRCSTSNNR
Found:    NDPRCSTSNNR

Expected: KGVSMSLPSSP
Found:    KGVSMSLPSSP

Expected: LGSTKSLNHSK
Found:    LGSTKSLNHSK