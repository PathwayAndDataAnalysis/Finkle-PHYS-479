from scipy.stats import hypergeom


def column_letter_counts(column):
    """Return an array of counts of letters in an amino acid column."""
    counts = [0 for _ in range(ord("Z"))]
    for letter in column:
        counts[ord(letter)] += 1
    return counts


def letter_counts(columns):
    """Return a 2d array of counts of letters in columns of aligned proteins."""
    width = len(columns)
    middle = width // 2
    relevant = [i for i in range(0, middle)] + [i for i in range(middle, width)]
    return [column_letter_counts(columns[column]) for column in relevant]
