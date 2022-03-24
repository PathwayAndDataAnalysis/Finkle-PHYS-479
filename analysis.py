from scipy.stats import hypergeom


def p_value(selected_column,
            remaining_column,
            selected_column_count,
            remaining_column_count):
    """
    Return the sum of the hypergeometric probabilities of more 
    favorable results in the same number of trials.
    """
    selected, remaining = len(selected_column), len(remaining_column)
    favorable = selected_column_count + remaining_column_count
    remaining_unfavorable = remaining - remaining_column_count
    trials = selected + remaining
    limit = favorable + max(remaining_unfavorable - favorable, 0)
    possibilities = [i for i in range(favorable, limit)]
    return sum(hypergeom.pmf(possibilities, trials, favorable, selected))


def column_letter_counts(column):
    """Return an array of counts of letters in an amino acid column."""
    letter_counts = [0 for _ in range(ord("Y"))]
    for letter in column: letter_counts[ord(letter)] += 1
    return letter_counts


def letter_counts(sequences):
    """Return a 2d array of counts of letters in aligned proteins."""
    width = len(sequences)
    middle = length // 2
    columns = [i for i in range(0, middle)] + [i for i in range(middle, width)]
    letter_counts = [column_letter_counts(sequences[column] for column in columns]
    return letter_counts
