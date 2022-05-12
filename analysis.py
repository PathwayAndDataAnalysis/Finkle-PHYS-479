from scipy.stats import hypergeom


def column_letter_counts(column):
    """
    Return a list of counts of letters in an amino acid column.
    
    :param column: a list of amino acids to be counted
    :type column: list[str]
    :return: list of counts of letters in an amino acid column
    :rtype: tuple[int]
    """
    counts = [0 for _ in range(ord("Z"))]
    for letter in column: counts[ord(letter)] += 1
    return tuple(counts)


def letter_counts(columns):
    """
    Return a 2d array of counts of letters in columns of aligned proteins.
    
    :param columns: a list of lists of amino acids to be counted
    :type columns: list[list[str]]
    :return: list of counts of letters in each column
    :rtype: tuple(tuple(int))
    """
    return tuple(column_letter_counts(column) for column in columns)
    
