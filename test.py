import p_values


def hypergeometric_test():
    """
    Assume the letters a, b, c, and d represent the counts as displayed below.

      - +
    - a b
    + c d

    Then if a = 4, b = 6, c = 6, d = 4,
    The enrichment p-value should be 0.9105522960012121
    The deficiency p-value should be 0.3281408993483296

    If a = 20, b = 30, c = 25, d = 30,
    The enrichment p-value should be 0.7766662300662146
    The deficiency p-value should be 0.357198690677301

    If a = 15, b = 8, c = 20, d = 42
    The enrichment p-value should be 0.006445865568610187
    The deficiency p-value should be 0.9985899806396821

    If a = 21, b = 20, c = 34, d = 13
    The enrichment p-value should be 0.9883420938210076
    The deficiency p-value should be 0.0341612031176084
    """
    print("yp")
    tables = [[4, 6, 6, 4],
              [20, 30, 25, 20],
              [15, 8, 20, 42],
              [21, 20, 34, 13]]
    values = [[0.9105522960012121, 0.3281408993483296],
              [0.7766662300662146, 0.357198690677301],
              [0.006445865568610187, 0.9985899806396821],
              [0.9883420938210076, 0.0341612031176084]]
    for table, p_value in zip(tables, values):
        value = p_values.p_value(*table)
        print(value)
        if value != round(p_value[0], 5):
            print("Failed")
        else: print("Passed")
hypergeometric_test()
        
                           
