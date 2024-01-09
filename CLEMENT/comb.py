import numpy as np
import pandas as pd
import itertools

def comb2 (t):     # Input : 2D ndarray (mixture) -> Combinatory sum of two columns
    comb = []
    combi = list(itertools.combinations(range(t.shape[1]), 2))
    for a, b in combi:
        comb.append(list(t[:,a] + t[:,b]))
    return np.array(pd.DataFrame(comb).sort_values(0))

def comb3 (t):
    comb = []
    combi = list(itertools.combinations(range(t.shape[1]), 3))
    for a, b, c in combi:
        comb.append(list(t[:,a] + t[:,b] + t[:,c]))
    return np.array(pd.DataFrame(comb).sort_values(0))

def comball (ss, t):
    combi =  itertools.chain(*map(lambda x: itertools.combinations(ss, x), range(2, len(ss)+1)))
    subset_list_acc = []
    subset_ndarray_acc = []
    sum_acc = []
    for subset in combi:
        subset_list = list (subset)

        subset_ndarray = []
        for i in subset_list:
            subset_ndarray.append (list(t[:,i]))

        subset_list_acc.append(subset_list)
        subset_ndarray_acc.append (subset_ndarray)
        sum_acc.append (list(np.sum(subset_ndarray, axis = 0)))

    return (  subset_list_acc, subset_ndarray_acc , sum_acc )


def main(t, n):
    if n == 2:
        return comb2(t)
    elif n == 3:
        return comb3(t)
