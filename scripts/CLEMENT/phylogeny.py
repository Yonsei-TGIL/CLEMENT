import scipy.stats
import numpy as np
import pandas as pd
import itertools
import math
import graph, comb


def main (membership_child, membership_outside, mixture_total, **kwargs):
    subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(list( membership_child ), mixture_total)   # Listing all combination

    output_file = open(kwargs["PHYLOGENY_DIR"], "w")

    g = graph.UnidirectedGraph()
    completed_outside = []

    print (subset_list_acc)

    while len(completed_outside) < len( membership_outside ):          # Iterate for Outside cluster
        if len(subset_list_acc) == 0:  # If child clone is less than the outside cluster ( = putative parent clone)
            break

        
        p_maxmax = float("-inf"); subset_list_maxmax = []; subset_mixture_maxmax = []; sum_mixture_maxmax = []
        for j1 in sorted( list( membership_outside ) ) :           # Iterating outside cluster
            if j1 in completed_outside:
                continue

            outside_element_mixture = mixture_total[:,j1]
            p_max = float("-inf"); subset_list_max = []; subset_mixture_max = []; sum_mixture_max = []

            for j2 in range(len(subset_mixture_acc)):      # Iterating other cluster
                subset_list = subset_list_acc[j2]
                subset_mixture = subset_mixture_acc[j2]
                sum_mixture = sum_mixture_acc[j2]

                p = 0
                for i in range (kwargs["NUM_BLOCK"]):
                    depth = 100
                    a = int(sum_mixture[i] * 100 / 2) 
                    b = depth - a
                    target_a = int (outside_element_mixture[i] * 100/ 2)
                    try:
                        p = p + math.log10(scipy.stats.betabinom.pmf(target_a, depth, a + 1, b+1))
                    except:
                        p = p - 400
                        
                if p > p_max:
                    p_max = p
                    subset_list_max = subset_list
                    subset_mixture_max = subset_mixture
                    sum_mixture_max = sum_mixture

            if p_max > p_maxmax:              # Select the outside cluster with highest probability
                p_maxmax = p_max
                j_maxmax = j1 
                subset_list_maxmax = subset_list_max
                subset_mixture_maxmax = subset_mixture_max
                sum_mixture_maxmax = sum_mixture_max

        # Now, exclude the combination from the list
        completed_outside.append (j_maxmax)
        subset_list_acc.remove(subset_list_maxmax)
        subset_mixture_acc.remove(subset_mixture_maxmax)
        sum_mixture_acc.remove(sum_mixture_maxmax)

        print ("outside No = {0}, outside_mixture = {1}, sum_mixture = {2}, subset_list = {3},  p = {4}".format(j_maxmax,  mixture_total[:,j_maxmax], sum_mixture_maxmax, subset_list_maxmax, p_maxmax ), file = output_file )
        print ("outside No = {0}, outside_mixture = {1}, sum_mixture = {2}, subset_list = {3},  p = {4}".format(j_maxmax,  mixture_total[:,j_maxmax], sum_mixture_maxmax, subset_list_maxmax, p_maxmax ) )
        g.intervene (j_maxmax, subset_list_maxmax)            # If father - grandfather?
        for proband_clone_index in subset_list_maxmax:      # 4 -6,  4- 7
            g.addEdge(j_maxmax, proband_clone_index)

    print ("", file = output_file)
    print ("")
    for root in completed_outside:
        if g.findparent(root) == "None":   # Run if only  root node.  If not, pass
            g.dfs(root, kwargs["PHYLOGENY_DIR"])
        print ("", file = output_file)
        print ("")


    output_file.close()
    return g