import scipy.stats
import numpy as np
import pandas as pd
import itertools
import math
import graph, comb


def main (membership_child, membership_parent, mixture_total, **kwargs):
    subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(list( membership_child ), mixture_total)   # Find all combinatory sums

    output_file = open(kwargs["PHYLOGENY_DIR"], "w")

    g = graph.UnidirectedGraph()
    completed_parent = []

    print ("\n")

    while len(completed_parent) < len( membership_parent ):  
        if len(subset_list_acc) == 0:  
            break

        
        p_maxmax = float("-inf"); subset_list_maxmax = []; subset_mixture_maxmax = []; sum_mixture_maxmax = []
        for j1 in sorted( list( membership_parent ) ) :           # Iterating parent clusters
            if j1 in completed_parent:
                continue

            parent_element_mixture = mixture_total[:,j1]
            p_max = float("-inf"); subset_list_max = []; subset_mixture_max = []; sum_mixture_max = []

            for j2 in range(len(subset_mixture_acc)):
                subset_list = subset_list_acc[j2]
                subset_mixture = subset_mixture_acc[j2]
                sum_mixture = sum_mixture_acc[j2]

                p = 0
                for i in range (kwargs["NUM_BLOCK"]):
                    depth = 100
                    a = int(sum_mixture[i] * 100 / 2) 
                    b = depth - a
                    target_a = int (parent_element_mixture[i] * 100/ 2)
                    try:
                        p = p + math.log10(scipy.stats.betabinom.pmf(target_a, depth, a + 1, b+1))
                    except:
                        p = p - 400
                        
                if p > p_max:
                    p_max = p
                    subset_list_max = subset_list
                    subset_mixture_max = subset_mixture
                    sum_mixture_max = sum_mixture

            if p_max > p_maxmax:              
                p_maxmax = p_max
                j_maxmax = j1 
                subset_list_maxmax = subset_list_max
                subset_mixture_maxmax = subset_mixture_max
                sum_mixture_maxmax = sum_mixture_max

        # Now, exclude from the list
        completed_parent.append (j_maxmax)
        subset_list_acc.remove(subset_list_maxmax)
        subset_mixture_acc.remove(subset_mixture_maxmax)
        sum_mixture_acc.remove(sum_mixture_maxmax)

        print ("\t\tparent No = {0}, parent_mixture = {1}, sum_mixture = {2}, subset_list = {3},  p = {4}".format(j_maxmax,  mixture_total[:,j_maxmax], np.round(sum_mixture_maxmax, 2), subset_list_maxmax, round (p_maxmax, 2) ), file = output_file )
        print ("\t\tparent No = {0}, parent_mixture = {1}, sum_mixture = {2}, subset_list = {3},  p = {4}".format(j_maxmax,  mixture_total[:,j_maxmax], np.round(sum_mixture_maxmax, 2), subset_list_maxmax, round (p_maxmax, 2) ) )
        g.intervene (j_maxmax, subset_list_maxmax)           
        for proband_clone_index in subset_list_maxmax:      # 4 -6,  4- 7
            g.addEdge(j_maxmax, proband_clone_index)

    print ("\t\t", file = output_file)
    print ("\t\t")
    for root in completed_parent:
        if g.findparent(root) == "None":   # Run only if root node
            g.dfs(root, kwargs["PHYLOGENY_DIR"])
        print ("\t\t", file = output_file)
        print ("\t\t")


    output_file.close()
    return g