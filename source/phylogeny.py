import scipy.stats
import numpy as np
import pandas as pd
import itertools
import math
import graph, comb


def main (membership_child, membership_parent, mixture_total, **kwargs):
    subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(list( membership_child ), mixture_total)   # 모든 덧셈 조합을 구하기

    output_file = open(kwargs["PHYLOGENY_DIR"], "w")

    g = graph.UnidirectedGraph()
    completed_parent = []

    print ("\n")

    while len(completed_parent) < len( membership_parent ):          # parent cluster의 개수만큼 돈다 (OUtlier가 있는 데이터일 경우 당연히 맨 마지막에 있는 outlier는 안 돈다)
        if len(subset_list_acc) == 0:  # parent cluster는 넘쳐나도 proband clone이 적어서 조합에 한계가 있으면
            break

        
        p_maxmax = float("-inf"); subset_list_maxmax = []; subset_mixture_maxmax = []; sum_mixture_maxmax = []
        for j1 in sorted( list( membership_parent ) ) :           # parent cluster를 돈다
            if j1 in completed_parent:
                continue

            parent_element_mixture = mixture_total[:,j1]
            p_max = float("-inf"); subset_list_max = []; subset_mixture_max = []; sum_mixture_max = []

            for j2 in range(len(subset_mixture_acc)):      #여러 조합을 돈다
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

            if p_max > p_maxmax:              # parent cluster중에 가장 신뢰도 있는 값을 보인 애를 선택
                p_maxmax = p_max
                j_maxmax = j1 
                subset_list_maxmax = subset_list_max
                subset_mixture_maxmax = subset_mixture_max
                sum_mixture_maxmax = sum_mixture_max

        # 제일 높은 애는 이제 list에서 제외
        completed_parent.append (j_maxmax)
        subset_list_acc.remove(subset_list_maxmax)
        subset_mixture_acc.remove(subset_mixture_maxmax)
        sum_mixture_acc.remove(sum_mixture_maxmax)

        print ("parent No = {0}, parent_mixture = {1}, sum_mixture = {2}, subset_list = {3},  p = {4}".format(j_maxmax,  mixture_total[:,j_maxmax], sum_mixture_maxmax, subset_list_maxmax, p_maxmax ), file = output_file )
        print ("parent No = {0}, parent_mixture = {1}, sum_mixture = {2}, subset_list = {3},  p = {4}".format(j_maxmax,  mixture_total[:,j_maxmax], sum_mixture_maxmax, subset_list_maxmax, p_maxmax ) )
        g.intervene (j_maxmax, subset_list_maxmax)           # 혹시 아버지 -할아버지의 관계까 아닌지 찾는다
        for proband_clone_index in subset_list_maxmax:      # 4 -6,  4- 7
            g.addEdge(j_maxmax, proband_clone_index)

    print ("", file = output_file)
    print ("")
    for root in completed_parent:
        if g.findparent(root) == "None":   # root node일 때에만 돌자
            g.dfs(root, kwargs["PHYLOGENY_DIR"])
        print ("", file = output_file)
        print ("")


    output_file.close()
    return g