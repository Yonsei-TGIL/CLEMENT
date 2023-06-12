import numpy as np
import scipy
import copy
from scipy.special import expit
import random
import math

def phred_to_percentile (phred_score):
    return round ( 10 ** (-phred_score / 10) , 3)  #  Converts a Phred score to a percentile.


def calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs):
    mixture = step.mixture

    max_prob = float("-inf")
    max_clone = -1

    SEQ_ERROR2 = kwargs["TN_CONFIDENTIALITY"]

    prob = np.zeros(kwargs["NUM_CLONE"], dtype="float64")

    check = 0

    for j in range(kwargs["NUM_CLONE"]): 
        
        if (j == step.fp_index):   # Calculating the probability if this variants is false positive
            for i in range(kwargs["NUM_BLOCK"]):
                SEQ_ERROR1 = phred_to_percentile ( np_BQ[k][i] )
                depth_calc, alt_calc = int(df[k][i]["depth"] ), int(df[k][i]["depth"] * mixture[i][j] * 0.5)
                depth_obs, alt_obs = int(df[k][i]["depth"]), int(df[k][i]["alt"])

            
                SEQ_ERROR = (1 - SEQ_ERROR2) if alt_obs == 0 else SEQ_ERROR1
                try:
                    p = math.log10(scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs))
                    prob[j] =  prob[j] + p
                except:
                    prob[j] = -399
                    

        else:  #  # Calculating the probability if this variants is in other clusters
            for i in range(kwargs["NUM_BLOCK"]):
                SEQ_ERROR1 = phred_to_percentile ( np_BQ[k][i] )
                depth_calc, alt_calc = int(df[k][i]["depth"] ), int(df[k][i]["depth"] * mixture[i][j] * 0.5)
                depth_obs, alt_obs = int(df[k][i]["depth"]), int(df[k][i]["alt"])

                # Beta binomial distribution
                a = alt_calc              # alt_expected
                b = depth_obs - a            # ref_expected

                
                SEQ_ERROR = (1 - SEQ_ERROR2) if alt_obs == 0 else SEQ_ERROR1
                
                try:
                    p1 = math.log10(1 - scipy.stats.binom.pmf(n=depth_obs, p=SEQ_ERROR, k=alt_obs))
                    p2 = math.log10(scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1))
                    prob[j] = prob[j] + p1 + p2
                except:
                    prob[j] = prob[j] - 400


        if prob[j] > max_prob:
            max_prob = prob[j]
            max_prob_clone_candidate = [j]
        elif prob[j] == max_prob:
            max_prob_clone_candidate.append(j)

    if ( int (df[k][i]["alt"] == 0)  ) & (kwargs["STEP"] == 1) & (k < 50) & (check == 1): 
        print ("\t\t\t  k = {}, prob = {}".format(k, prob))

    max_clone = random.choice(max_prob_clone_candidate)




    if kwargs["OPTION"] in ["Hard", "hard"]:
        return list(prob), max_prob, max_clone

    elif kwargs["OPTION"] in ["Soft", "soft"]:
        weight = np.zeros(len(list(prob)), dtype="float")
        for j in range(len(list(prob))):
            weight[j] = math.pow(10, prob[j])

        new_likelihood = round(np.average(prob, weights=weight), 4)       # Likelihood in Soft clustering
        #print ("Total likelihood : {}\tSoft weighted likelihood : {}".format(max_prob, new_likelihood))

        return list(prob), new_likelihood, max_clone 



def main (df, np_vaf, np_BQ, step, **kwargs):
    total_prob = 0

    if (step.fp_index != -1) & (kwargs["OPTION"] in ["Soft", "soft"]):
        temp_fp_index = step.fp_index

        #1. FP is FP
        step.fp_index, total_prob1 = temp_fp_index, 0
        temp1_membership = np.zeros ( kwargs["NUM_MUTATION"]  ,dtype = "int")
        temp1_membership_p = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"])  ,dtype = "float")
        for k in range(kwargs["NUM_MUTATION"]):
            temp1_membership_p[k], max_prob, temp1_membership[k] = calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs)
            total_prob1 = total_prob1 + max_prob

        #2. FP is indepenent clone
        step.fp_index, total_prob2 = -1, 0
        temp2_membership = np.zeros ( kwargs["NUM_MUTATION"]  ,dtype = "int")
        temp2_membership_p = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"])  ,dtype = "float")
        for k in range(kwargs["NUM_MUTATION"]):
            temp2_membership_p[k], max_prob, temp2_membership[k] = calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs)
            total_prob2 = total_prob2 + max_prob


        if total_prob1  > total_prob2:
            if kwargs["VERBOSE"] >= 2:
                print ("#1 is favorable :  #1 = {} > #2 = {}".format( round(total_prob1, 2) , round (total_prob2, 2) ))
            step.fp_index = temp_fp_index
            step.likelihood = total_prob1
            step.membership = copy.deepcopy  ( temp1_membership  )
            step.membership_p = copy.deepcopy  ( temp1_membership_p  )
        else:
            if kwargs["VERBOSE"] >= 2:
                print ("#2 is favorable :  #1 = {} < #2 = {}".format( round(total_prob1, 2) , round (total_prob2, 2) ))
            step.makeone_index = sorted ( step.makeone_index + [temp_fp_index] )       # rescue and include it in makeone_index
            step.fp_index = -1
            step.likelihood = total_prob2
            step.membership = copy.deepcopy  ( temp2_membership  )
            step.membership_p = copy.deepcopy  ( temp2_membership_p  )

    else:  # Most of the cases
        total_prob = 0
        FPorCluster3_prob = 0
        for k in range(kwargs["NUM_MUTATION"]):
            step.membership_p[k], max_prob, step.membership[k] = calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs)
            total_prob = total_prob + max_prob
        step.likelihood = total_prob
        step.likelihood_record[kwargs["STEP"]] = total_prob




    # membership_p : extermely low value if the variant is  fp
    for k in range(kwargs["NUM_MUTATION"]):
        if step.membership[k] == step.fp_index:
            step.membership_p[k] = [-999] * len(step.membership_p[k])

    return step
