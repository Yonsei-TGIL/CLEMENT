import numpy as np
import scipy
from scipy.special import expit
import random, math


def TN_prior_cal(x):
    return (1 - expit( 30*x - 5)) * 0.7
    #return (1 - expit( 50*x - 5)) * 0.5

def situation1 (step, depth_calc, alt_calc, depth_obs, alt_obs, a, b, prob, j, np_vaf_nonzero, SEQ_ERROR):   
    # In case of mutation in AXIS  →  decide whether allocate to FP or True clone
    FN_prior = 1
    for i2 in range (len (np_vaf_nonzero)) :
        FN_prior = FN_prior * (1 - TN_prior_cal ( np_vaf_nonzero [i2] * 2  ))

    prob[j] = prob[j] + FN_prior * math.log10(scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1))

    return prob




def calc_likelihood(df,  np_vaf, step, k, **kwargs):       # Probability of each mutation belong to each cluster (soft : prob[], hard : max_prob)
    mixture = step.mixture

    max_prob = float("-inf")
    max_clone = -1
    SEQ_ERROR = 0.015

    prob = np.zeros(kwargs["NUM_CLONE"] , dtype="float64")       

    for j in range(kwargs["NUM_CLONE"]):
        np_vaf_nonzero = []
        for i in range(kwargs["NUM_BLOCK"]):
            if int (df[k][i]["alt"]) != 0:
                np_vaf_nonzero.append ( np_vaf[k][i] )            #  if, NUM_BLOCK == 3 & df[alt] =  (0, 0, 2)   →  save np_vaf of  block3

        if  ( j  == step.fp_index ) :    # Situation2 : Probability should be structed based on SEQ_ERROR conditional probability
            FP_prob = math.pow(SEQ_ERROR , len(np_vaf_nonzero))
            for i2 in range ( len (np_vaf_nonzero) ) :
                FP_prob = FP_prob * TN_prior_cal ( np_vaf_nonzero [i2]  * 2 )
            try:
                prob[j] = math.log10(FP_prob)
            except:
                prob[j] = prob[j] - 400
    


        else:
            for i in range(kwargs["NUM_BLOCK"]):
                depth_calc, alt_calc = int(df[k][i]["depth"] * mixture[i][j]), int(df[k][i]["depth"] * mixture[i][j] * 0.5)
                depth_obs, alt_obs = int(df[k][i]["depth"]), int(df[k][i]["alt"])

                # Beta binomial distribution
                a = df[k][i]["depth"] * mixture[i][j] * 0.5              # alt_expected
                b = depth_obs - a            # ref_expected
                try:
                    if alt_obs != 0:
                        prob[j] = prob[j] + math.log10(scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1)) + math.log10 ( 1 -  TN_prior_cal ( (alt_obs / depth_obs) * 2) ) 
                    elif alt_obs == 0:     # In case of mutation in AXIS  →  decide whether allocate to FP or True clone
                        prob = situation1 ( step, depth_calc, alt_calc, depth_obs, alt_obs, a, b, prob, j, np_vaf_nonzero, SEQ_ERROR )

                except:
                    prob[j] = prob[j] - 400

            # print("{0}th mutation : {1}th clone, {2}th sample → alt_expected : {3},  alt_observed : {4}, depth_observed : {5}, likelihood : {6}"
            #           .format(k, j, i, round(a, 1), alt_obs, depth_obs,  scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1)))

        if prob[j] > max_prob:
            max_prob = prob[j]
            max_prob_clone_candidate = [j]
        elif prob[j] == max_prob:
            max_prob_clone_candidate.append(j)

    # if np.argmax (prob) == step.fp_index:
    #     print ("{}th mutation : fp ({}) is the most likely than any other clones   ({})".format(k, step.fp_index, prob))

    max_clone = random.choice(max_prob_clone_candidate)

    # print("mutation_{0} : {1}th clone is the most appropriate (log p = {2})\n".format(k, max_clone, max_prob))
    if kwargs["OPTION"] in ["Hard", "hard", "Outlier", "outlier", "Start", "start"]:
        return list(prob), max_prob, max_clone

    elif kwargs["OPTION"] in ["Soft", "soft"]:
        weight = np.zeros( len(list(prob)), dtype="float")
        for j in range ( len(list(prob)) ):
            weight[j] = math.pow(10, prob[j])
    
        new_likelihood = round( np.average(prob , weights=weight), 4)
        #print ("Total likelihood : {}\tSoft weighted likelihood : {}".format(max_prob, new_likelihood))

        return list(prob), new_likelihood, max_clone



def main(df, np_vaf, step, **kwargs): 

    total_prob = 0
    for k in range(kwargs["NUM_MUTATION"]):       # Allocate "k"th mutation to which clone?
        temp_membership_p, max_prob, step.membership[k] = calc_likelihood(df,  np_vaf, step, k, **kwargs)
        total_prob = total_prob + max_prob

        if step.membership[k] == step.fp_index:
            step.membership_p[k] = [-999] * len(step.membership_p[k])
        else:  # Most of the case
            step.membership_p[k] = temp_membership_p


    step.likelihood = total_prob
    step.likelihood_record [kwargs["STEP"]] = total_prob

     #print("\tTotal likelihood : {0}".format(round(total_prob, 2)))

    return step
