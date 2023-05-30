import numpy as np
import scipy
import copy
from scipy.special import expit
import random
import math

def phred_to_percentile (phred_score):
    return round ( 10 ** (-phred_score / 10) , 3)  #  Converts a Phred score to a percentile.



# 해당 mutation이 각 cluster에 속할 확률 (soft : prob[], hard : max_prob)
def calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs):
    mixture = step.mixture

    max_prob = float("-inf")
    max_clone = -1
    #SEQ_ERROR1 = 0.01     # 한 read가 FP가 생길 확률
    #SEQ_ERROR2 = 0.995   # 한 read가 TN일 확률
    SEQ_ERROR2 = kwargs["TN_CONFIDENTIALITY"]

    prob = np.zeros(kwargs["NUM_CLONE"], dtype="float64")

    check = 0

    for j in range(kwargs["NUM_CLONE"]):      # 몇 번째 clone으로 편입되는 게 가장 좋으려나?
        
        if (j == step.fp_index):    # FP로 편입되는걸 고려할 떄에는 SEQ_ERROR를 이용한 probability
            # if ( int (df[k][1]["alt"] == 0)  ) & (kwargs["STEP"] == 1) & (k < 50): 
            #     check = 1
                
            for i in range(kwargs["NUM_BLOCK"]):
                SEQ_ERROR1 = phred_to_percentile ( np_BQ[k][i] )
                depth_calc, alt_calc = int(df[k][i]["depth"] ), int(df[k][i]["depth"] * mixture[i][j] * 0.5)
                depth_obs, alt_obs = int(df[k][i]["depth"]), int(df[k][i]["alt"])

                # 1번 (FP)와 4번 (TN)에서의 SEQ_ERROR 값은 구분하는 것이 좋겠다
                SEQ_ERROR = (1 - SEQ_ERROR2) if alt_obs == 0 else SEQ_ERROR1
                try:
                    p = math.log10(scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs))
                    prob[j] =  prob[j] + p
                except:
                    prob[j] = -399
                    
                # if check == 1: 
                #     print ("\t\tk = {} (black 가정) (i = {}): {},{}\tw = {}\tp = {}".format( k,  i,   df[k][i]["alt"], df[k][i]["depth"],SEQ_ERROR, round (p, 2) ))

        else:  # 일반 clone으로 편입되는 걸 계산
            # if ( int (df[k][1]["alt"] == 0)  )  & (kwargs["STEP"] == 1) & (k < 50) & ( step.mixture[0][j] < 0.2) & ( step.mixture[1][j] < 0.2)  : 
            #     check = 1
                
            for i in range(kwargs["NUM_BLOCK"]):
                SEQ_ERROR1 = phred_to_percentile ( np_BQ[k][i] )
                depth_calc, alt_calc = int(df[k][i]["depth"] ), int(df[k][i]["depth"] * mixture[i][j] * 0.5)
                depth_obs, alt_obs = int(df[k][i]["depth"]), int(df[k][i]["alt"])

                # Beta binomial distribution
                a = alt_calc              # alt_expected
                b = depth_obs - a            # ref_expected

                # 굳이 space에 있는 mutation과 axis에 있는 mutation을 구별해야 하나? (SEQ_ERROR1과 SEQ_ERROR2를 구분하는 게 좋겠다)
                SEQ_ERROR = (1 - SEQ_ERROR2) if alt_obs == 0 else SEQ_ERROR1
                # 일단 FP가 아니여야 하고 (binomial) + 이 cluster에 포함되어야 하고 (beta binomial)
                try:
                    p1 = math.log10(1 - scipy.stats.binom.pmf(n=depth_obs, p=SEQ_ERROR, k=alt_obs))
                    p2 = math.log10(scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1))
                    prob[j] = prob[j] + p1 + p2
                except:
                    prob[j] = prob[j] - 400

                # if check == 1:
                #     print ("\t\tk = {} (j = {}, i = {}): {},{}\tw = {}\tp1 = {}\tp2 = {}".format( k, j, i, df[k][i]["alt"], df[k][i]["depth"],SEQ_ERROR,  round(p1, 2) , round( p2, 2)  ))

        if prob[j] > max_prob:
            max_prob = prob[j]
            max_prob_clone_candidate = [j]
        elif prob[j] == max_prob:
            max_prob_clone_candidate.append(j)

    if ( int (df[k][i]["alt"] == 0)  ) & (kwargs["STEP"] == 1) & (k < 50) & (check == 1): 
        print ("\t\t\t  k = {}, prob = {}".format(k, prob))

    max_clone = random.choice(max_prob_clone_candidate)




    if kwargs["OPTION"] in ["Hard", "hard"]:
        # if np.argmax (prob) == step.fp_index:
        #     print ("{}번째 mutation : fp ({})로 배정되고 싶어함   ({p =})   (다른 cluster : {})".format(k, step.fp_index, round (np.max(prob), 2), np.round (prob, 2) ))
        #print("mutation_{0} : {1}번 clone이 가장 적절 (log p = {2})".format(k, max_clone, max_prob))
        return list(prob), max_prob, max_clone

    elif kwargs["OPTION"] in ["Soft", "soft"]:
        weight = np.zeros(len(list(prob)), dtype="float")
        for j in range(len(list(prob))):
            weight[j] = math.pow(10, prob[j])

        new_likelihood = round(np.average(prob, weights=weight), 4)       # Soft clustering에서의 likelihood계산법
        #print ("Total likelihood : {}\tSoft weighted likelihood : {}".format(max_prob, new_likelihood))

        return list(prob), new_likelihood, max_clone  # 어쨌든 max_clone 을 배정하기는 한다



def main (df, np_vaf, np_BQ, step, **kwargs):
    total_prob = 0

    # (230122 수정)  Soft step인데, fp가 있긴 있을 때, 얘를 rescue 해줄 수 있는 지 한번 살펴보자
    if (step.fp_index != -1) & (kwargs["OPTION"] in ["Soft", "soft"]):
        temp_fp_index = step.fp_index

        #1. 원래대로 fp를 fp로 칠 경우
        step.fp_index, total_prob1 = temp_fp_index, 0
        temp1_membership = np.zeros ( kwargs["NUM_MUTATION"]  ,dtype = "int")
        temp1_membership_p = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"])  ,dtype = "float")
        for k in range(kwargs["NUM_MUTATION"]):
            temp1_membership_p[k], max_prob, temp1_membership[k] = calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs)
            total_prob1 = total_prob1 + max_prob

        #2. fp를 독립된 clone으로 볼 경우
        step.fp_index, total_prob2 = -1, 0
        temp2_membership = np.zeros ( kwargs["NUM_MUTATION"]  ,dtype = "int")
        temp2_membership_p = np.zeros ( (kwargs["NUM_MUTATION"], kwargs["NUM_CLONE"])  ,dtype = "float")
        for k in range(kwargs["NUM_MUTATION"]):
            temp2_membership_p[k], max_prob, temp2_membership[k] = calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs)
            total_prob2 = total_prob2 + max_prob

        #1. 1번 경우가 더 나을 경우
        if total_prob1  > total_prob2:
            if kwargs["VERBOSE"] >= 2:
                print ("1번 경우가 더 낫다 :  원래대로 fp를 fp로 칠 경우 = {} > fp를 독립된 clone으로 rescue할 경우 = {}".format( round(total_prob1, 2) , round (total_prob2, 2) ))
            step.fp_index = temp_fp_index
            step.likelihood = total_prob1
            step.membership = copy.deepcopy  ( temp1_membership  )
            step.membership_p = copy.deepcopy  ( temp1_membership_p  )
        else:
            if kwargs["VERBOSE"] >= 2:
                print ("2번 경우가 더 낫다 :  원래대로 fp를 fp로 칠 경우 = {} < fp를 독립된 clone으로 rescue할 경우 = {}".format( round(total_prob1, 2) , round (total_prob2, 2) ))
            step.makeone_index = sorted ( step.makeone_index + [temp_fp_index] )       # rescue해서 makeone_index로 소생시켜줘야지
            step.fp_index = -1
            step.likelihood = total_prob2
            step.membership = copy.deepcopy  ( temp2_membership  )
            step.membership_p = copy.deepcopy  ( temp2_membership_p  )

    else:  #  (대부분의 경우))
        # k번째 mutation은 membership을 누구에게 주는 게 좋을까?
        total_prob = 0
        FPorCluster3_prob = 0
        for k in range(kwargs["NUM_MUTATION"]):
            step.membership_p[k], max_prob, step.membership[k] = calc_likelihood(df,  np_vaf, np_BQ, step, k, **kwargs)
            total_prob = total_prob + max_prob
        step.likelihood = total_prob
        step.likelihood_record[kwargs["STEP"]] = total_prob

        # if (kwargs["STEP"] == 1):
        #     print ("\t\tFP prob 합 = {}".format (round (FPorCluster3_prob, 2))) if step.fp_index != -1 else print ("\t\t구석탱이 cluster 합 = {}".format (round (FPorCluster3_prob, 2)))



    # membership_p의 경우 fp로 배정된다면 굉장히 낮은 값을 줘야
    for k in range(kwargs["NUM_MUTATION"]):
        if step.membership[k] == step.fp_index:
            step.membership_p[k] = [-999] * len(step.membership_p[k])

    return step
