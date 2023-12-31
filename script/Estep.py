import numpy as np
import scipy
import copy
from scipy.special import expit
import random
import math

def phred_to_percentile (phred_score):
    return round ( 10 ** (-phred_score / 10) , 5)  #  Converts a Phred score to a percentile.


def expected_calculator (  i, j, k, mixture, df, input_containpos, **kwargs ):
    import re

    # if (kwargs["SEX"] == "M") & ( bool(re.search(r'X|Y', input_containpos.iloc[k]["pos"]))  == True  ) :
    #     #depth_calc, alt_calc = int(df[k][i]["depth"] ), int( round ( df[k][i]["depth"] * mixture[i][j]) )  
    #     depth_calc, alt_calc = int(df[k][i]["depth"] ), int( df[k][i]["depth"] * mixture[i][j]) 
    #     if mixture[i][j] > 1:
    #         alt_calc = int( df[k][i]["depth"] * 1 ) 
    # else:
        #depth_calc, alt_calc = int(df[k][i]["depth"] ), int( round ( df[k][i]["depth"] * mixture[i][j] * 0.5) )
    depth_calc, alt_calc = int(df[k][i]["depth"] ), int( df[k][i]["depth"] * mixture[i][j] * 0.5) 

    depth_obs, alt_obs = int(df[k][i]["depth"]), int(df[k][i]["alt"])

    if (kwargs["MAKEONE_STRICT"] in [1,2, 3] )  & (int(df[k][i]["depth"] < 100) ):
        depth_calc = depth_obs = 100
        alt_calc = int ( 50 * mixture[i][j] ) 
        alt_obs = int ( alt_obs * (100 / int(df[k][i]["depth"])) )

    a, b = alt_calc, depth_obs - alt_calc              # alt_expected, ref_expected

    return (depth_calc, alt_calc, depth_obs, alt_obs, a, b)

def equivalent_prior (alt_obs, i, j, mixture, TF, **kwargs):
    equivalent = []

    if TF in ["TN", "FN"]:
        TN_clone_index = np.where ( mixture[i] == 0 )[0]
        if TF == "TN":
            return kwargs["TN_PRIOR"] / len (TN_clone_index)
        elif TF == "FN":
            return ( 1 - kwargs ["TN_PRIOR"] / (kwargs ["NUM_CLONE"] - len(TN_clone_index)))
        
    elif TF in ["FP", "TP"]:
        if TF == "FP":
            return kwargs ["FP_PRIOR"]
        elif TF == "TP":
            FP_clone_index = np.where ( mixture[i] == 0 )[0]
            return ( 1 -  ( kwargs ["FP_PRIOR"] * len (FP_clone_index) ) )  / (kwargs["NUM_CLONE"] - len(FP_clone_index))



def calc_posterior(input_containpos, df,  np_vaf, np_BQ, step, k, **kwargs):
    import re
    global debug_k

    mixture = step.mixture

    max_posterior_allsample = max_likelihood_allsample = float("-inf")
    max_posterior_clone_candidate = []
    max_clone = -1

    posterior_allsample = np.zeros(kwargs["NUM_CLONE"], dtype="float64")            # 합치면 1 되게
    likelihood_allsample = np.zeros(kwargs["NUM_CLONE"], dtype="float64")          # 순수한 beta binomal만..

    check = 0

    if kwargs["DEBUG"] == True:
        if (k in debug_k):
            print ("\t\t[DEBUG] k = {}\tnp_vaf*2 = {}\tmixture = {}".format(k, np_vaf[k]*2,  "\t".join(str(np.round(row, 2)) for row in step.mixture)  )  )
    else:
        debug_k = []

    proportion = np.unique(step.membership, return_counts = True)[1] / kwargs["NUM_MUTATION"]   # clone 당 variant 개수


    # Denominator (분모)를 한번에 계산해놓자
    L_denominator = np.zeros ( (kwargs["NUM_BLOCK"]  ), dtype = "float")
    for i in range(kwargs["NUM_BLOCK"]):
        for j in range(kwargs["NUM_CLONE"]): 
            SEQ_ERROR = phred_to_percentile ( np_BQ[k][i] )
            depth_calc, alt_calc, depth_obs, alt_obs, a, b = expected_calculator (  i, j, k, mixture, df, input_containpos, **kwargs )

            if mixture[i][j] == 0:    # TN or FP   (Binomial)
                TF = "TN" if alt_obs == 0 else "FP"
            else:    # TN or FP   (Binomial)
                TF = "FN" if alt_obs == 0 else "TP"
            # Prior 계산
            equivalent = equivalent_prior (alt_obs, i, j, mixture, TF, **kwargs)

            try:
                if (alt_obs == 0) & (mixture [i][j] == 0):    # TN
                    #L_denominator[i] += 0
                    L_denominator[i] += (equivalent * scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = 0))
                elif (alt_obs != 0) & (mixture [i][j] == 0):    # FP
                    L_denominator[i] += (equivalent * scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs))          # 여기에 Prior를 뭘 곱해야 할까?
                else:   # FN or TN
                    if  kwargs["MODEL"] == "betabinomial":
                        L_denominator[i] += (equivalent * scipy.stats.betabinom.pmf( alt_obs, depth_obs, a+1, b+1 )  )
                    elif kwargs["MODEL"] == "binomial":
                        L_denominator[i] += (equivalent * scipy.stats.binom.pmf( alt_obs , int (depth_obs *  mixture[i][j]  ) , 0.5)  )
            except:
                L_denominator[i] += 0

        if(kwargs ["DEBUG"] == True ) :
            if  ( k in debug_k)  :            
                print ("\t\t\t\t\ti = {}\tL_denominator[i] = {}".format ( i, np.round (L_denominator, 7) ) )

        if L_denominator[i] <= 0 :  # 전혀 말도 안되는 경우
            print ("L_denominator[{}] = {}\tk = {}\talt_obs = {}\tdepth_obs = {}".format(i, L_denominator[i], k, alt_obs, depth_obs))
            L_denominator[i] = 1


    # Likelihood estimation & Assignment
    for j in range(kwargs["NUM_CLONE"]): 
        SEQ_ERROR = phred_to_percentile ( np_BQ[k][i] )

        for i in range(kwargs["NUM_BLOCK"]):
            depth_calc, alt_calc, depth_obs, alt_obs, a, b = expected_calculator (  i, j, k, mixture, df, input_containpos, **kwargs )

            if mixture[i][j] == 0:    # TN or FP   (Binomial)
                TF = "TN" if alt_obs == 0 else "FP"
                equivalent = equivalent_prior (alt_obs, i, j,  mixture, TF, **kwargs)

                try:
                    p_likelihood = math.log10( scipy.stats.binom.pmf(n = depth_obs, p = SEQ_ERROR, k = alt_obs) )
                except:
                    p_likelihood, p_numerator, posterior = -400, -400, -400

                if kwargs ["PRIOR"] == "equivalent":
                    p_numerator = math.log10( equivalent ) + p_likelihood
                    posterior = p_numerator - math.log10 (  L_denominator[i] ) 

                try:
                    if(kwargs ["DEBUG"] == True ) :
                        if  ( k in debug_k)  :            
                            np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
                            if TF == "TN":
                                print ( "\t\t\t\t\t\tj = {}(TN)\ti = {}\tmixture = {}\talt_obs = {}\tdepth_obs = {}\tSEQ_ERROR = {}\tPRIOR = {}\tlog (likelihood) = {}\tlog(numerator) = {}\tlog(posterior) = {}".format( j, i,  round(mixture[i][j],2) , alt_obs, depth_obs, SEQ_ERROR, round (equivalent, 2)  , round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round ( posterior, 2 ) ) )
                            else:
                                print ( "\t\t\t\t\t\tj = {}(FP)\ti = {}\tmixture = {}\talt_obs = {}\tdepth_obs = {}\tSEQ_ERROR = {}\tPRIOR = {}\tlog (likelihood) = {}\tlog(numerator) = {}\tlog(posterior) = {}".format( j, i,  round(mixture[i][j],2) , alt_obs, depth_obs, SEQ_ERROR, round (equivalent, 2) , round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round ( posterior, 2 ) ) )
                except:
                    p_likelihood, p_numerator, posterior = -400, -400, -400
                    if(kwargs ["DEBUG"] == True ) :
                        if  ( k in debug_k)  :            
                            np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
                            if TF == "TN":
                                print ( "\t\t\t\t\t\tj = {}(TN)\ti = {}\tmixture = {}\talt_obs = {}\tdepth_obs = {}\tSEQ_ERROR = {}\tPRIOR = {}\tlog (likelihood) = {}\tlog(numerator) = {}\tlog(posterior) = {}".format( j, i,  round(mixture[i][j],2) , alt_obs, depth_obs, SEQ_ERROR, round(equivalent, 2)  , round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round ( posterior, 2 ) ) )
                            else:
                                print ( "\t\t\t\t\t\tj = {}(FP)\ti = {}\tmixture = {}\talt_obs = {}\tdepth_obs = {}\tSEQ_ERROR = {}\tPRIOR = {}\tlog (likelihood) = {}\tlog(numerator) = {}\tlog(posterior) = {}".format( j, i,  round(mixture[i][j],2) , alt_obs, depth_obs, SEQ_ERROR, round(equivalent, 2) , round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round ( posterior, 2 ) ) )


            else: # FN or TP   (Beta binomial)
                TF = "FN" if alt_obs == 0 else "TP"
                equivalent = equivalent_prior (alt_obs, i, j, mixture, TF, **kwargs)

                try:
                    if  kwargs["MODEL"] == "betabinomial":  
                        p_likelihood = math.log10 ( scipy.stats.betabinom.pmf(alt_obs, depth_obs, a+1, b+1) ) 
                    elif kwargs["MODEL"] == "binomial":  
                        p_likelihood = math.log10 ( scipy.stats.binom.pmf( alt_obs , int (depth_obs *  mixture[i][j]  ) , 0.5) ) 

                    p_numerator = p_likelihood + math.log10( equivalent )
                    posterior = p_numerator - math.log10 (  L_denominator[i] ) 

                    if(kwargs ["DEBUG"] == True ) :
                        if  ( k in debug_k)  :            
                            np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
                            if  kwargs["MODEL"] == "betabinomial":  
                                print ( "\t\t\t\t\t\tj = {}({})\ti = {}\tmixture = {}\talt_obs = {}\talt_calc = {}\tdepth_obs = {}\tPRIOR = {}\tlog(likelihood) = {}\tlog (numerator) = {}\tlog(posterior) = {}".format( j, TF,  i,  round(mixture[i][j],2) , alt_obs, alt_calc, depth_obs, round(equivalent, 2), round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round(posterior, 2 ) ) )
                            elif kwargs["MODEL"] == "binomial":  
                                print ( "\t\t\t\t\t\tj = {}({})\ti = {}\tmixture = {}\talt_obs = {}\talt_calc = {}\tdepth_obs = {}\tPRIOR = {}\tlog(likelihood) = {}\tllog (numerator) = {}\tlog(posterior) = {}".format( j, TF, i,  round(mixture[i][j],2) , alt_obs, alt_calc, depth_obs, round(equivalent, 2), round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round(posterior, 2 ) ) )
                except:
                    p_likelihood, p_numerator, posterior = -400, -400, -400
                    if(kwargs ["DEBUG"] == True ) :
                        if  ( k in debug_k)  :            
                            np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
                            if  kwargs["MODEL"] == "betabinomial":  
                                print ( "\t\t\t\t\t\tj = {}({})\ti = {}\tmixture = {}\talt_obs = {}\talt_calc = {}\tdepth_obs = {}\tPRIOR = {}\tlog(likelihood) = {}\tlog (numerator) = {}\tlog(posterior) = {}".format( j, TF,  i,  round(mixture[i][j],2) , alt_obs, alt_calc, depth_obs, round(equivalent, 2), round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round(posterior, 2 ) ) )
                            elif kwargs["MODEL"] == "binomial":  
                                print ( "\t\t\t\t\t\tj = {}({})\ti = {}\tmixture = {}\talt_obs = {}\talt_calc = {}\tdepth_obs = {}\tPRIOR = {}\tlog(likelihood) = {}\tllog (numerator) = {}\tlog(posterior) = {}".format( j, TF, i,  round(mixture[i][j],2) , alt_obs, alt_calc, depth_obs, round(equivalent, 2), round ( p_likelihood, 2 ), round ( p_numerator, 2 ), round(posterior, 2 ) ) )




            # 각 sample마다 더해준다  (log니까)
            posterior_allsample[j] += posterior
            likelihood_allsample[j] += p_numerator

        if posterior_allsample[j] > max_posterior_allsample:
            max_posterior_clone_candidate = [j]
            max_posterior_allsample = posterior_allsample[ j ]
            max_likelihood_allsample = likelihood_allsample[ j ]
        elif posterior_allsample[j] == max_posterior_allsample :     # 동점이 나오는 경우도 있어서
            max_posterior_clone_candidate.append (j)
    
    try:
        max_clone = random.choice( max_posterior_clone_candidate )
    except:
        print ( "max clone choice가 이상하다", k, max_posterior_clone_candidate, posterior_allsample )
    max_posterior_allsample = posterior_allsample[ max_clone ]
    max_likelihood_allsample = likelihood_allsample[ max_clone ]

    posterior_allsample_normalized = np.power (10, posterior_allsample) 
    likelihood_allsample_normalized =  np.power (10, likelihood_allsample) 
    posterior_allsample_normalized = posterior_allsample_normalized / np.sum (posterior_allsample_normalized)
    likelihood_allsample_normalized = likelihood_allsample_normalized / np.sum (likelihood_allsample_normalized)
    max_posterior_allsample_normalized = posterior_allsample_normalized [ max_clone ]
    max_likelihood_allsample_normalized = likelihood_allsample_normalized [ max_clone ]

    # if max_posterior_allsample_normalized <= 0:
    #     print (max_posterior_allsample_normalized, k, max_clone)

    

    if(kwargs ["DEBUG"] == True ) :
        if  ( k in debug_k)  :
            np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
            print ( "\t\t\t\t\t❯❯❯ max_clone = {}\tlog(likelihood_all) = {}\tlog(posterior_all) = {}\tposterior_allsample_normalized = {}".format( max_clone, np.round (likelihood_allsample, 2) , np.round (posterior_allsample, 2) , np.round ( posterior_allsample_normalized, 2)   ))


    if kwargs["OPTION"] in ["Hard", "hard"]:
        return list(posterior_allsample), max_posterior_allsample, max_likelihood_allsample, posterior_allsample_normalized, likelihood_allsample_normalized, max_posterior_allsample_normalized, max_likelihood_allsample_normalized, max_clone

    elif kwargs["OPTION"] in ["Soft", "soft"]:
        weight_posterior = np.power (10, posterior_allsample)  
        soft_posterior_allsample = round( np.average(posterior_allsample, weights = weight_posterior), 3)       # Posterior in Soft clustering
        soft_posterior_allsample_normalized = np.average(posterior_allsample_normalized, weights = weight_posterior)
        #soft_posterior_allsample_normalized = posterior_allsample_normalized [ max_clone ]

        weight_likelihood = np.power (10, likelihood_allsample)  
        soft_likelihood_allsample = round( np.average(likelihood_allsample, weights = weight_likelihood), 3)       # Likelihood in Soft clustering
        soft_likelihood_allsample_normalized = np.average(likelihood_allsample_normalized, weights = weight_likelihood)
        #soft_likelihood_allsample_normalized = likelihood_allsample_normalized [ max_clone ]
        
        if(kwargs ["DEBUG"] == True ) :
            if  ( k in debug_k)  :
                np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
                

        return list(posterior_allsample), soft_posterior_allsample, soft_likelihood_allsample, posterior_allsample_normalized, likelihood_allsample_normalized, soft_posterior_allsample_normalized, soft_likelihood_allsample_normalized,  max_clone      # mixture 구할때에는 posterior를 이용하고, termination 조건 구할 때에는 likelihood를 이용하자



def main (input_containpos, df, np_vaf, np_BQ, step, **kwargs):
    import math
    global debug_k
    total_posterior_allsample = total_likelihood_allsample = 0
    total_posterior_allsample_normalized = total_likelihood_allsample_normalized = 0

    if kwargs["DEBUG"] == True:
        #debug_k = random.choice ( np.where(  ( (np_vaf[:, 0]  > 0.3 ) &  (np_vaf[:, 0] < 0.14 )) )  [0] )
        # print ( "debug_k = {}".format(debug_k))
        debug_k = [410]
    else:
        debug_k = -1

    kwargs ["PRIOR"] = "equivalent"   # equivalent, centroid, proportion

    #max_posterior_allsample_list, max_likelihood_allsample_list = [], []
    for k in range(kwargs["NUM_MUTATION"]):
        step.membership_p[k], max_posterior_allsample, max_likelihood_allsample, posterior_allsample_normalized, likelihood_allsample_normalized, max_posterior_allsample_normalized, max_likelihood_allsample_normalized, step.membership[k] = calc_posterior(input_containpos, df,  np_vaf, np_BQ, step, k, **kwargs)
        if math.isnan(max_posterior_allsample) == False:
            total_posterior_allsample += max_posterior_allsample
        if math.isnan(max_likelihood_allsample) == False:    
            total_likelihood_allsample += max_likelihood_allsample
        if math.isnan(max_posterior_allsample_normalized) == False:           
            total_posterior_allsample_normalized += math.log10 ( max_posterior_allsample_normalized )
        if math.isnan(max_likelihood_allsample_normalized) == False:           
            total_likelihood_allsample_normalized += math.log10 ( max_likelihood_allsample_normalized )


        if kwargs["OPTION"] in ["Soft", "soft"]:
            if step.membership[k] == kwargs["NUM_CLONE_ITER"]:  # FP clone
                step.membership_p_normalize[k] = np.zeros( kwargs["NUM_CLONE"], dtype="float64")  # Set  1 (FP_index) 0 0 0 0 0    
                step.membership_p_normalize[k][ kwargs["NUM_CLONE_ITER"] ] = 1
            else:
                step.membership_p_normalize[k] = np.round ( posterior_allsample_normalized, 2)
                step.membership_p_normalize[k][ kwargs["NUM_CLONE_ITER"] ] = 0

        # if k % 100 == 0:
        #     print ("k = {}, vaf * 2 = {}\tp_a_n = {}".format ( k, round ( int(df[k][0]["alt"] ) / int(df[k][0]["depth"]) * 2, 2), np.round ( posterior_allsample_normalized, 2)  ) )


    step.posterior, step.posterior_record[kwargs["STEP"]] = total_posterior_allsample, total_posterior_allsample
    step.likelihood, step.likelihood_record[kwargs["STEP"]] = total_likelihood_allsample, total_likelihood_allsample

    step.posterior_normalized =  step.posterior_normalized_record[kwargs["STEP"]] =   total_posterior_allsample_normalized 
    step.likelihood_normalized =  step.likelihood_normalized_record[kwargs["STEP"]] =  total_likelihood_allsample_normalized 


    if  step.fp_index in set(step.membership):   # FP를 선택한 variant가 하나라도 있어야 비로소 includefp = True로 변경
        step.includefp = True
        step.fp_member_index = list(np.where(step.membership == step.fp_index)[0])
    else:
        step.includefp = False
        step.fp_member_index = []


    if (kwargs["VERBOSE"] >= 1):
        np.set_printoptions(suppress=True)   # Scientific expression이 싫어요
        #print ("\t\t\tEstep.py : set(step.membership) = {}\tcounts = {}\tincludefp = {}\tstep.likelihood = {}\tstep.likelihood_normalized = {}\tstep.posterior = {}\tstep.posterior_normalized = {}".format ( set(step.membership), np.unique(step.membership  , return_counts=True)[1], step.includefp, round(step.likelihood), round( step.likelihood_normalized, 2), round(step.posterior),  round( step.posterior_normalized, 2) ) )
        print ("\t\t\tEstep.py : set(step.membership) = {}\tcounts = {}\tincludefp = {}\tstep.likelihood = {}\tstep.posterior_normalized = {}".format ( set(step.membership), np.unique(step.membership  , return_counts=True)[1], step.includefp, round(step.likelihood), round( step.posterior_normalized, 2) ) )


    # membership_p : extermely low value if the variant is  fp
    for k in range(kwargs["NUM_MUTATION"]):
        if step.membership[k] == step.fp_index:
            step.membership_p[k] = [-999] * len(step.membership_p[k])

    return step


