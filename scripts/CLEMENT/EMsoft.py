import os
import numpy as np
import Estep, Mstep, Bunch, miscellaneous

import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)
warnings.filterwarnings("ignore")


def iszerocolumn (step, **kwargs):
    for j in range ( kwargs ["NUM_CLONE"] ) :
        if  np.array_equal (step.mixture[: , j],  np.zeros ( kwargs["NUM_BLOCK"]) )  == True :
            #print ("\t\tNull clone found.\t{}".format (  step.mixture  ))
            return True
    return False



def main (cluster, df, np_vaf, mixture_kmeans, **kwargs):
    NUM_CLONE = kwargs ["NUM_CLONE"]
    NUM_BLOCK, kwargs["NUM_BLOCK"]= len(df[0]), len(df[0])
    NUM_MUTATION, kwargs["NUM_MUTATION"]  = kwargs["RANDOM_PICK"], kwargs["RANDOM_PICK"]
    kwargs["STEP_NO"] = 30
    kwargs["OPTION"] = "soft"

    cluster = Bunch.Bunch2(**kwargs)
    
    print("\n\n\tNow Soft clustering trial begun  \n")
    trial = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["TRIAL_NO"])
    
    trial_index, failure_num = 0, 0
    while trial_index < kwargs["TRIAL_NO"]:
        kwargs["TRIAL"] = trial_index
        print("\t#{0}th trial".format(trial_index))

        step = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["STEP_NO"])
        step.mixture = miscellaneous.set_initial_parameter(df, mixture_kmeans, **kwargs)


        for step_index in range(0, kwargs["STEP_NO"]):
            kwargs["STEP"], kwargs["STEP_TOTAL"] = step_index, step_index
            kwargs["OPTION"] = "soft"
            
            
            step = Estep.main(df, np_vaf, step, **kwargs)                 # E step: Setting new membership under new mixture
            step = Mstep.main(df, np_vaf, step, "Soft", **kwargs)    # M step: Setting new mixture under new membrship

            step.acc(step.mixture, step.membership, step.likelihood, step.membership_p, step.membership_p_normalize, step.makeone_index, step.fp_index, step_index, step.outlier_index, step.includeoutlier, kwargs["STEP"], kwargs["STEP"])       # step_index,  max_step_index : not important  in this line
            print ("\t\t{}th step\tFP_index : {}\tChild_index : {}\tstep.likelihood : {}".format(step_index, step.fp_index, step.makeone_index , round(step.likelihood, 1)))

            if (miscellaneous.GoStop(step, **kwargs) == "Stop") | ( iszerocolumn (step, **kwargs) == True):
                i =  step.find_max_likelihood(3, kwargs["STEP"]) 
                trial.acc ( step.mixture_record [i],  step.membership_record [i], step.likelihood_record [i], step.membership_p_record [i], step.membership_p_normalize_record [i], step.makeone_index_record[i], step.fp_index_record[i],  step_index + 1, step.outlier_index_record[i], step.includeoutlier_record[i], i, kwargs["TRIAL"] )
                trial_index = trial_index + 1
                failure_num = 0
                break

        if failure_num >= 2:  # Give up and advance to next trial
            trial.acc ( step.mixture,  np.random.randint (0, NUM_CLONE - 1, NUM_MUTATION)  , float("-inf"), step.membership_p, step.membership_p_normalize, step.makeone_index, step.fp_index, step_index + 1, step.outlier_index, step.includeoutlier, 0, kwargs["TRIAL"] )
            failure_num = 0
            trial_index = trial_index + 1

        
    i =  trial.find_max_likelihood(0, kwargs["TRIAL_NO"]) 

    if trial.likelihood_record[i] <= -9999999:
        print ("\n본격적으로 TRIAL_NO까지 돌렸는데도 모두 망해서 soft clustering의 결과값 없음 (trial.likelihood_record : {})".format(trial.likelihood_record) )
        cluster.acc ( trial.mixture_record [0], trial.membership_record [0], float("-inf"), trial.membership_p_record [0], trial.membership_p_normalize_record [0], trial.stepindex_record [0], 0, trial.max_step_index_record [0], trial.makeone_index_record[0], trial.fp_index_record[0], trial.includeoutlier_record[0], trial.outlier_index_record [0], **kwargs )  
    else:
        print ("\nWe chose {}th trial: (trial.likelihood_record : {})".format(i, trial.likelihood_record) )
        print ("\t\tsoft max : {} →  clone{}.{}-{}(soft).{} to candidate.{}".format ( trial.max_step_index_record[i],  NUM_CLONE, i , trial.max_step_index_record[i], kwargs["IMAGE_FORMAT"], kwargs["IMAGE_FORMAT"] ) )

        os.system ("cp " + kwargs["MYEM_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE"]) + "." + str( i ) + "-"  + str(  trial.max_step_index_record [i]  ) + "\(soft\).jpg" + " " + 
        kwargs["MYEM_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE"]) + ".\(soft\).jpg"  ) 
        
        cluster.acc ( trial.mixture_record [i], trial.membership_record [i], trial.likelihood_record [i], trial.membership_p_record [i], trial.membership_p_normalize_record [i], trial.stepindex_record [i], i, trial.max_step_index_record [i], trial.makeone_index_record[i], trial.fp_index_record[i], trial.includeoutlier_record[i], trial.outlier_index_record [i], **kwargs )  

    return cluster
