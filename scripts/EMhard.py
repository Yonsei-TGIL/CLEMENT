import os
import numpy as np
import  Estep, Mstep, Bunch, miscellaneous

import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)
warnings.filterwarnings("ignore")




def main ( df, np_vaf, mixture_kmeans, **kwargs):
    NUM_BLOCK, kwargs["NUM_BLOCK"]= len(df[0]), len(df[0])
    NUM_MUTATION =  kwargs["RANDOM_PICK"]
    kwargs["STEP_NO"] = 30

    cluster = Bunch.Bunch2(**kwargs)

    for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        kwargs["NUM_CLONE"], kwargs["NUM_CLONE_NOMINAL"] = NUM_CLONE, NUM_CLONE
        print("\n#########################################################################\n\nNUM_CLONE = {0}".format(NUM_CLONE))
        trial = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["TRIAL_NO"])
        
        trial_index, failure_num = 0, 0
        while trial_index < kwargs["TRIAL_NO"]:
            kwargs["TRIAL"] = trial_index
            if kwargs ["VERBOSE"] >= 1:
                print("\t#{0}th trial".format(trial_index))

            step = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["STEP_NO"])
            step.mixture = miscellaneous.set_initial_parameter(df, mixture_kmeans, **kwargs)


            for step_index in range(0, kwargs["STEP_NO"]):
                kwargs["STEP"], kwargs["STEP_TOTAL"] = step_index, step_index
                kwargs["OPTION"] = "hard"
                #print("\t\t#{}th step :  fp_index {}\tmakeone_index {}".format(step_index, step.fp_index, step.makeone_index))

                step = Estep.main(df, np_vaf, step, **kwargs)                   # E step: Setting new membership under new mixture
             
                if ( len ( set (step.membership) ) < NUM_CLONE ) |  ( np.min( np.unique(step.membership, return_counts=True)[1] ) < kwargs["MIN_CLUSTER_SIZE"]  )  :
                    failure_num = failure_num + 1
                    if ( np.min( np.unique(step.membership, return_counts=True)[1] ) < kwargs["MIN_CLUSTER_SIZE"]  ):
                        if kwargs ["VERBOSE"] >= 2:
                            print ("\t\t{}th step\tTerminated due to number of cluster undersized --MIN_CLUSTER_SIZE in {}-{} :  {} ( < {}) ({})".format(step_index, kwargs["TRIAL"], kwargs["STEP"], np.min( np.unique(step.membership, return_counts=True)[1] ),  kwargs["MIN_CLUSTER_SIZE"], np.unique(step.membership, return_counts=True)[1] ))
                    else:
                        if kwargs ["VERBOSE"] >= 2:
                            print ("\t\t{}th step\tTerminated due to a null clone was found\t{}".format ( step_index, np.unique(step.membership, return_counts=True)  ))
                    i =  step.find_max_likelihood(0, kwargs["STEP"] - 1) 
                    if step.likelihood_record [i] > trial.likelihood_record [ kwargs["TRIAL"]]  : 
                        trial.acc ( step.mixture_record [i],  step.membership_record [i], step.likelihood_record [i], step.membership_p_record [i], step.membership_p_normalize_record [i], 
                                        step.makeone_index_record[i], step.fp_index_record[i],  step_index, step.outlier_index_record[i], step.includeoutlier_record[i], i, kwargs["TRIAL"] )              
                    break
                
                step = Mstep.main(df, np_vaf, step, "Hard", **kwargs)     # M step: Setting new mixture under new membrship

                step.acc(step.mixture, step.membership, step.likelihood, step.membership_p, step.membership_p_normalize, step.makeone_index, step.fp_index, step_index, step.outlier_index, step.includeoutlier, kwargs["STEP"], kwargs["STEP"])       #여기서  step_index,  max_step_index 저장은 별로 안 중요함
                if kwargs ["VERBOSE"] >= 2:
                    print ("\t\t{}th step\tFP_index : {}\tChild_index : {}\tstep.likelihood : {}".format(step_index, step.fp_index, step.makeone_index , round(step.likelihood, 1)))
        

                if miscellaneous.GoStop(step, **kwargs) == "Stop":
                    i =  step.find_max_likelihood(1, kwargs["STEP"]) 
                    if kwargs ["VERBOSE"] >= 2:
                        print ("\t\t\tmax_step : {}".format(i ))
                    trial.acc ( step.mixture_record [i],  step.membership_record [i], step.likelihood_record [i], step.membership_p_record [i], step.membership_p_normalize_record [i], step.makeone_index_record[i], step.fp_index_record[i],  step_index + 1, step.outlier_index_record[i], step.includeoutlier_record[i], i, kwargs["TRIAL"] )
                    trial_index = trial_index + 1
                    failure_num = 0
                    break

            if failure_num >= 2:  # Give up and advance to next trial
                #print ("\t---Give up and advance to next trial")
                
                failure_num = 0
                trial_index = trial_index + 1

            
        i =  trial.find_max_likelihood(0, kwargs["TRIAL_NO"]) 
        if kwargs ["VERBOSE"] >= 1:
            print ("\nWe chose {}th trial, {}th step\n\t(trial.likelihood_record : {})\n\tFP_index : {}\n\tlen(Outlier_index) : {}".format(i, trial.max_step_index_record[i], trial.likelihood_record, trial.fp_index_record[i],  len (trial.outlier_index_record[i] ) ) )

        os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE"]) + "." + str( i ) + "-"  + str(  trial.max_step_index_record [i]  ) + "\(hard\)." + kwargs["IMAGE_FORMAT"] + " " + 
            kwargs["CLEMENT_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE"]) + ".\(hard\)." + kwargs["IMAGE_FORMAT"]  ) 
        
        cluster.acc ( trial.mixture_record [i], trial.membership_record [i], trial.likelihood_record [i], trial.membership_p_record [i], trial.membership_p_normalize_record [i], trial.stepindex_record [i], i, trial.max_step_index_record [i], trial.makeone_index_record[i], trial.fp_index_record[i], trial.includeoutlier_record[i], trial.outlier_index_record [i], **kwargs )  

    return cluster

    #print ("cluster_hard.makeone_index_record : {}".format(cluster_hard.makeone_index_record))
