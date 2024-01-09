import os
import numpy as np
import pandas as pd
import  Estep, Mstep, Bunch, miscellaneous, subprocess
import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)
warnings.filterwarnings("ignore")


def whether_trial_acc (max_step, step_index, step, trial, **kwargs):
    if max_step != -1:
        if step.likelihood_record [max_step] > trial.likelihood_record [ kwargs["TRIAL"]] : 
            print ("\t\t\t✓ (EMhard.py_early termination) max_step : #{}th step\t\tcheckall_strict = {}\t\tstep.likelihood_record [max_step] = {}".format( max_step , step.checkall_strict_record[max_step], np.round (step.likelihood_record [max_step] , 2) ))    
            print ("wheter_trial_acc\t\t\t {}".format (step.checkall_strict_record[max_step] ) )
            trial.acc ( step.mixture_record [max_step],  step.membership_record [max_step], step.posterior_record [max_step], step.likelihood_record [max_step], step.membership_p_record [max_step], step.membership_p_normalize_record [max_step], 
                            step.makeone_index_record[max_step], step.tn_index_record[max_step],   step.fp_index_record[max_step],  step_index, step.fp_member_index_record[max_step], step.lowvafcluster_index_record[max_step], step.includefp_record[max_step], step.fp_involuntary_record[max_step], step.checkall_strict_record[max_step], step.checkall_lenient_record[max_step], max_step, kwargs["TRIAL"] )
    
    return step, trial


def main (input_containpos, df, np_vaf, np_BQ, mixture_kmeans, **kwargs):
    NUM_BLOCK, kwargs["NUM_BLOCK"]= len(df[0]), len(df[0])
    NUM_MUTATION =  kwargs["RANDOM_PICK"]

    cluster = Bunch.Bunch2(**kwargs)

    for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        kwargs["NUM_CLONE"] = kwargs["NUM_CLONE_ITER"] =  NUM_CLONE
        kwargs = miscellaneous.meandepth (**kwargs)
            

        if kwargs["VERBOSE"] >= 1:
            print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nNUM_CLONE = {0}".format(NUM_CLONE))
        trial = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["TRIAL_NO"])
        

        if kwargs["KMEANS_CLUSTERNO"] < kwargs["NUM_CLONE"]:  
            continue

        else:  # Most of the cases
            trial_index, failure_num = 0, 0
            while trial_index < kwargs["TRIAL_NO"]:
                kwargs["TRIAL"] = trial_index

                step = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["STEP_NO"])

                step, kwargs = miscellaneous.set_initial_parameter(np_vaf, mixture_kmeans, 
                                                                kwargs["CLEMENT_DIR"] + "/trial/clone" + str(NUM_CLONE) + "." + str(kwargs["TRIAL"]) + "-0.initial_kmeans(hard)." + kwargs["IMAGE_FORMAT"] ,
                                                                step, trial, **kwargs)
                #subprocess.run(["cp -rf " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str(NUM_CLONE) + "." + str(kwargs["TRIAL"]) + "-0.initial_kmeans\(hard\)." + kwargs["IMAGE_FORMAT"] + " " + kwargs["COMBINED_OUTPUT_DIR"] + "/trial/clone" + str(NUM_CLONE) + "." + str(kwargs["TRIAL"]) + "-0.initial_kmeans\(hard\)." + kwargs["IMAGE_FORMAT"]], shell=True)
                

                if kwargs["VERBOSE"] >= 1:
                    if (kwargs["MAXIMUM_NUM_PARENT"] == 0) & (kwargs["TRIAL"] % 2 == 1):
                        print("\tTrial #{} (make strictly 1)\t{}".format(trial_index, ", ".join(str(np.round(row, 3)) for row in step.mixture )   ) )
                    else: 
                        print("\tTrial #{}\t{}".format(trial_index, ", ".join(str(np.round(row, 3)) for row in step.mixture )   ) )
                                                

                if step.initial_sampling == False:
                    break


                for step_index in range(0, kwargs["STEP_NO"]):
                    kwargs["STEP"], kwargs["STEP_TOTAL"] = step_index, step_index
                    kwargs["OPTION"] = "hard"
                    if (step_index == (kwargs["STEP_NO"] - 1)):  # 맨 뒤까지 오면 종료
                        trial_index += 1
                        continue
                    if kwargs["VERBOSE"] >= 1:
                        if step_index < kwargs ["COMPULSORY_NORMALIZATION"]:
                            print ("\t\tStep #{} (NORMALIZATION & lenient)".format(step_index))
                        else:
                            print ("\t\tStep #{}".format(step_index))


                    # if ( kwargs ["NUM_CLONE_ITER"] == 3 )  & (kwargs["TRIAL"] in [0, 1] ) :
                    #     kwargs["DEBUG"]  = True
                    # else:
                    #     kwargs["DEBUG"] = False

                    kwargs["DEBUG"]  = False

                    step = Estep.main(input_containpos, df, np_vaf, np_BQ, step, **kwargs)  
                    cluster.df.loc [ len(cluster.df.index) ] = [ kwargs["NUM_BLOCK"], kwargs["TRIAL"], kwargs["STEP"], kwargs["STEP_TOTAL"], step.posterior_normalized, step.likelihood, "hard", NUM_CLONE ]   # 맨 끝에 하나씩 추가
                    
                    
                    ################################ Early terminating condition  (NUM_PARENT, MIN_CLUSTER_SIZE) ################################################
                    if (  kwargs["NUM_CLONE"] -  len (step.makeone_index) - 1  > kwargs["MAXIMUM_NUM_PARENT"] ):     #  1st early terminating condition
                        failure_num = failure_num + 1
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t\t♣ STOP:  {}th step,  because in E step →  NUM_CHILD = {}\tNUM_PARENT = {}".format( step_index, len (step.makeone_index), kwargs["NUM_CLONE"] -  len (step.makeone_index) - 1  ))    
                        if kwargs["STEP"] >= 1:
                            print ("kwargs[STEP] = {}\t{}".format(kwargs["STEP"], step.likelihood) )
                            max_step =  step.find_max_likelihood_strictfirst(1, kwargs["STEP"] - 1)   # 이상, 이하
                            step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs) 
                        break
                
                    if  ( len( set ( range (0, kwargs["NUM_CLONE_ITER"]) ) - set(step.membership) )  != 0  ):   # 빈 clone이 있을 때에는 바로 종료
                        failure_num = failure_num + 1
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t\t♣ STOP: {}th step, because in E step →  Clone {}  is empty".format(step_index, set ( range (0, kwargs["NUM_CLONE_ITER"]) ) - set(step.membership) ) )
                        max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ]  =  step.find_max_likelihood_strictfirst(1, kwargs["STEP"] - 1)
                        step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs)
                        break

                    if  ( np.min( np.unique(step.membership, return_counts=True)[1][np.arange(len(set(step.membership))) != step.fp_index] ) < kwargs["MIN_CLUSTER_SIZE"]  ) :          #  2nd early terminating condition  (except for FP index (the last index))
                        if step.less_than_min_cluster_size == True:     # If it has previous history 
                            failure_num = failure_num + 1
                            extincted_clone_index = np.argmin( np.unique(step.membership, return_counts=True)[1][np.arange(len(set(step.membership))) != step.fp_index] )
                            extincted_clone_count = np.min( np.unique(step.membership, return_counts=True)[1][np.arange(len(set(step.membership))) != step.fp_index] )

                            if kwargs["VERBOSE"] >= 1:
                                print ("\t\t\t♣ STOP: {}th step, because in E step →  The number of variants in clone {}  is {}개 ( < {}). ({})".format(step_index, extincted_clone_index, extincted_clone_count,  kwargs["MIN_CLUSTER_SIZE"], np.unique(step.membership, return_counts=True)[1] ))
                            max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ]  =  step.find_max_likelihood_strictfirst(1, kwargs["STEP"] - 2)            # 2번 용서해줬으니 그 전까지 봐야 함
                            step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs)
                            break
                        else:    # Just give a pardon in the first
                            step.less_than_min_cluster_size = True
                            
                    if np.all(step.mixture[:, range(step.mixture.shape[1] - 1)] == 0, axis=0).any():      # centroid가 (0) or (0, 0)  or (0, 0, 0) 이 나오는 경우
                        failure_num = failure_num + 1
                        zero_column_index = np.where( np.all( step.mixture[:, range (step.mixture.shape[1] - 1) ] == 0, axis = 0) )[0]
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t\t♣ STOP: {}th step, because before M step →  clone {}  is no other than zero point ".format(step_index, zero_column_index  ))
                        max_step, max_step_bool =  step.find_max_likelihood_strictfirst(1, kwargs["STEP"] - 1)          
                        step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs)
                        break
                    ###########################################################################################
                    
                    
                    step = Mstep.main(input_containpos, df, np_vaf, np_BQ, step, "Hard", **kwargs)   # M step  (Draw figure + Select makeone )
                    step.acc(step.mixture, step.membership, step.posterior, step.likelihood, step.membership_p, step.membership_p_normalize, step.makeone_index, step.tn_index, step.fp_index, step_index, step.fp_member_index, step.lowvafcluster_index, step.includefp, step.fp_involuntary, step.checkall_strict, step.checkall_lenient, kwargs["STEP"], kwargs["STEP"]) 
                    if kwargs["VERBOSE"] >= 1:
                        print ("\t\t\t\t▶ makeone_index : {}\tparent_index : {}\tfp_index : {}".format( step.makeone_index , sorted ( list ( set( list (range(0, kwargs["NUM_CLONE"] )) ) - set( step.makeone_index ) - set ( [step.fp_index] ) )),  step.fp_index, step.checkall_lenient ) )

                    

                    if  ( (kwargs["STEP"] <  kwargs["COMPULSORY_NORMALIZATION"]) & (step.checkall_lenient == False) ) | ( (kwargs["STEP"] >= kwargs["COMPULSORY_NORMALIZATION"]) & (step.checkall_strict == False) ) | (miscellaneous.GoStop(step, **kwargs) == "Stop"):
                        #max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ] =  step.find_max_likelihood_strictfirst ( 1, kwargs["STEP"] - 1 )     # Excluding 0th step,  Whether this trial's best step is prenormalized or not
                        max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ] =  step.find_max_likelihood_strictfirst_fromtheend ( 1, kwargs["STEP"] - 1 )     # 되도록이면 맨 마지막을 신뢰
                        if max_step == -1: #그래도 안 되면
                            max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ] =  step.find_max_likelihood_strictfirst_fromtheend ( 1, kwargs["STEP"]  )     # 되도록이면 맨 마지막을 신뢰
                        if max_step == -1: #그래도 안 되면
                            max_step, trial.checkall_strict_record[ kwargs["TRIAL"] ] =  step.find_max_likelihood_strictfirst_fromtheend ( 0, kwargs["STEP"] - 1 )     # 되도록이면 맨 마지막을 신뢰

                        trial.acc ( step.mixture_record [max_step],  step.membership_record [max_step], step.posterior_record [max_step], step.likelihood_record [max_step], step.membership_p_record [max_step], step.membership_p_normalize_record [max_step], step.makeone_index_record[max_step], step.tn_index_record[max_step],  step.fp_index_record[max_step],  step_index + 1, step.fp_member_index_record[max_step], step.lowvafcluster_index_record[max_step], step.includefp_record[max_step], step.fp_involuntary_record[max_step], step.checkall_strict_record[max_step], step.checkall_lenient_record[max_step], max_step, kwargs["TRIAL"] )
                        if kwargs["VERBOSE"] >= 1:
                            print ("\t\t✓ (EMhard.py) max_step : #{}th step\t\tcheckall_strict = {}\t\tstep.likelihood_record [max_step] = {}".format( max_step , trial.checkall_strict_record[ kwargs["TRIAL"] ]  , round (step.likelihood_record [max_step] , 2 )  ))
                        
                        trial_index += 1
                        failure_num = 0
                        break

                if failure_num >= 1: 
                    if kwargs["VERBOSE"] >= 3:
                        print ("\t\t\t\tfailure_num = 1  → Give up and pass to the next trial")
                    trial_index += 1
                    failure_num = 0
                

            i, cluster.checkall_strict_record [ kwargs["NUM_CLONE_ITER"] ] =  trial.find_max_likelihood_strictfirst ( 0, kwargs["TRIAL_NO"] )             # Best trial을 찾되, strict == True에 우선권을 줌
        
            if i == -1:      # 이 NUM_CLONE에서 죄다 실패할 경우
                if kwargs["VERBOSE"] >= 1:
                    print ("\n\n\t(EMhard.py) In NUM_CLONE = {}, we couldn't chose any trial".format(kwargs["NUM_CLONE_ITER"], i, trial.max_step_index_record [i], np.round ( trial.likelihood_record ), cluster.checkall_strict_record[ kwargs["NUM_CLONE_ITER"] ] ,trial.fp_index_record[i],  len (trial.fp_member_index_record[i] ) ) )
                cluster.acc ( np.zeros ( (NUM_BLOCK, NUM_CLONE + 1), dtype = "float"), 
                                                    np.zeros ( (NUM_MUTATION), dtype = "int"), 
                                                    float("-inf"), 
                                                    float("-inf"), 
                                                    np.zeros ( (NUM_MUTATION, NUM_CLONE + 1) , dtype = "float"),   #membership_p
                                                    np.zeros ( (NUM_MUTATION, NUM_CLONE + 1), dtype = "float"),    #membership_p_normalize
                                                    -1,  #step_index
                                                    -1,  #trial_index
                                                    -1,   # max_step_index
                                                    [],     # makeone_index
                                                    False,  False,  # checkall_strict_record,  checkall_lenient_record
                                                    [],  -1,   #tn_index, fp_index
                                                    False, False, [], **kwargs )  

            else:  # Most of the cases
                if kwargs["VERBOSE"] >= 1:
                    print ("\n\n\t(EMhard.py) In NUM_CLONE = {}, we chose {}th trial, {}th step\tcheckall_strict = {}, checkall_lenient = {}\n".format(kwargs["NUM_CLONE_ITER"], i, trial.max_step_index_record [i], trial.checkall_strict_record [i], trial.checkall_lenient_record[i]  ) )
                if trial.max_step_index_record [i]  != -1:   # If available in this trial
                    os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str( i ) + "-"  + str(  trial.max_step_index_record [i]  ) + "\(hard\)." + kwargs["IMAGE_FORMAT"] + " " + 
                        kwargs["CLEMENT_DIR"] + "/candidate/clone" + str(kwargs["NUM_CLONE_ITER"]) + ".\(hard\)." + kwargs["IMAGE_FORMAT"]  ) 
                cluster.acc ( trial.mixture_record [i], trial.membership_record [i], trial.posterior_record [i], trial.likelihood_record [i], trial.membership_p_record [i], trial.membership_p_normalize_record [i], trial.stepindex_record [i], i, trial.max_step_index_record [i], trial.makeone_index_record[i], trial.checkall_strict_record[i], trial.checkall_lenient_record[i],   trial.tn_index_record[i],  trial.fp_index_record[i], trial.includefp_record[i], trial.fp_involuntary_record[i], trial.fp_member_index_record [i], **kwargs )  
            

    return cluster

    #print ("cluster_hard.makeone_index_record : {}".format(cluster_hard.makeone_index_record))
