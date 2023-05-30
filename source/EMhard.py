import os
import numpy as np
import  Estep, Mstep, Bunch, miscellaneous
import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)
warnings.filterwarnings("ignore")


def whether_trial_acc (max_step, step_index, step, trial, **kwargs):
    if step.likelihood_record [max_step] > trial.likelihood_record [ kwargs["TRIAL"]] : 
        print ("\t\t\t✓ max_step : #{}번째 step\t\tstep.likelihood_record [max_step] = {}".format( max_step , round (step.likelihood_record [max_step] , 2) ))
        trial.acc ( step.mixture_record [max_step],  step.membership_record [max_step], step.likelihood_record [max_step], step.membership_p_record [max_step], step.membership_p_normalize_record [max_step], 
                        step.makeone_index_record[max_step], step.fp_index_record[max_step],  step_index, step.fp_member_index_record[max_step], step.includefp_record[max_step], step.fp_involuntary_record[max_step], step.makeone_prenormalization_record[max_step], max_step, kwargs["TRIAL"] )
    
    return step, trial

def checkall (step, **kwargs):
    sum_mixture = np.zeros ( kwargs["NUM_BLOCK"], dtype = "float")
    for i in range (kwargs["NUM_BLOCK"]):
        for j in range (kwargs ["NUM_CLONE"]):
            if j in step.makeone_index:
                sum_mixture[i] += step.mixture[i][j]

    #print ("\t\tstep.makone_index = {}\tsum_mixture = {}".format (step.makeone_index, sum_mixture))
    if kwargs["MAKEONE_STRICT"] == 1:    # BioData, CellData
        makeone_standard = np.array ( [ [0.93, 1.07], [0.84, 1.15] ],dtype = float)  # 앞은 1D, 뒤는 2D, 3D
    elif kwargs["MAKEONE_STRICT"] == 2:   # Moore
        makeone_standard = np.array ( [ [0.77, 1.25], [0.77, 1.25] ],dtype = float) # 앞은 1D, 뒤는 2D, 3D

    if len(step.makeone_index) == 1:  # 단 1개의 clone만을 makeone이라고 designate했을 때는 좀 봐주는게 어떨까?
        makeone_standard = np.array ( [ [0.7, 1.3], [0.7, 1.3] ],dtype = float)

    if kwargs["NUM_BLOCK"] == 1:      # 1D
        if (sum_mixture[0] < makeone_standard[0][0]) | (sum_mixture[0] > makeone_standard[0][1]):
            return False, sum_mixture
        else:
            return True, sum_mixture
    else:
        for i in range( kwargs["NUM_BLOCK"] ):  # 한 block이라도 1에 걸맞지 않은 게 있으면 False를 return함
            if (sum_mixture[i] < makeone_standard[1][0]) | (sum_mixture[i] > makeone_standard[1][1]): 
                return False, sum_mixture
        return True, sum_mixture
    


def main ( df, np_vaf, np_BQ, mixture_kmeans, **kwargs):
    NUM_BLOCK, kwargs["NUM_BLOCK"]= len(df[0]), len(df[0])
    NUM_MUTATION =  kwargs["RANDOM_PICK"]
    kwargs["STEP_NO"] = 30

    cluster = Bunch.Bunch2(**kwargs)

    for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        kwargs["NUM_CLONE"], kwargs["NUM_CLONE_NOMINAL"] = NUM_CLONE, NUM_CLONE
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nNUM_CLONE = {0}".format(NUM_CLONE))
        trial = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["TRIAL_NO"])
        
        if kwargs["KMEANS_CLUSTERNO"] < kwargs["NUM_CLONE"]:    # Kmeans 로 뽑은 개수가 NUM_CLONE보다 작은 경우 뭘 할수가 없다
            continue

        else:  # 정상적인 경우
            trial_index, failure_num = 0, 0
            while trial_index < kwargs["TRIAL_NO"]:
                kwargs["TRIAL"] = trial_index
                print("\t#{0}번째 trial".format(trial_index))

                step = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, kwargs["STEP_NO"])

                step.mixture = miscellaneous.set_initial_parameter(np_vaf, mixture_kmeans, 
                                                                kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-0.initial_kmeans(hard).pdf"  , **kwargs)
                
                if np.any(step.mixture < 0) == True: # 하나라도 음수가 나오는 경우, 그걸 0으로 바꿔주고 일단 시작해보자
                    print ("\t\tset_initial_parameter : 음수가 나옴", end = "\t")
                    print(",".join(str(row) for row in step.mixture ))

                    step.mixture[:, -1][step.mixture[:, -1] < 0] = 0
                    print ("\t\tset_initial_parameter : 그 음수를 0으로 바꾸고 돌림", end = "\t")
                    print(",".join(str(row) for row in step.mixture ))



                for step_index in range(0, kwargs["STEP_NO"]):
                    kwargs["STEP"], kwargs["STEP_TOTAL"] = step_index, step_index
                    kwargs["OPTION"] = "hard"
                    print ("\t\t#{}번째 step".format(step_index))

                    step = Estep.main(df, np_vaf, np_BQ, step, **kwargs)                   # 주어진 mixture 내에서 새 membership 정하기
                    
                    
                    ################################ 탈락조건 2개  (PARENT 개수, MIN_CLUSTER_SIZE) ################################################
                    if (  NUM_CLONE -  len (step.makeone_index) - int(step.includefp)  > kwargs["MAXIMUM_NUM_PARENT"] ):  # NUM_PARENT가 허용된 것보다 더 많을 시
                        failure_num = failure_num + 1
                        print ("\t\t\t♣ STOP:  {}번째 step 도중 E step에서 →  NUM_CHILD = {}\tNUM_PARENT = {}\tFP 여부 = {}\t이라서".format( step_index, len (step.makeone_index), NUM_CLONE -  len (step.makeone_index) - int(step.includefp), step.includefp ))    
                        if kwargs["STEP"] >= 1:
                            max_step =  step.find_max_likelihood_fp_voluntary(0, kwargs["STEP"] - 1)       # 하나 전까지의 max_likelihood를 찾아서 trial에 넣어줄까 말까 고민
                            step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs)   # 그래도 지금까지의 노력을 인정해서 trial에 등록할 수 있는지 보자
                        break
                
                    if ( len ( set (step.membership) ) < NUM_CLONE ) |  ( np.min( np.unique(step.membership, return_counts=True)[1] ) < kwargs["MIN_CLUSTER_SIZE"]  )  :          #  Cluster SIZE가 작을 때
                        failure_num = failure_num + 1
                        extincted_clone_index = np.argmin ( np.unique(step.membership, return_counts=True)[1]  )
                        if ( np.min( np.unique(step.membership, return_counts=True)[1] ) < kwargs["MIN_CLUSTER_SIZE"]  ):
                            if extincted_clone_index == step.fp_index:    
                                print ("\t\t\t♣ STOP:  {}번째 step 도중 E step에서 →  clone {} ( = FP clone) 의 원소가 {}개 ( < {})라서. ({})".format(step_index, extincted_clone_index, np.min( np.unique(step.membership, return_counts=True)[1] ),  kwargs["MIN_CLUSTER_SIZE"], np.unique(step.membership, return_counts=True)[1] ))
                                # fp_involuntary였는데 다음 E step에서 배정도 못 받는다는 것은, 그 전에 강제로 FP로 designate 한게 아예 말이 안된다는 소리.  Best result로 뽑으면 안된다. 즉 fp_voluntary에서만 뽑자
                                max_step =  step.find_max_likelihood_fp_voluntary(0, kwargs["STEP"] - 1)     
                            else:
                                print ("\t\t\t♣ STOP: {}번째 step 도중 E step에서 →  clone {} 의 원소가 {}개 ( < {})라서. ({})".format(step_index, extincted_clone_index, np.min( np.unique(step.membership, return_counts=True)[1] ),  kwargs["MIN_CLUSTER_SIZE"], np.unique(step.membership, return_counts=True)[1] ))
                                max_step =  step.find_max_likelihood_fp_voluntary(0, kwargs["STEP"] - 1)       # 하나 전까지의 max_likelihood를 찾아서 trial에 넣어줄까 말까 고민
                        else:
                            print ("\t\t\t♣ STOP:  {}번째 step 도중 E step에서 →  빈 clone이 생겨서.\t{}".format ( step_index, np.unique(step.membership, return_counts=True)  ))
                            max_step =  step.find_max_likelihood_fp_voluntary(0, kwargs["STEP"] - 1)       # 하나 전까지의 max_likelihood를 찾아서 trial에 넣어줄까 말까 고민
                            
                        step, trial = whether_trial_acc (max_step, step_index, step, trial, **kwargs)     # 그래도 지금까지의 노력을 인정해서 trial에 등록할 수 있는지 보자
                        break
                    ###########################################################################################
                    
                    
                    step = Mstep.main(df, np_vaf, np_BQ, step, "Hard", **kwargs)   # 새 memberhsip에서 새 mixture구하기

                    if kwargs["STEP"] >= 5: #  step 5이상이라고 해서 각 sample당 합친게 너무 얼토당토없으면 기각시킨다
                        if checkall (step, **kwargs) == False:
                            if kwargs["VERBOSE"] >= 2:
                                print ("\t\t➨ STEP {} : mixture가 checkall 기준에 안 맞아서 사이에 있지 않아서 stop\t{}".format( kwargs["STEP"], step.mixture))
                            step.likelihood = -9999999

                    if step.likelihood > -9999990: #  대부분의 경우는 그냥 acc해주면 됨
                        step.acc(step.mixture, step.membership, step.likelihood, step.membership_p, step.membership_p_normalize, step.makeone_index, step.fp_index, step_index, step.fp_member_index, step.includefp, step.fp_involuntary, step.makeone_prenormalization, kwargs["STEP"], kwargs["STEP"])       #여기서  step_index,  max_step_index 저장은 별로 안 중요함
                        if step.fp_involuntary == True:
                            print ("\t\t\t▶ fp_index : {}, makeone_index : {}, parent_index : {}\tlikelihood : {}\tfp_involuntary : {}".format(step.fp_index, step.makeone_index , sorted( list ( set( list (range(0, NUM_CLONE )) ) - set( step.makeone_index ) - set ( [step.fp_index] ) )) ,  round(step.likelihood, 1), step.fp_involuntary ) )
                        else:  # 대 부분
                            print ("\t\t\t▶ fp_index : {}, makeone_index : {}, parent_index : {}\tlikelihood : {}".format(step.fp_index, step.makeone_index , sorted ( list ( set( list (range(0, NUM_CLONE )) ) - set( step.makeone_index ) - set ( [step.fp_index] ) )),  round(step.likelihood, 1)  ) )

                    

                    if miscellaneous.GoStop(step, **kwargs) == "Stop":
                        max_step =  step.find_max_likelihood_fp_voluntary(0, kwargs["STEP"])     # 0 번째 step도 포함하는게 좋겠다. 다만 
                        print ("\t\t\t✓ max_step : #{}번째 step\t\tstep.likelihood_record [max_step] = {}".format( max_step , round (step.likelihood_record [max_step] , 2) ))
                        
                        trial.acc ( step.mixture_record [max_step],  step.membership_record [max_step], step.likelihood_record [max_step], step.membership_p_record [max_step], step.membership_p_normalize_record [max_step], step.makeone_index_record[max_step], step.fp_index_record[max_step],  step_index + 1, step.fp_member_index_record[max_step], step.includefp_record[max_step], step.fp_involuntary_record[max_step], step.makeone_prenormalization_record[max_step], max_step, kwargs["TRIAL"] )
                        trial_index = trial_index + 1
                        failure_num = 0
                        break

                if failure_num >= 1:  # 그냥 포기하고 넘어간다  (1? 2?)
                    if kwargs["VERBOSE"] >= 3:
                        print ("\t\t\t failure_num = 1  → 그냥 포기하고 다음 trial로 넘어간다")
                    
                    # 어차피 default로 float("-inf")가 차 있으니까 딱히 뭘 더 할 필요는 없다
                    failure_num = 0
                    trial_index = trial_index + 1
                
            i =  trial.find_max_likelihood(0, kwargs["TRIAL_NO"]) 
            print ("\n\n\t{}번째 trial, {}번째 step을 선택함\n\t(trial.likelihood_record : {})\n\tFP_index : {}\n\tlen(fp_member_index) : {}".format(i, trial.max_step_index_record[i], np.round ( trial.likelihood_record ), trial.fp_index_record[i],  len (trial.fp_member_index_record[i] ) ) )

            if trial.max_step_index_record [i]  != -1:   # 아예 이 trial에서 못 찾은 경우 blank
                os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE"]) + "." + str( i ) + "-"  + str(  trial.max_step_index_record [i]  ) + "\(hard\).pdf" + " " + 
                    kwargs["CLEMENT_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE"]) + ".\(hard\).pdf"  ) 
            
            cluster.acc ( trial.mixture_record [i], trial.membership_record [i], trial.likelihood_record [i], trial.membership_p_record [i], trial.membership_p_normalize_record [i], trial.stepindex_record [i], i, trial.max_step_index_record [i], trial.makeone_index_record[i], trial.fp_index_record[i], trial.includefp_record[i], trial.fp_involuntary_record[i], trial.fp_member_index_record [i], **kwargs )  

    return cluster

    #print ("cluster_hard.makeone_index_record : {}".format(cluster_hard.makeone_index_record))
