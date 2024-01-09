def recursive (input_containpos, df, np_vaf, np_BQ, subdim, nonzero_dim, index_interest, index_interest_nonzero, kwargs):
    import os, re, subprocess, sys, datetime, time, copy
    import numpy as np
    import pandas as pd   
    import EMhard, Estep, Mstep, Bunch, miscellaneous

    pd.options.mode.chained_assignment = None

    if os.path.exists( kwargs["CLEMENTREC_DIR"] ) == False:
        os.system("rm -rf " + kwargs["CLEMENTREC_DIR"] )
        os.system("mkdir -p " + kwargs["CLEMENTREC_DIR"] )

    print( "\n--- recursive CLEMENT started" )
    print( "\t\tindex_interest = {}\tindex_interest_nonzero = {}".format (len(index_interest), len(index_interest_nonzero) ) )

    if kwargs["VERBOSE"] >= 1:
        print("\n\n--- step #1.  data extracation from the answer set")

    input_containpos_new = input_containpos.loc [ index_interest_nonzero ]
    df_new = [[None] * kwargs["NUM_BLOCK"] for i in range( len(index_interest_nonzero) )]
    for row in range ( len(index_interest_nonzero )):
        for col in range ( kwargs["NUM_BLOCK"] ):
            df_new[row][col] = df[ index_interest_nonzero[ row ] ][col]
    np_vaf_new = np_vaf [index_interest_nonzero, :]
    np_BQ_new =  np_BQ [index_interest_nonzero, :]


    NUM_MUTATION = kwargs["NUM_MUTATION"] = kwargs["RANDOM_PICK"] = len(index_interest_nonzero)
    NUM_BLOCK = kwargs["NUM_BLOCK"]


    if kwargs["VERBOSE"] >= 1:
        print("\tNUM_BLOCK = {}".format(NUM_BLOCK))
        print("\tRANDOM_PICK = {}".format(kwargs["RANDOM_PICK"]))
        print ("\tMIN_CLUSTER_SIZE = {}".format (kwargs["MIN_CLUSTER_SIZE"]))
        print ("\tnp_vaf_new = {}".format(np_vaf_new.shape))
    

    START_TIME = datetime.datetime.now()

    if kwargs["VERBOSE"] >= 1:
        print ("\n\n --- step #2. initial_kmeans ")
    mixture_kmeans, kwargs = miscellaneous.initial_kmeans (input_containpos_new, df_new, np_vaf_new, np_BQ_new, kwargs["CLEMENT_DIR"] + "/trial/0.inqitial_kmeans." + kwargs["IMAGE_FORMAT"], **kwargs)

    cluster_hard = Bunch.Bunch2(**kwargs)
    cluster_soft = Bunch.Bunch2(**kwargs)



    if kwargs["VERBOSE"] >= 1:
        print ("\n --- step #3.   em hard ")

    kwargs["DEBUG"] = False
    cluster_hard = EMhard.main (input_containpos_new, df_new, np_vaf_new, np_BQ_new, mixture_kmeans, **kwargs)
    # subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/candidate  " + kwargs["COMBINED_OUTPUT_DIR"]  ], shell = True)
    # subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/trial  " + kwargs["COMBINED_OUTPUT_DIR"]  ], shell = True)



    if kwargs["VERBOSE"] >= 1:
        print ("\n --- step $4. em soft ")

    for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if kwargs["VERBOSE"] >= 1:
            print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nNUM_CLONE = {0}".format(NUM_CLONE))
        kwargs["NUM_CLONE_ITER"] = NUM_CLONE
        kwargs["NUM_CLONE"] = NUM_CLONE + 1
        kwargs["OPTION"] = "soft"

        if cluster_hard.likelihood_record[ NUM_CLONE ] !=  float("-inf"):
            if kwargs["VERBOSE"] >= 1:
                print("\n\n\tSequential Soft clustering (TRIAL_NO = {}, HARD_STEP = {})".format ( cluster_hard.trialindex_record[ NUM_CLONE ], cluster_hard.stepindex_record [ NUM_CLONE ] ))
            step_soft = Bunch.Bunch1(kwargs["NUM_MUTATION"] , NUM_BLOCK, NUM_CLONE, cluster_hard.stepindex_record [ NUM_CLONE ] + kwargs["STEP_NO"])
            step_soft.copy (cluster_hard, 0, NUM_CLONE)  # 0번 step에 cluster_hard를 복사한다


            for step_index in range(1, kwargs["STEP_NO"]):   # 0번은 채웠으니 1번부터 시작
                kwargs["STEP"], kwargs["TRIAL"] = step_index, cluster_hard.trialindex_record[ NUM_CLONE ]
                kwargs["STEP_TOTAL"] = step_index + cluster_hard.stepindex_record [ NUM_CLONE ] - 1
                
                if kwargs["VERBOSE"] >= 1:
                    print("\t\tStep #{} ( = TOTAL Step #{})".format(kwargs["STEP"], kwargs["STEP_TOTAL"]) )

                step_soft = Estep.main(input_containpos_new, df_new, np_vaf_new, np_BQ_new, step_soft, **kwargs)                   # 주어진 mixture 내에서 새 membership 정하기
                #print ( "\t\t\tEstep.py : {}\tmakeone_index : {}".format( np.unique(step_soft.membership  , return_counts=True), step_soft.makeone_index ) )
                step_soft = Mstep.main(input_containpos_new, df_new, np_vaf_new, np_BQ_new, step_soft, "Soft", **kwargs)     # 새 memberhsip에서 새 mixture구하기
                if kwargs["VERBOSE"] >= 1:
                    print("\t\t\tMstep.py : makeone_index : {}\tfp_index : {}\ttn_index : {}".format( step_soft.makeone_index, step_soft.fp_index, step_soft.tn_index  ))

                if step_soft.makeone_index == []:   # if failed  (첫 3개는 웬만하면 봐주려고 하지만, 그래도 잘 안될 때)
                    if kwargs["VERBOSE"] >= 1:
                        print ("\t\t\t\t→ checkall == False라서 종료\t{}".format(step_soft.mixture))
                    break

            
                step_soft.acc(step_soft.mixture, step_soft.membership, step_soft.posterior, step_soft.likelihood, step_soft.membership_p, step_soft.membership_p_normalize, step_soft.makeone_index, step_soft.tn_index,  step_soft.fp_index, step_index + 1, step_soft.fp_member_index, step_soft.lowvafcluster_index, step_soft.includefp, step_soft.fp_involuntary, step_soft.checkall_strict, step_soft.checkall_lenient, step_index, step_index)

                if (miscellaneous.GoStop(step_soft, **kwargs) == "Stop")  :
                    break
                if ( miscellaneous.iszerocolumn (step_soft, **kwargs) == True) :
                    if kwargs["VERBOSE"] >= 1:
                        print ("\t\t\t\t→ 빈 mixture가 있어서 종료\t{}".format(step_soft.mixture))
                    break
                if ( len ( set (step_soft.membership) ) < kwargs["NUM_CLONE_ITER"] ) :
                    if kwargs["VERBOSE"] >= 1:
                        print ("\t\t\t\t→ 빈 clone이 있어서 종료")
                    break


            i = step_soft.max_step_index = step_soft.find_max_likelihood_strictfirst(1, step_soft.stepindex - 2 )   # 합쳐서 무조건 1이 되게 한다면 현실과 안 맞을수도 있음...
            i = step_soft.max_step_index = step_index - 1

            # soft clustering에서 아예 답을 못 찾을 경우
            if i == 0:
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t\t1번째 soft step부터 망해서 이번 clone은 망함")
            elif  (step_soft.likelihood_record [i]  <= -9999999) :
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t\t모든 step에서 망해서 (-9999999) 이번 clone은 망함")

            else:  # 대부분의경우:  Soft clustering에서 답을 찾은 경우
                if kwargs["VERBOSE"] >= 1:
                    print ("\t\t✓ max_step : Step #{} ( = TOTAL Step #{})\t\tstep.likelihood_record [max_step] = {}".format( i , i + cluster_hard.stepindex_record [ kwargs["NUM_CLONE_ITER"] ] - 1 , round (step_soft.likelihood_record [i] )  ))
                os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE_ITER"]) + "." + str( kwargs["TRIAL"] ) + "-"  + str(step_soft.max_step_index  + cluster_hard.stepindex_record [ kwargs["NUM_CLONE_ITER"] ] - 1) + "\(soft\)." + kwargs["IMAGE_FORMAT"] + "  " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE_ITER"])  + ".\(soft\)." + kwargs["IMAGE_FORMAT"]  )
                cluster_soft.acc ( step_soft.mixture_record [i], step_soft.membership_record [i], step_soft.posterior_record [i], step_soft.likelihood_record [i], step_soft.membership_p_record [i], step_soft.membership_p_normalize_record [i], step_soft.stepindex_record[i], cluster_hard.trialindex, step_soft.max_step_index_record[i], step_soft.makeone_index_record[i], step_soft.checkall_strict_record[i], step_soft.checkall_lenient_record[i], step_soft.tn_index_record[i], step_soft.fp_index_record[i], step_soft.includefp_record[i], step_soft.fp_involuntary_record[i], step_soft.fp_member_index_record[i]   ,**kwargs )


        else:   # hard clustering에서 아예 답을 못 찾은 경우
            if kwargs["VERBOSE"] >= 1:
                print ("Hard clustering에서조차 모두 망해서 이번 clone은 더 돌리기 싫다")




    if kwargs["VERBOSE"] >= 1:
        print ("\n\n\n\n==================================== STEP #5.  OPTIMAL K DETERMINATION  =======================================")

    NUM_CLONE_hard , NUM_CLONE_soft = [], []    # Hard clustering에서의 order, Soft clustering에서의 order

    if kwargs["VERBOSE"] >= 1:
        print ("\n\n★★★ Gap Statistics method (Hard clustering)\n")


    NUM_CLONE_hard = miscellaneous.decision_gapstatistics (cluster_hard, np_vaf_new, np_BQ_new, **kwargs)

    if kwargs["MODE"] in ["Soft", "Both"]:
        if NUM_BLOCK >= 1:
            if kwargs["VERBOSE"] >= 1:
                print ("\n\n\n★★★ XieBeni index method (2D, 3D Soft clustering)\n")
            NUM_CLONE_soft = miscellaneous.decision_XieBeni (cluster_soft, np_vaf_new, **kwargs)

        if NUM_BLOCK == 1:
            if kwargs["VERBOSE"] >= 1:
                print ("\n\n\n★★★ Max likelihood method (1D Soft clustering)\n")
            NUM_CLONE_soft = miscellaneous.decision_max (cluster_soft, np_vaf_new, **kwargs)





    if kwargs["VERBOSE"] >= 1:
        print ("\n\n\n\n=============================================== STEP #6.  SCORING :  EM HARD  =======================================")

    #subprocess.run (["cp -r " +  kwargs["CLEMENT_DIR"]+ "/candidate  " + kwargs["COMBINED_OUTPUT_DIR"] ], shell = True)
    DECISION = "hard_1st"
    i, priority = 0, "1st"

    if kwargs["NUM_BLOCK"] == 1:       # 1D에서는 웬만하면 soft를 신뢰하는게 어떨까?

        if cluster_soft.mixture_record [ NUM_CLONE_soft[0] ] != []:
            DECISION = "soft_1st"
            soft_std = float("inf")
        else:
            DECISION = "hard_1st"

        if kwargs["VERBOSE"] >= 1:
            print ( "DECISION : {}\nhard_mixture : {}\nsoft_mixture : {}".format( DECISION, cluster_hard.mixture_record [NUM_CLONE_hard[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]] ))
            with open ( kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_decision.evidence.txt"  , "w", encoding = "utf8") as output_file:
                print ( "DECISION : {}\nhard_mixture : {}\n\tsoft_mixture : {}".format( DECISION, cluster_hard.mixture_record [NUM_CLONE_hard[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]] ), file = output_file)
            print ("\nNUM_CLONE_soft (by order) : {}\n".format(NUM_CLONE_soft))

    else:  # 2D, 3D에서는 moved column을 본다
        if kwargs["VERBOSE"] >= 1:
            print ("\nNUM_CLONE_hard (by order) : {}\n".format(NUM_CLONE_hard))

        if cluster_soft.mixture_record [ NUM_CLONE_soft[0] ] == []:   # Soft가 다 망했을 경우
            soft_std = float("inf")
            DECSION  = "hard_1st"

        else: # Most of the cases
            moved_col_list = miscellaneous.movedcolumn ( cluster_hard, cluster_soft,  NUM_CLONE_hard[i]  )

            if len (moved_col_list) == 0:
                DECSION  = "hard_1st"
                hard_std = soft_std = float("inf")
            else: # Most of the cases
                hard_std = miscellaneous.std_movedcolumn ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] , moved_col_list )
                if ( cluster_soft.mixture_record [NUM_CLONE_hard[i]] == []): 
                    soft_std = float("inf")
                else:
                    soft_std = miscellaneous.std_movedcolumn ( cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] , moved_col_list )
                    if  ( ( soft_std / hard_std )<  kwargs["DECISION_STANDARD"]):
                        DECISION = "soft_1st"
                    else:
                        DECISION = "hard_1st"

                    if kwargs["VERBOSE"] >= 1:
                        print ( "Moved column : {}".format(moved_col_list) )
                        print ( "Hard (n = {}) : std = {}\tstep = {}\nhard_mixture = {}\n".format( cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( hard_std, 3) ,  cluster_hard.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]   )   )
                        print ( "Soft (n = {}) : std = {}\tstep = {}\nhard_mixture = {}".format( cluster_soft.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( soft_std, 3) ,  cluster_soft.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ]   )   )
                        print ( "ratio : {}".format ( round(soft_std, 3) / round(hard_std, 3) ) )
                        print ("\nsoft 선택 기준 :  < {}\nDECISION\t{}".format( round (kwargs["DECISION_STANDARD"], 2) , DECISION)  )


    if DECISION == "hard_1st":
        NUM_CLONE_recursive =  NUM_CLONE_hard[0]
        mixture_recursive = cluster_hard.mixture_record [ NUM_CLONE_hard[0]  ]
        membership_recursive = cluster_hard.membership_record [ NUM_CLONE_hard[0]  ]
    elif DECISION == "soft_1st":
        NUM_CLONE_recursive =  NUM_CLONE_soft[0]
        mixture_recursive = cluster_soft.mixture_record [ NUM_CLONE_soft[0]  ]
        membership_recursive = cluster_soft.membership_record [ NUM_CLONE_soft[0]  ]


    if np.all(mixture_recursive == 0) == True:     # 전체가 (0,0)이면 무시하고 바로 return
        #return 0, mixture_recursive  # 전체가 (0,0)이면 무시하고 바로 return
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters = kwargs["NUM_CLONE_TRIAL_END"]   , init='k-means++', max_iter=100, random_state=0)  # model generation
        kmeans.fit ( np_vaf_new )  
        membership_kmeans = kmeans.labels_     
        mixture_recursive = kmeans.cluster_centers_.T * 2   
        # 맨 뒤에 FP (0도 넣어줘야 한다. )
        mixture_recursive = np.hstack((mixture_recursive, np.zeros((mixture_recursive.shape[0], 1))))
        print ( "\t실패해서 그냥 Simple Kmeans 돌림 → mixture_recursive = {}".format (  ", ".join(str(np.round(row, 2)) for row in mixture_recursive  ) ) )



    ########### 정사영을 내리지 말고 E step을 한번 더 돌아보고, membership이 부족한 centroid는 빼 준다 (mixture 복구는 돌아가서)
    #if kwargs["VERBOSE"] >= 1:
    print ("\nDECISION = {}".format(DECISION))
    print ("\t(before) mixture_recursive = {}\tshape = {}".format ( ", ".join(str(np.round(row, 2)) for row in mixture_recursive ) , mixture_recursive.shape))
    pd.DataFrame(mixture_recursive).to_csv(kwargs["CLEMENTREC_DIR"] + "/mixture_pre.tsv", index=False, header=False,  sep="\t")
    pd.DataFrame(membership_recursive).to_csv(kwargs["CLEMENTREC_DIR"] + "/membership_pre.txt", index=False, header=False,  sep="\t")
    pd.DataFrame( np.unique( membership_recursive, return_counts=True) ).to_csv(kwargs["CLEMENTREC_DIR"] + "/count_pre.tsv", index=False, header=False,  sep="\t")

    kwargs["NUM_MUTATION"] = len(index_interest)
    kwargs["NUM_CLONE"] = mixture_recursive.shape[1] 
    kwargs["NUM_CLONE_ITER"] = mixture_recursive.shape[1]  - 1
    step_final  = Bunch.Bunch1( kwargs["NUM_MUTATION"]  , kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"] - 1, 2)   # 어차피 2번만 돌아볼 거니까
    step_final.mixture  = copy.deepcopy ( mixture_recursive )
    input_containpos = input_containpos.reset_index(drop = True) 
    kwargs["STEP"] = 0

    df_extract = [[None] * kwargs["NUM_BLOCK"] for i in range( len(index_interest) )]
    for row in range ( len(index_interest) ):
        for col in range ( kwargs["NUM_BLOCK"] ):
            df_extract[row][col] = df[ index_interest[row] ][col]

    step_final = Estep.main( input_containpos = input_containpos.iloc [index_interest].reset_index(drop = True),
                                             df = df_extract,
                                            np_vaf = np_vaf[index_interest, : ],
                                            np_BQ = np_BQ[index_interest,: ], step = step_final, **kwargs)                   # 주어진 mixture 내에서 새 membership 정하기
    
    
    print ( "\t\treturn_counts = {}".format( np.unique(step_final.membership, return_counts=True)[1] ))

    unique_values, counts = np.unique( step_final.membership, return_counts=True)
    pd.DataFrame( np.unique( step_final.membership, return_counts=True) ).to_csv(kwargs["CLEMENTREC_DIR"] + "/count_middle.tsv", index=False, header=False,  sep="\t")
    values_with_count_more_than_MIN_CLUSTER_SIZE = unique_values[counts >= kwargs["MIN_CLUSTER_SIZE"]]
    step_final.mixture = step_final.mixture [ :, values_with_count_more_than_MIN_CLUSTER_SIZE ]
    NUM_CLONE_recursive = step_final.mixture.shape[1]

    #if kwargs["VERBOSE"] >= 1:
    #print ( values_with_count_more_than_MIN_CLUSTER_SIZE )
    print ( "\t(after) mixture_recursive = {}".format (  ", ".join(str(np.round(row, 2)) for row in step_final.mixture  ) ) )


    ### Print results
    
    pd.DataFrame(step_final.mixture).to_csv(kwargs["CLEMENTREC_DIR"] + "/mixture_post.tsv", index=False, header=False,  sep="\t")
    pd.DataFrame(step_final.membership).to_csv(kwargs["CLEMENTREC_DIR"] + "/membership_post.txt", index=False, header=False,  sep="\t")
    pd.DataFrame(np_vaf_new).to_csv(kwargs["CLEMENTREC_DIR"] + "/np_vaf.tsv", index=False, header=False,  sep="\t")
    pd.DataFrame(np_BQ_new).to_csv(kwargs["CLEMENTREC_DIR"] + "/np_BQ.tsv", index=False, header=False,  sep="\t")
    input_containpos.to_csv (kwargs["CLEMENTREC_DIR"] + "/input_containpos.tsv" , sep = "\t", index = False, header = False)
    subprocess.run([ "cp -rf " + kwargs["CLEMENT_DIR"] + "/candidate  "  +  kwargs["CLEMENTREC_DIR"] + "/candidate" ], shell=True)
    subprocess.run([ "cp -rf " + kwargs["CLEMENT_DIR"] + "/trial  "  +  kwargs["CLEMENTREC_DIR"] + "/trial" ], shell=True)
    if DECISION == "hard_1st":
        subprocess.run ([ "cp -rf " +  kwargs["CLEMENTREC_DIR"]+ "/candidate/clone{}.\(hard\).{}".format ( str(NUM_CLONE_hard[0]) ,  kwargs["IMAGE_FORMAT"] )  + " " + kwargs["CLEMENTREC_DIR"]  ], shell = True)
    elif DECISION == "soft_1st":
        subprocess.run ([ "cp -rf " +  kwargs["CLEMENTREC_DIR"]+ "/candidate/clone{}.\(soft\).{}".format ( str(NUM_CLONE_soft[0]) ,  kwargs["IMAGE_FORMAT"] )  + " " + kwargs["CLEMENTREC_DIR"]  ], shell = True)

    return NUM_CLONE_recursive, step_final.mixture




#python3 /data/project/Alzheimer/YSscript/cle/CLEMENT_bm.py --INPUT_TSV /data/project/Alzheimer/CLEMENT/01.INPUT_TSV/3.BioData/Moore_2D/2.all_woMosaic/PD28690/adrenal_gland_zona/fasciculata_L1_glomerulosa_L1_input.txt --NUM_CLONE_TRIAL_START 5 --NUM_CLONE_TRIAL_END 7 --RANDOM_PICK -1 --AXIS_RATIO -1 --PARENT_RATIO 0 --NUM_PARENT 0 --FP_RATIO 0 --FP_USEALL False --DEPTH_CUTOFF 10 --MIN_CLUSTER_SIZE 5 --VERBOSE 1 --KMEANS_CLUSTERNO 8 --RANDOM_SEED 0 --NPVAF_DIR /data/project/Alzheimer/CLEMENT/02.npvaf/3.BioData/Moore_2D_AG/fasciculata_L1_glomerulosa_L1 --SIMPLE_KMEANS_DIR /data/project/Alzheimer/YSscript/cle/data/SIMPLE_KMEANS/3.BioData/Moore_2D_AG/fasciculata_L1_glomerulosa_L1 --CLEMENT_DIR /data/project/Alzheimer/YSscript/cle/data/CLEMENT/3.BioData/Moore_2D_AG/fasciculata_L1_glomerulosa_L1 --SCICLONE_DIR /data/project/Alzheimer/YSscript/cle/data/sciclone/3.BioData/Moore_2D_AG/fasciculata_L1_glomerulosa_L1 --PYCLONEVI_DIR /data/project/Alzheimer/YSscript/cle/data/pyclone-vi/3.BioData/Moore_2D_AG/adrenal_gland_zona/PD28690-fasciculata_L1_glomerulosa_L1 --QUANTUMCLONE_DIR /data/project/Alzheimer/YSscript/cle/data/quantumclone/3.BioData/Moore_2D_AG/fasciculata_L1_glomerulosa_L1 --COMBINED_OUTPUT_DIR /data/project/Alzheimer/CLEMENT/03.combinedoutput/3.BioData/Moore_2D_AG/fasciculata_L1_glomerulosa_L1 --MODE Both --SCORING False --MAKEONE_STRICT 3 --MAXIMUM_NUM_PARENT 1 --TRIAL_NO 8