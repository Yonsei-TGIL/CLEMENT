
def main():
    import os, subprocess, sys, datetime, time, io, contextlib, argparse
    import numpy as np
    import pandas as pd
    if os.getcwd() not in sys.path:
        sys.path.append  ( os.getcwd() )
    import EMhard, EMsoft, Estep, Mstep, Bunch, miscellaneous, datapreparation, phylogeny, visualizationsingle, visualizationsinglesoft, filetype


    pd.options.mode.chained_assignment = None

    kwargs = {}

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument('--INPUT_TSV', type=str, default="/data/project/Alzheimer/YSscript/EM_MRS/CLEMENT/example/2.CellData/MRS_2D/M1-3_M1-8/M1-3_M1-8_input.txt",  help="Path where TSV format file locate. (Important : DIRECTORY (X), File path (O) )")
    parser.add_argument('--CLEMENT_DIR', default="/data/project/Alzheimer/YSscript/EM_MRS/CLEMENT/example/2.CellData/MRS_2D/M1-3_M1-8",   help="Directory where input and output of CLEMENT deposits (Important : DIRECTORY (O), File path (X) )")
    parser.add_argument('--MODE', type=str, choices=["Hard", "Both"], default="Both", help="Selection of clustering method. Default : Both")
    parser.add_argument('--RANDOM_PICK', type=int, default=-1,  help="The number of mutations to alleviate the computational load. Default : -1 (load the whole data)")
    parser.add_argument('--KMEANS_CLUSTERNO',  type=int, default=8,  choices=range(5, 20), help="Number of initial K-means cluster. Recommendation : 5~8 for one-sample, 8-15 for larger-sample. Default : 8")
    parser.add_argument('--NUM_CLONE_TRIAL_START', type=int, default=3,  help="Minimum number of expected cluster_hards (initation of K). Default : 3")
    parser.add_argument('--NUM_CLONE_TRIAL_END', type=int, default=5, choices=range(1, 11), help="Maximum number of expected cluster_hards (termination of K). Default : 5")
    parser.add_argument('--TRIAL_NO', default=5, type=int, choices=range(1, 21),  help="The number of trials in each candidate clone number. DO NOT recommend over 20. Default : 5")
    parser.add_argument('--MAXIMUM_NUM_PARENT',  default=1, type=int,  help="The maximum number of parents in the given samples. Recommendation : 0-2. Default : 1")
    parser.add_argument('--MIN_CLUSTER_SIZE', type=int, default=9, help="The minimum number of membersip in single cluster. Default : 9")
    parser.add_argument('--RANDOM_SEED', type=int, default=1,  help="random_seed for regular random sampling")
    parser.add_argument('--MAKEONE_STRICT', type=int,  choices=[1, 2], default = 1, help="1:strict, 2:lenient. Default : 1")
    parser.add_argument('--TN_CONFIDENTIALITY', default=0.995, type=float, help="Confidentiality that negative being negative (TN). Recommendation : > 0.99. Default : 0.995")
    parser.add_argument('--VERBOSE', type=int, choices=[0, 1, 2, 3], default=2, help="0: Verbose, 3: Concise. Default : 2")


    args = parser.parse_args()
    kwargs["INPUT_TSV"] = args.INPUT_TSV
    kwargs["INPUT_FILETYPE"], kwargs["NUM_BLOCK"] = filetype.main( kwargs["INPUT_TSV"] )
    kwargs["NUM_BLOCK_INPUT"] =  kwargs["NUM_BLOCK"]
    kwargs["MODE"] = args.MODE
    kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] = args.NUM_CLONE_TRIAL_START, args.NUM_CLONE_TRIAL_END
    kwargs["RANDOM_PICK"] = int(args.RANDOM_PICK)
    kwargs["NUM_MUTATION"] = kwargs["RANDOM_PICK"]
    kwargs["MAXIMUM_NUM_PARENT"] = int(args.MAXIMUM_NUM_PARENT)
    kwargs["TN_CONFIDENTIALITY"] = args.TN_CONFIDENTIALITY
    kwargs["TRIAL_NO"] = int(args.TRIAL_NO)
    kwargs["VERBOSE"] = int(args.VERBOSE)
    kwargs["MIN_CLUSTER_SIZE"] = int(args.MIN_CLUSTER_SIZE)
    kwargs["KMEANS_CLUSTERNO"] = args.KMEANS_CLUSTERNO
    kwargs["RANDOM_SEED"] = int(args.RANDOM_SEED)
    kwargs["SCORING"] = False
    kwargs["MAKEONE_STRICT"] = int(args.MAKEONE_STRICT)
    kwargs["CLEMENT_DIR"] = args.CLEMENT_DIR
    if kwargs["CLEMENT_DIR"][-1] == "/":
        kwargs["CLEMENT_DIR"] = kwargs["CLEMENT_DIR"][0:-1]
    kwargs["method"] = "gap+normal"
    kwargs["adjustment"] = "half"
    kwargs["STEP_NO"] = 30
    kwargs["DECISION_STANDARD"]= 0.8 + (0.03 * kwargs ["NUM_BLOCK"])


    print("\n\n\n\nNOW RUNNING IS STARTED  :  {}h:{}m:{}s\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec)))
    print("NUMBER OF INPUT SAMPLES = {}\n\n\n".format( kwargs["NUM_BLOCK"] ))



    print("============================== STEP #1.   DATA EXTRACTION FROM THE ANSWER SET  ==============================")

    for DIR in [ kwargs["CLEMENT_DIR"] + "/trial",  kwargs["CLEMENT_DIR"] + "/Kmeans",   kwargs["CLEMENT_DIR"] + "/candidate",   kwargs["CLEMENT_DIR"] + "/result" ]:
        if os.path.isdir(DIR) == True:
            os.system("rm -rf  " + DIR)
        if os.path.isdir(DIR) == False:
            os.system("mkdir -p " + DIR)

    with open (kwargs["CLEMENT_DIR"] + "/0.commandline.txt", "w" , encoding="utf8" ) as output_logfile:
        print ("python3 CLEMENTDNA.py  --INPUT_TSV {} --NUM_CLONE_TRIAL_START {} --NUM_CLONE_TRIAL_END {}  --RANDOM_PICK {} --MIN_CLUSTER_SIZE {}   --KMEANS_CLUSTERNO {}   --CLEMENT_DIR {}  --MODE {} \
                    --TN_CONFIDENTIALITY {} --MAKEONE_STRICT {} --MAXIMUM_NUM_PARENT {} --TRIAL_NO {}  --VERBOSE {}".format(kwargs["INPUT_TSV"], kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"],  kwargs["RANDOM_PICK"], kwargs["MIN_CLUSTER_SIZE"], kwargs["KMEANS_CLUSTERNO"],   kwargs["CLEMENT_DIR"], kwargs["MODE"], kwargs["TN_CONFIDENTIALITY"], kwargs["MAKEONE_STRICT"], kwargs["MAXIMUM_NUM_PARENT"], kwargs["TRIAL_NO"],   kwargs["VERBOSE"]  ),         file = output_logfile)


    # membership_answer: ['V1', 'V2', 'FP', 'S0', 'V2', 'S0', 'V2', ....
    # membership_answer_numerical : [0 1 2 3 1 3 1 ...

    inputdf, df, np_vaf, np_BQ, membership_answer, mutation_id, samplename_dict_CharacterToNum, kwargs  = datapreparation.main( **kwargs)

    print (kwargs["RANDOM_PICK"])

    membership_answer_numerical = np.zeros( kwargs["RANDOM_PICK"], dtype="int")
    membership_answer_numerical_nofp_index = []



    if type(inputdf) != type(False):
        samplename_dict_NumToCharacter = {v: k for k, v in samplename_dict_CharacterToNum.items()}   # {0: 'FP', 1: 'V2', 2: 'S0', 3: 'V1'}

        print ("samplename_dict_CharacterToNum = {}\nsamplename_dict_NumToCharacter = {}".format( samplename_dict_CharacterToNum, samplename_dict_NumToCharacter ))

        if kwargs["NUM_BLOCK"] == 1:
            x_median = miscellaneous.VAFdensitogram(np_vaf, "INPUT DATA", kwargs["CLEMENT_DIR"] + "/0.inputdata.pdf", **kwargs)
            visualizationsingle.drawfigure_1d(membership_answer_numerical, "ANSWER_SET (n={})".format(kwargs["RANDOM_PICK"]), kwargs["CLEMENT_DIR"] + "/0.inputdata.pdf", np_vaf, samplename_dict_NumToCharacter, False, -1, list (set (membership_answer_numerical)))
        elif kwargs["NUM_BLOCK"] == 2:
            visualizationsingle.drawfigure_2d(membership_answer, "ANSWER_SET (n={})".format(kwargs["RANDOM_PICK"]), kwargs["CLEMENT_DIR"] + "/0.inputdata.pdf", np_vaf, samplename_dict_CharacterToNum, False, -1)
        elif kwargs["NUM_BLOCK"] >= 3:
            visualizationsingle.drawfigure_2d(membership_answer, "ANSWER_SET (n={})".format(kwargs["RANDOM_PICK"]), kwargs["CLEMENT_DIR"] + "/0.inputdata.pdf", np_vaf, samplename_dict_CharacterToNum, False, -1, "SVD")
        subprocess.run(["cp " + kwargs["CLEMENT_DIR"] + "/0.inputdata.pdf  " +  kwargs["CLEMENT_DIR"] + "/candidate/0.inputdata.pdf"], shell=True)
        subprocess.run(["cp " + kwargs["CLEMENT_DIR"] + "/0.inputdata.pdf  " +  kwargs["CLEMENT_DIR"] + "/trial/0.inputdata.pdf"], shell=True)


    START_TIME = datetime.datetime.now()


    np_vaf = miscellaneous.np_vaf_extract(df)
    mixture_kmeans, kwargs["KMEANS_CLUSTERNO"] = miscellaneous.initial_kmeans (np_vaf, kwargs["KMEANS_CLUSTERNO"], kwargs["CLEMENT_DIR"] + "/trial/0.inqitial_kmeans.pdf")  

    cluster_hard = Bunch.Bunch2(**kwargs)
    cluster_soft = Bunch.Bunch2(**kwargs)


    ################### Step 2. Hard clustering ##################
    print ("\n\n\n\n========================================== STEP #2.   EM HARD  ==========================================")

    cluster_hard = EMhard.main (df, np_vaf, np_BQ, mixture_kmeans, **kwargs)



    ################### Step 3. Hard -> Soft clustering ##################
    if kwargs["MODE"] in ["Soft", "Both"]:
        print ("\n\n\n======================================== STEP #3.   EM SOFT  ========================================")

        for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
            print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nNUM_CLONE = {0}".format(NUM_CLONE))
            kwargs["NUM_CLONE"], kwargs["NUM_CLONE_NOMINAL"] = NUM_CLONE, NUM_CLONE
            kwargs["OPTION"] = "soft"

            if cluster_hard.likelihood_record[ NUM_CLONE ] !=  float("-inf"):
                print("\n\n\tSequential Soft clustering (TRIAL_NO = {}, STEP_NO = {})".format ( cluster_hard.trialindex_record[ NUM_CLONE ], cluster_hard.stepindex_record [ NUM_CLONE ] ))
                step_soft = Bunch.Bunch1( kwargs["NUM_MUTATION"] , kwargs["NUM_BLOCK"], NUM_CLONE, cluster_hard.stepindex_record [ NUM_CLONE ] + kwargs["STEP_NO"])
                step_soft.copy (cluster_hard, 0, NUM_CLONE) 

                for step_index in range(1, kwargs["STEP_NO"]): 
                    kwargs["STEP"], kwargs["TRIAL"] = step_index, cluster_hard.trialindex_record[ NUM_CLONE ]
                    kwargs["STEP_TOTAL"] = step_index + cluster_hard.stepindex_record [ NUM_CLONE ] - 1
                    
                    print("\t\tstep #{} ( = TOTAL step #{})".format(kwargs["STEP"], kwargs["STEP_TOTAL"]) )

                    step_soft = Estep.main(df, np_vaf, np_BQ, step_soft, **kwargs)               # E step
                    print ( "\t\t\tAfter the E step: {}\tmakeone_index : {}".format( np.unique(step_soft.membership  , return_counts=True), step_soft.makeone_index ) )
                    step_soft = Mstep.main(df, np_vaf, np_BQ, step_soft, "Soft", **kwargs)     # M step
                    print("\t\t\tAfter the M step : fp_index {}\tmakeone_index {}\tlikelihood : {}".format( step_soft.fp_index, step_soft.makeone_index, round (step_soft.likelihood, 1), step_soft.mixture ))


                    step_soft.acc(step_soft.mixture, step_soft.membership, step_soft.likelihood, step_soft.membership_p, step_soft.membership_p_normalize, step_soft.makeone_index, step_soft.fp_index, step_index + 1, step_soft.fp_member_index, step_soft.includefp, step_soft.fp_involuntary, step_soft.makeone_prenormalization, step_index, step_index)

                    if (miscellaneous.GoStop(step_soft, **kwargs) == "Stop")  :
                        break
                    if ( EMsoft.iszerocolumn (step_soft, **kwargs) == True) :
                        print ("\t\t\t\t→ Terminated due to empty mixture\t{}".format(step_soft.mixture))
                        break
                    if ( len ( set (step_soft.membership) ) < NUM_CLONE ) :
                        print ("\t\t\t\t→ Terminated due to empty clone")
                        break


                step_soft.max_step_index =  step_soft.find_max_likelihood(1, step_soft.stepindex - 2 )  
                i = step_soft.max_step_index

                # If unavaliable to find the possiible answer in soft clustering
                if i == -1:
                    print ("\t\t\tNot available this clone")
                elif  (step_soft.likelihood_record [i]  <= -9999999) :
                    print ("\t\t\tNot available this clone")

                else:  # (Most of the case) Available 
                    os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE"]) + "." + str( kwargs["TRIAL"] ) + "-"  + str(step_soft.max_step_index  + cluster_hard.stepindex_record [ NUM_CLONE ] - 1) + "\(soft\).pdf" + "  " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE"])  + ".\(soft\).pdf"  )
                    cluster_soft.acc ( step_soft.mixture_record [i], step_soft.membership_record [i], step_soft.likelihood_record [i], step_soft.membership_p_record [i], step_soft.membership_p_normalize_record [i], step_soft.stepindex_record[i], cluster_hard.trialindex, step_soft.max_step_index_record[i], step_soft.makeone_index_record[i], step_soft.fp_index_record[i], step_soft.includefp_record[i], step_soft.fp_involuntary_record[i], step_soft.fp_member_index_record[i]   ,**kwargs )


            else:   # If unavaliable to find the possiible answer in hard clustering
                print ("Cant' find possible combination even in Hard clustering.")






    print ("\n\n\n\n==================================== STEP #4.  OPTIMAL K DETERMINATION  =======================================")

    NUM_CLONE_hard , NUM_CLONE_soft = [], []    

    print ("\n\n★★★ Gap Statistics method (Hard clustering)\n")

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        NUM_CLONE_hard, Gap_list = miscellaneous.decision_gapstatistics (cluster_hard, np_vaf, **kwargs)
    print ( f.getvalue() )
    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard.gapstatistics.txt", "w", encoding = "utf8") as gap_myEM:
        print (f.getvalue(), file = gap_myEM)

    if kwargs["MODE"] in ["Soft", "Both"]:
        if kwargs["NUM_BLOCK"] >= 1:
            print ("\n\n\n★★★ XieBeni index method (2D, 3D Soft clustering)\n")

            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                NUM_CLONE_soft = miscellaneous.decision_XieBeni (cluster_soft, np_vaf, **kwargs)

            print ( f.getvalue() )
            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft.xiebeni.txt", "w", encoding = "utf8") as xiebeni_myEM:
                print (f.getvalue(), file = xiebeni_myEM)
            

        if kwargs["NUM_BLOCK"] == 1:
            print ("\n\n\n★★★ Max likelihood method (1D Soft clustering)\n")

            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                NUM_CLONE_soft = miscellaneous.decision_max (cluster_soft, np_vaf, **kwargs)

            print ( f.getvalue() )
            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft.maxlikelihood.txt", "w", encoding = "utf8") as maxlikelihood_myEM:
                print (f.getvalue(), file = maxlikelihood_myEM)
            


    print ("\nCurrent Time : {}h:{}m:{}s    (Time consumed : {})\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))






    if all(x == float("-inf") for x in Gap_list):
        print ("\n\n!!!!!!!!!! Can't determine the clusters by CLEMENT !!!!!!!!!")

    else:
        print ("\n\n\n\n================================ STEP #5-1.  PRINT & VISUALIZING HARD CLUSTERING  ====================================")
        
        DECISION = "hard_1st"


        for i, priority in enumerate(["1st"]):
            if i >= len (NUM_CLONE_hard):
                break
            
            print ( "NUM_CLONE_hard (by order): {}".format(NUM_CLONE_hard))

            if cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] == []:
                if (priority == "1st") & (kwargs["MODE"] in ["Both"]):
                    print ("DECISION\t{}".format(DECISION))
                    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_vs_fuzzy.txt", "w", encoding = "utf8") as output_hard_vs_fuzzy:
                        print ("DECISION\t{}".format(DECISION) , file = output_hard_vs_fuzzy )            
                
            else:
                if (priority == "1st") &  (kwargs["MODE"] in ["Both"]):
                    moved_col_list = miscellaneous.movedcolumn ( cluster_hard, cluster_soft,  NUM_CLONE_hard[i]  )
                    hard_std = np.std(  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
                    soft_std = np.std(  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
                    if ( len (moved_col_list) == 1 )  & (  kwargs["NUM_BLOCK"] == 1 ) :  # Incalculable std in this condition
                        not_moved_col_list = [i for i in list(range ( NUM_CLONE_hard[i] )) if i not in moved_col_list]
                        print ( "moved_col_list = {}\nnot_moved_col_list = {}\ncluster_hard.mixture_record = {}".format (moved_col_list, not_moved_col_list, cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] ))
                        not_moved_col_mean  =  np.mean ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]  [0] [ not_moved_col_list ] ) 
                        moved_col_mean = np.mean ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] [ 0 , moved_col_list ]  ) 
                        print ( "not_moved_col_mean = {}\nmoved_col_mean = {}".format (not_moved_col_mean, moved_col_mean) )

                        if abs (moved_col_mean - not_moved_col_mean) > 0.1:  
                            DECISION = "hard_1st"
                        else:
                            DECISION = "soft_1st"

                        

                    if  (soft_std < hard_std * kwargs["DECISION_STANDARD"]):
                        DECISION = "soft_1st"
                    
                    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_vs_fuzzy.txt", "w", encoding = "utf8") as output_hard_vs_fuzzy:
                        print ( "Moved column : {}".format(moved_col_list), file = output_hard_vs_fuzzy)
                        print ( "Hard (n = {}) : std = {}\tstep = {}\nhard_mixture = {}\n".format( cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( hard_std, 3) ,  cluster_hard.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]   ) , file = output_hard_vs_fuzzy  )
                        print ( "Soft (n = {}) : std = {}\tstep = {}\nhard_mixture = {}".format( cluster_soft.mixture_record [NUM_CLONE_hard[i]].shape[1],  round( soft_std, 3) ,  cluster_soft.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ]   )  , file = output_hard_vs_fuzzy )
                        print ( "ratio : {}".format ( round(soft_std, 3) / round(hard_std, 3) ), file = output_hard_vs_fuzzy )
                        print ("\nsoft 선택 기준 :  < {}\nDECISION\t{}".format(kwargs["DECISION_STANDARD"], DECISION) , file = output_hard_vs_fuzzy )



            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt", "w", encoding = "utf8") as output_myEM:
                print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nmakeone_index\t{}".
                        format(cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1] - int (cluster_hard.includefp_record [ NUM_CLONE_hard[i] ]) , len (cluster_hard.makeone_index_record [NUM_CLONE_hard[i]]),   round((datetime.datetime.now() - START_TIME).total_seconds()), cluster_hard.includefp_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]] , cluster_hard.makeone_index_record [NUM_CLONE_hard[i]]  ), file = output_myEM)
        

            pd.DataFrame(cluster_hard.membership_record [NUM_CLONE_hard[i]]).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_hard.mixture_record [NUM_CLONE_hard[i]],).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_hard.membership_record [NUM_CLONE_hard[i]], return_counts = True ) ).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )

            samplename_dict = {k:k for k in range(0, np.max(cluster_hard.membership_record [NUM_CLONE_hard[i]])+ 1)}


            subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str(NUM_CLONE_hard[i]) +".\(hard\).pdf " + kwargs["CLEMENT_DIR"]+  "/result/CLEMENT_hard_" + priority + ".pdf"], shell = True)
            print ("\n→ Hard {} results printed".format(priority))
                


    if kwargs["MODE"] in ["Soft", "Both"]:       
        print ("\n\n\n\n================================ STEP #5-2.  PRINT & VISUALIZING SOFT CLUSTERING  ====================================")
        for i, priority in enumerate(["1st"]):
            if i >= len (NUM_CLONE_soft):
                break
            
            if cluster_soft.mixture_record [NUM_CLONE_soft[i]] == []:
                print ("Soft : Empty")
                break

            print ( "NUM_CLONE_soft (by order): {}".format(NUM_CLONE_soft))

            with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt", "w", encoding = "utf8") as output_myEM:
                try:
                    print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}".
                        format(cluster_soft.mixture_record [NUM_CLONE_soft[i]].shape[1] - int (cluster_soft.includefp_record [ NUM_CLONE_soft[i] ]) , len (cluster_soft.makeone_index_record [NUM_CLONE_soft[i]]),   round((datetime.datetime.now() - START_TIME).total_seconds()),  cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]]  ), file = output_myEM)
                except:
                    print ("NUM_CLONE\t-1\nNUM_CHILD\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}", file = output_myEM)
        


            pd.DataFrame(cluster_soft.membership_record [NUM_CLONE_soft[i]]).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame(cluster_soft.mixture_record [NUM_CLONE_soft[i]],).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
            pd.DataFrame( np.unique( cluster_soft.membership_record [NUM_CLONE_soft[i]], return_counts = True ) ).to_csv (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )

            samplename_dict = {k:k for k in range(0, np.max(cluster_soft.membership_record [NUM_CLONE_soft[i]])+ 1)}

            if kwargs["NUM_BLOCK"] == 1:
                visualizationsinglesoft.drawfigure_1d (cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]],
                                                                        "CLEMENT_soft : {}".format( round (cluster_soft.likelihood_record [NUM_CLONE_soft[i]] ) ), kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_soft_" + priority + ".pdf", np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]] , cluster_soft.fp_index_record[NUM_CLONE_soft[i]] , cluster_soft.makeone_index_record[NUM_CLONE_soft[i]] )
            elif kwargs["NUM_BLOCK"] == 2:
                visualizationsinglesoft.drawfigure_2d (cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]],
                                                                        "CLEMENT_soft : {}".format( round (cluster_soft.likelihood_record [NUM_CLONE_soft[i]] ) ), kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_soft_" + priority + ".pdf", np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]] , cluster_soft.fp_index_record[NUM_CLONE_soft[i]] , cluster_soft.makeone_index_record[NUM_CLONE_soft[i]] )
            else:
                visualizationsinglesoft.drawfigure_2d (cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.mixture_record [NUM_CLONE_soft[i]], cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]],
                                                                        "CLEMENT_soft : {}".format( round (cluster_soft.likelihood_record [NUM_CLONE_soft[i]]) ), kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_soft_" + priority + ".pdf", np_vaf, samplename_dict, cluster_soft.includefp_record [NUM_CLONE_soft[i]] , cluster_soft.fp_index_record[NUM_CLONE_soft[i]] , cluster_soft.makeone_index_record[NUM_CLONE_soft[i]], "SVD" )

            subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str(NUM_CLONE_soft[i]) +".\(soft\).pdf " + kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_soft_" + priority + ".pdf"], shell = True)
        
            print ("\n→ Soft {} results printed".format(priority))



        
        print ("\n\n\n======================================= STEP #6. DECISION:  HARD VS SOFT  =======================================")
        if soft_std != 0 :
            print ( "DECISION : {}\t\thard_std : {}\tsoft_std : {}\tratio : {}\tcriteria : < {}".format( DECISION, round(hard_std, 3), round(soft_std, 3), round ( round(soft_std, 3) / round(hard_std, 3), 2) , round (kwargs["DECISION_STANDARD"], 2) ))
        else:
            try:
                print ( "DECISION : {}\t\tmoved_col_mean : {}\tnot_moved_col_mean : {}\tcriteria : > 0.1".format( DECISION, moved_col_mean, not_moved_col_mean ) )
            except:
                print ( "DECISION : {}\t\tmoved_col_list : {}".format( DECISION, moved_col_list) )
        
        if "hard" in DECISION:
            NUM_CLONE_CLEMENT = cluster_hard.mixture_record [NUM_CLONE_hard[0]].shape[1] - int (cluster_hard.includefp_record [ NUM_CLONE_hard[0] ])
            NUM_CHILD_CLEMENT = len (cluster_hard.makeone_index_record [NUM_CLONE_hard[0]])
        elif "soft" in DECISION:
            NUM_CLONE_CLEMENT = cluster_soft.mixture_record [NUM_CLONE_soft[0]].shape[1] - int (cluster_soft.includefp_record [ NUM_CLONE_soft[0] ])
            NUM_CHILD_CLEMENT = len (cluster_soft.makeone_index_record [NUM_CLONE_soft[0]])
            

        subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".pdf  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_decision.pdf" ], shell = True)
        subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".results.txt  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_decision.results.txt" ], shell = True)
        subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".membership.txt  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_decision.membership.txt" ], shell = True)
        subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".membership_count.txt  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_decision.membership_count.txt" ], shell = True)
        subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".mixture.txt  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_decision.mixture.txt" ], shell = True)
        with open (kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_decision.results.txt", "a", encoding = "utf8") as output_myEM:
            print ( "DECISION\t{}\nhard_std\t{}\nsoft_std\t{}".format( DECISION, round(hard_std, 3), round(soft_std, 3) ), file = output_myEM)
                



    print("\nCurrent Time : {}h:{}m:{}s     (Time consumed : {})".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME))


if  __name__ == '__main__':
    main()