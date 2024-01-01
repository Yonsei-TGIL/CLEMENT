def main():
    import os, re, subprocess, sys, datetime, time, io, contextlib, argparse, math
    import numpy as np
    import pandas as pd

    print ( "Package directory : {}\n".format (  os.path.dirname(__file__) ) )
    SCRIPT_DIR = os.path.dirname(__file__)
    if os.path.dirname(__file__) not in sys.path:
        sys.path.append  ( os.path.dirname(__file__) )
        
    import EMhard, Estep, Mstep, Bunch, miscellaneous, datapreparation, phylogeny, visualizationsingle, visualizationpair, visualizationsinglesoft, filetype, result, scoring

    pd.options.mode.chained_assignment = None

    kwargs = {}

    # python3 /data/project/Alzheimer/YSscript/cle/TOOL/CLEMENT/script/CLEMENT_bm.py	--INPUT_TSV "/data/project/Alzheimer/YSscript/cle/TOOL/CLEMENT/example/2.CellData/MRS_2D/M1-8_M2-4/M1-8_M2-4_input.txt"  --CLEMENT_DIR "/data/project/Alzheimer/YSscript/cle/TOOL/CLEMENT/example/2.CellData/MRS_2D/M1-8_M2-4"   --NUM_CLONE_TRIAL_START 3 		--NUM_CLONE_TRIAL_END 6  --RANDOM_PICK 500 --SCORING True --NUM_PARENT  1

    parser = argparse.ArgumentParser(description='The below is usage direction.')
    parser.add_argument('--INPUT_TSV', type=str, default="/data/project/Alzheimer/EM_cluster/EM_input/MRS_2_sample/M1-3_M1-8_input.txt",  help="Input data whether TSV or VCF. The tool automatically detects the number of samples")
    parser.add_argument('--CLEMENT_DIR', default=None,   help="Directory where input and output of CLEMENT deposits")
    parser.add_argument('--MODE', type=str, choices=["Hard", "Soft", "Both"], default="Both")
    parser.add_argument('--RANDOM_PICK', type=int, default=-1,  help="The number of mutations that are randomly selected in each trials")
    parser.add_argument('--KMEANS_CLUSTERNO',  type=int, default=8,  choices=range(5, 20), help="Number of initial K-means cluster_hard")
    parser.add_argument('--NUM_CLONE_TRIAL_START',  type=int, default=5,  help="Minimum number of expected cluster_hards (initation of K)")
    parser.add_argument('--NUM_CLONE_TRIAL_END', type=int, default=6, choices=range(1, 11), help="Maximum number of expected cluster_hards (termination of K)")
    parser.add_argument('--TRIAL_NO', default=5, type=int, choices=range(1, 31),  help="Trial number in each candidate cluster_hard number. DO NOT recommend over 30")
    parser.add_argument('--MAXIMUM_NUM_PARENT',  default=1, type=int,  help="The maximum number of parents in the given samples.")
    parser.add_argument('--MIN_CLUSTER_SIZE', type=int, default = 9)
    parser.add_argument('--RANDOM_SEED', type=int, default=1,  help="random_seed for regular random sampling")
    parser.add_argument('--MAKEONE_STRICT', type=int,  choices=[1, 2, 3], default = 1)
    parser.add_argument('--FP_PRIOR', default=0.01, type=float, help="Prior of false positive (FP) being generated. Recommendation :  <= 0.1. Default : 0.01")
    parser.add_argument('--TN_PRIOR', default=0.99, type=float, help="Confidentiality that negative being negative (TN). Recommendation : > 0.99. Default : 0.99")
    parser.add_argument('--ALT_THRESHOLD', default = 4, type = int, help = "Default : 3")
    parser.add_argument('--COMPULSORY_NORMALIZATION', default = 3, type = int, help = "Default : 3")
    parser.add_argument('--SOFT_LENIENT', default = 3, type = int, help = "Default : 3")
    parser.add_argument('--FONT_FAMILY', type=str, default="arial", help="Font family that displayed in the plots. Default : arial")
    parser.add_argument('--VISUALIZATION', type=str, default="True", help="Produce image per step")
    parser.add_argument('--IMAGE_FORMAT', type=str, default="jpg", choices = ["jpg", "pdf"], help="Image format that displayed in the plots. Default : jpg")
    parser.add_argument('--MODEL', type=str, default="betabinomial", choices = ["binomial", "betabinomial"], help="")
    parser.add_argument('--VERBOSE', type=int, choices=[0, 1, 2, 3, 4], default=2, help="0: Verbose, 3: Concise. Default : 2")

    parser.add_argument('--AXIS_RATIO', default=-1,  type=float,  help="The fraction of the mutations not shared at least one sample")
    parser.add_argument('--PARENT_RATIO', default=0,  type=float,  help="The fraction of parent clone mutations. If this values is designated, do not set NUM_PARENT")
    parser.add_argument('--NUM_PARENT',  default=0, type=int,  help="The fraction of parent clones being inserted by large order. If this values is designated, do not set PARENT_RATIO")
    parser.add_argument('--FP_RATIO', default=0,  type=float, help="The fraction of false positive mutations regardless of all samples")
    parser.add_argument('--FP_USEALL', default="False", choices=["True", "False"], help="True : extract ALL FPs,   False : Do not extact FPs")
    parser.add_argument('--DEPTH_CUTOFF', default=100, type=int, help="The mutation of which depth below this values is abandoned")
    parser.add_argument('--SCORING', type=str, choices=["True", "False"], default="True", help="True : comparing with the answer set,  False : just visualization")

    args = parser.parse_args()


    kwargs["INPUT_TSV"] = args.INPUT_TSV
    INPUT_TSV = kwargs["INPUT_TSV"]
    INPUT_FILETYPE, NUM_BLOCK = filetype.main(INPUT_TSV)
    kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK
    kwargs["SAMPLENAME"] = SAMPLENAME = INPUT_TSV.split("/")[-1].split("_input")[0]      # 'M1-5_M1-8_input' -> 'M1-5_M1-8'
    kwargs["SEX"] = filetype.sexdetermination(INPUT_TSV)
    kwargs["CLEMENT_DIR"] = args.CLEMENT_DIR
    kwargs["MODE"] = args.MODE
    kwargs["MODEL"] = args.MODEL
    kwargs["RANDOM_PICK"] = int(args.RANDOM_PICK)
    kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] = args.NUM_CLONE_TRIAL_START, args.NUM_CLONE_TRIAL_END
    kwargs["TRIAL_NO"] = int(args.TRIAL_NO)
    kwargs["MAXIMUM_NUM_PARENT"] = int(args.MAXIMUM_NUM_PARENT)
    kwargs["KMEANS_CLUSTERNO"] = args.KMEANS_CLUSTERNO
    kwargs["MIN_CLUSTER_SIZE"] = int(args.MIN_CLUSTER_SIZE)
    kwargs["RANDOM_SEED"] = int(args.RANDOM_SEED)
    kwargs["MAKEONE_STRICT"] = int(args.MAKEONE_STRICT)
    kwargs["FP_PRIOR"] = float(args.FP_PRIOR)
    kwargs["TN_PRIOR"] = float(args.TN_PRIOR)
    kwargs["COMPULSORY_NORMALIZATION"] = int(args.COMPULSORY_NORMALIZATION)
    kwargs["SOFT_LENIENT"] = int(args.SOFT_LENIENT)
    kwargs ["ALT_THRESHOLD"] = int(args.ALT_THRESHOLD)
    kwargs["FONT_FAMILY"] = str(args.FONT_FAMILY)
    kwargs["VISUALIZATION"] = True if args.VISUALIZATION in ["True", "T"] else False
    kwargs["IMAGE_FORMAT"] = str(args.IMAGE_FORMAT)
    kwargs["VERBOSE"] = int(args.VERBOSE)
    kwargs["NPVAF_DIR"] = kwargs["CLEMENT_DIR"]


    kwargs["AXIS_RATIO"] = float(args.AXIS_RATIO)
    kwargs["PARENT_RATIO"] = float(args.PARENT_RATIO)
    kwargs["NUM_PARENT"] = int(args.NUM_PARENT)
    kwargs["FP_RATIO"] = float(args.FP_RATIO)
    kwargs["FP_USEALL"] = args.FP_USEALL
    kwargs["DEPTH_CUTOFF"] = int(args.DEPTH_CUTOFF)
    if args.SCORING in ["False", "false"]:
        kwargs["SCORING"] = False
    elif args.SCORING in ["True", "true"]:
        kwargs["SCORING"] = True


    kwargs["DEBUG"] = False
    kwargs["STEP_NO"] = 30
    kwargs["DECISION_STANDARD"] = 0.15 + (0.03 * kwargs ["NUM_BLOCK"] ** 2)
    kwargs["SINGLE_MIGRATION"] = 0.10
    kwargs["MAKEWHAT"] = float(1.0)
    posterior_sim = 0





    print("\n\n\n\nNOW RUNNING IS STARTED  :  {}h:{}m:{}s\n\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec)))
    print ( "SAMPLENAME : {}".format ( kwargs["SAMPLENAME"]))
    print ( "NUMBER OF INPUT SAMPLES = {}".format(NUM_BLOCK))
    print ( "SEX = {}\n\n\n".format(kwargs["SEX"]))



    print("============================== STEP #1.   DATA EXTRACTION FROM THE ANSWER SET  ==============================")

    for DIR in [ kwargs["CLEMENT_DIR"] + "/trial",  kwargs["CLEMENT_DIR"] + "/Kmeans",   kwargs["CLEMENT_DIR"] + "/candidate",   kwargs["CLEMENT_DIR"]  + "/result" ]:
        if os.path.exists(DIR) == True:
            os.system("rm -rf  " + DIR)
        if os.path.exists(DIR) == False:
            os.system("mkdir -p " + DIR)

    output_logfile = open (kwargs["CLEMENT_DIR"] + "/0.commandline.txt", "w" , encoding="utf8" )
    print ("python3 /data/project/Alzheimer/YSscript/cle/CLEMENT_bm.py  --INPUT_TSV {}  --RANDOM_PICK {} --AXIS_RATIO {} --PARENT_RATIO {} --NUM_PARENT {} --FP_RATIO {} --FP_USEALL {}  --FP_PRIOR {} --TN_PRIOR {}  --DEPTH_CUTOFF {} --MIN_CLUSTER_SIZE {}   --KMEANS_CLUSTERNO {} --RANDOM_SEED {}  --CLEMENT_DIR {}  --IMAGE_FORMAT {} --VISUALIZATION {}  --MODE {} --SCORING {} --MAKEONE_STRICT {} --MAXIMUM_NUM_PARENT {} --MODEL {} --COMPULSORY_NORMALIZATION {} --SOFT_LENIENT {} --NUM_CLONE_TRIAL_START {} --NUM_CLONE_TRIAL_END {} --TRIAL_NO {}   --VERBOSE {}".format( kwargs["INPUT_TSV"],   kwargs["RANDOM_PICK"], kwargs["AXIS_RATIO"],  kwargs["PARENT_RATIO"],  kwargs["NUM_PARENT"],  kwargs["FP_RATIO"],  kwargs["FP_USEALL"],  kwargs ["FP_PRIOR"], kwargs ["TN_PRIOR"], kwargs["DEPTH_CUTOFF"],  kwargs["MIN_CLUSTER_SIZE"], kwargs["KMEANS_CLUSTERNO"],  kwargs["RANDOM_SEED"],    kwargs["CLEMENT_DIR"], kwargs["IMAGE_FORMAT"], kwargs["VISUALIZATION"] , kwargs["MODE"], kwargs["SCORING"], kwargs["MAKEONE_STRICT"], kwargs["MAXIMUM_NUM_PARENT"], kwargs["MODEL"], kwargs["COMPULSORY_NORMALIZATION"], kwargs["SOFT_LENIENT"],  kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"], kwargs["TRIAL_NO"],   kwargs["VERBOSE"]  ),         file = output_logfile)
    output_logfile.close()


    input_containpos, inputdf, df, np_vaf, np_BQ, membership_answer, mixture_answer, mutation_id, kwargs = datapreparation.main( **kwargs)
    NUM_MUTATION = kwargs["NUM_MUTATION"] = kwargs["RANDOM_PICK"] = input_containpos.shape[0]
    print("\n\nRANDOM_PICK = {}\n\nMEAN_DEPTH = {}\n\n".format(kwargs["RANDOM_PICK"], kwargs["MEAN_DEPTH"]))
    with open(kwargs["CLEMENT_DIR"] + "/0.input_mean_depth.txt", "w", encoding="utf8") as result_answer:  
        print( kwargs["MEAN_DEPTH"], file = result_answer)
    
    membership_answer_numerical = np.zeros(kwargs["NUM_MUTATION"], dtype="int")
    membership_answer_numerical_nofp_index = []

    # membership_answer: ['V1', 'V2', 'FP', 'S0', 'V2', 'S0', 'V2', ....
    # membership_answer_numerical : [0 1 2 3 1 3 1 ...


    if type(inputdf) != type(False):
        # {0: 'FP', 1: 'V2', 2: 'S0', 3: 'V1'}
        kwargs["samplename_dict_NumToCharacter"] = {v: k for k, v in kwargs["samplename_dict_CharacterToNum"].items()}

        print ("samplename_dict_CharacterToNum = {}\nsamplename_dict_NumToCharacter = {}".format( kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"] ))

        with open(kwargs["CLEMENT_DIR"] + "/0.input_membership_numerical.txt", "w", encoding="utf8") as result_answer:
            for k in range(kwargs["NUM_MUTATION"]):
                membership_answer_numerical[k] = kwargs["samplename_dict_CharacterToNum"][membership_answer[k]]
                if (membership_answer[k] != "FP"):
                    membership_answer_numerical_nofp_index.append(k)
                print(membership_answer_numerical[k], file=result_answer)

        with open(kwargs["CLEMENT_DIR"] + "/0.input_membership_letter.txt", "w", encoding="utf8") as result_answer:
            for k in range(kwargs["NUM_MUTATION"]):
                print(membership_answer[k], file=result_answer)

        pd.DataFrame(mixture_answer).to_csv (kwargs["CLEMENT_DIR"] + "/0.input_mixture.txt" , sep = "\t", index = False, header = False)
        pd.DataFrame(np_vaf).to_csv (kwargs["CLEMENT_DIR"] + "/0.input_npvaf.txt" , sep = "\t", index = False, header = False)
        pd.DataFrame(np_BQ).to_csv (kwargs["CLEMENT_DIR"] + "/0.input_npBQ.txt" , sep = "\t", index = False, header = False)
        input_containpos.to_csv (kwargs["CLEMENT_DIR"] + "/0.input_containpos.txt" , sep = "\t", index = False, header = False)
        
        with open(kwargs["CLEMENT_DIR"] + "/0.input_sumofmixture.txt", "w", encoding="utf8") as result_answer:
            for i in range(NUM_BLOCK):
                sum_mixture = 0
                for j in range(mixture_answer.shape[1]):
                    if "FP" not in kwargs["samplename_dict_CharacterToNum"]:
                        sum_mixture = sum_mixture + mixture_answer[i][j]
                    else:
                        if kwargs["samplename_dict_CharacterToNum"]["FP"] != j:
                            sum_mixture = sum_mixture + mixture_answer[i][j]
            print("Sum of mixture (except FP) in sample {} : {}".format(i, sum_mixture), file=result_answer)

        with open(kwargs["CLEMENT_DIR"] + "/0.input_membership_count.txt", "w", encoding="utf8") as result_answer:  
            print(np.unique(membership_answer_numerical, return_counts=True)[1], file=result_answer)

        if kwargs["NUM_BLOCK"] == 1:
            x_median = miscellaneous.VAFdensitogram(np_vaf, "INPUT DATA", kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], **kwargs)
            visualizationsingle.drawfigure_1d(membership_answer_numerical, "ANSWER_SET (n={})".format(kwargs["NUM_MUTATION"]), kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, kwargs["samplename_dict_NumToCharacter"], False, -1, [], **kwargs )
        elif kwargs["NUM_BLOCK"] == 2:
            visualizationsingle.drawfigure_2d(membership_answer, "ANSWER_SET (n={})".format(kwargs["NUM_MUTATION"]), kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, kwargs["samplename_dict_CharacterToNum"], False,  -1, "None", **kwargs)
        elif kwargs["NUM_BLOCK"] >= 3:
            visualizationsingle.drawfigure_3d(membership_answer, "ANSWER_SET (n={})".format(kwargs["NUM_MUTATION"]), kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, kwargs["samplename_dict_CharacterToNum"], False,  -1,  **kwargs)
            visualizationsingle.drawfigure_3d_SVD(membership_answer, "ANSWER_SET (n={})".format(kwargs["NUM_MUTATION"]), kwargs["CLEMENT_DIR"] + "/0.inputdata_SVD." + kwargs["IMAGE_FORMAT"], np_vaf, kwargs["samplename_dict_CharacterToNum"], False,  -1,  **kwargs)
            
        
        for COPYTOPATH in [ kwargs["CLEMENT_DIR"] + "/trial", kwargs["CLEMENT_DIR"] + "/candidate"  ]:
            subprocess.run(["cp  " + kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"] + " "  +  COPYTOPATH + "/0.inputdata." + kwargs["IMAGE_FORMAT"] ], shell=True)
            if kwargs["NUM_BLOCK"] >= 3:
                subprocess.run(["cp " + kwargs["CLEMENT_DIR"] + "/0.inputdata_SVD." + kwargs["IMAGE_FORMAT"] + " "  +  COPYTOPATH + "/0.inputdata_SVD." + kwargs["IMAGE_FORMAT"]], shell=True)
        
        if kwargs["SCORING"] == True:
            pd.DataFrame( [kwargs["samplename_dict_NumToCharacter"]] ).to_csv (kwargs["CLEMENT_DIR"] + "/0.input_samplename_dict_NumToCharacter.txt", sep = "\t", index = False)
            pd.DataFrame( [kwargs["samplename_dict_CharacterToNum"]] ).to_csv (kwargs["CLEMENT_DIR"] + "/0.input_samplename_dict_CharacterToNum.txt", sep = "\t", index = False)
        






    START_TIME = datetime.datetime.now()
    break_no = 0

    while True:
        print ("\n\n\n\n=============================================== STEP #2. INITIAL_KMEANS  =======================================")

        print ("\n\nKMEANS_CLUSTERNO : {}".format ( kwargs["KMEANS_CLUSTERNO"]))
        np_vaf = miscellaneous.np_vaf_extract(df)
        mixture_kmeans, kwargs = miscellaneous.initial_kmeans (input_containpos, df, np_vaf, np_BQ, kwargs["CLEMENT_DIR"] + "/trial/0.inqitial_kmeans." + kwargs["IMAGE_FORMAT"], **kwargs)
        for COPYTOPATH in [ kwargs["CLEMENT_DIR"] + "/candidate",  kwargs["CLEMENT_DIR"], kwargs["CLEMENT_DIR"] + "/trial", kwargs["CLEMENT_DIR"] + "/candidate" ]:
            subprocess.run(["cp " + kwargs["CLEMENT_DIR"] + "/trial/0.inqitial_kmeans." + kwargs["IMAGE_FORMAT"] + " "  +  COPYTOPATH + "/0.inqitial_kmeans." + kwargs["IMAGE_FORMAT"]], shell=True)

        for COPYTOPATH in [ kwargs["CLEMENT_DIR"] + "/trial", kwargs["CLEMENT_DIR"] + "/candidate",   kwargs["CLEMENT_DIR"] + "/trial", kwargs["CLEMENT_DIR"] + "/candidate"  ]:
            subprocess.run(["cp  " + kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"] + " "  +  COPYTOPATH + "/0.inputdata." + kwargs["IMAGE_FORMAT"] ], shell=True)
        if kwargs["NUM_BLOCK"] >= 3:
            subprocess.run(["cp " + kwargs["CLEMENT_DIR"] + "/0.inputdata_SVD." + kwargs["IMAGE_FORMAT"] + " "  +  COPYTOPATH + "/0.inputdata_SVD." + kwargs["IMAGE_FORMAT"]], shell=True)

        cluster_hard = Bunch.Bunch2(**kwargs)
        cluster_soft = Bunch.Bunch2(**kwargs)


        print ("\n\n\n======================================== STEP #3.   EM HARD  ========================================\n")
        print("checkall(strict) : 있으면 선정 과정에서 혜택을 받음" )
        print("{} 개 만큼은 강제 normalization + checkall(lenient) 기준으로 pass 따짐.\n그 이후는 normalization 없이 진행 + checkall(strict) 기준으로 pass 따짐".format ( kwargs["COMPULSORY_NORMALIZATION"]  ) )

        cluster_hard = EMhard.main (input_containpos, df, np_vaf, np_BQ, mixture_kmeans, **kwargs)
        subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/candidate  " + kwargs["CLEMENT_DIR"]  ], shell = True)
        subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/trial  " + kwargs["CLEMENT_DIR"]  ], shell = True)
        print ( cluster_hard.likelihood_record )


        if np.any( cluster_hard.likelihood_record != float("-inf") ):  #하나라도 답을 찾아야 break
            break
        if break_no == 1:   # 포기 조건
            break

        kwargs["KMEANS_CLUSTERNO"] = args.KMEANS_CLUSTERNO - 1       # 하나씩 빼면서 수행
        break_no += 1

        if break_no == 2:
            break





    ######################################################## Step 4. Hard -> Soft clustering ########################################################
    if kwargs["MODE"] in ["Soft", "Both"]:
        print ("\n\n\n======================================== STEP #4.   EM SOFT  ========================================\n")
        print("\tcheckall(strict) : 있으면 선정 과정에서 혜택을 받음\tcheckall(lenient) : 다음 step으로 넘어갈 수 있는 기준" )
        print("\t{} 개 만큼은 강제 normalization. 그 이후는 없이 진행".format ( kwargs["COMPULSORY_NORMALIZATION"]  ) )


        for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
            print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\nNUM_CLONE = {0}".format(NUM_CLONE))
            kwargs["NUM_CLONE_ITER"] = NUM_CLONE
            kwargs["NUM_CLONE"] = NUM_CLONE + 1
            kwargs["OPTION"] = "soft"

            kwargs = miscellaneous.meandepth (**kwargs)


            if cluster_hard.likelihood_record[ NUM_CLONE ] !=  float("-inf"):
                print("\n\n\tSequential Soft clustering (TRIAL_NO = {}, HARD_STEP = {})".format ( cluster_hard.trialindex_record[ NUM_CLONE ], cluster_hard.max_step_index_record [ NUM_CLONE ]  ))
                step_soft = Bunch.Bunch1(kwargs["NUM_MUTATION"] , NUM_BLOCK, NUM_CLONE, cluster_hard.stepindex_record [ NUM_CLONE ] + kwargs["STEP_NO"])
                step_soft.copy (cluster_hard, 0, NUM_CLONE)  # 0번 step에 cluster_hard를 복사한다
                step_soft.hard_or_soft = "soft"

                
                print("\t\tHard에서 받아온 것 : {}".format( step_soft.mixture ) )
                for step_index in range(1, kwargs["STEP_NO"]):   # 0번은 채웠으니 1번부터 시작
                    kwargs["STEP"], kwargs["TRIAL"] = step_index, cluster_hard.trialindex_record[ NUM_CLONE ]
                    kwargs["STEP_TOTAL"] = step_index + cluster_hard.stepindex_record [ NUM_CLONE ] - 1
                    
                    if kwargs["STEP"] <= (kwargs["SOFT_LENIENT"] ) :     
                        print("\t\tStep #{} ( = TOTAL Step #{})\tSOFT_LENEINT & observation".format(kwargs["STEP"], kwargs["STEP_TOTAL"]) )
                    else:
                        print("\t\tStep #{} ( = TOTAL Step #{})".format(kwargs["STEP"], kwargs["STEP_TOTAL"]) )

                    step_soft = Estep.main(input_containpos, df, np_vaf, np_BQ, step_soft, **kwargs)                   # 주어진 mixture 내에서 새 membership 정하기
                    cluster_hard.df.loc [ len(cluster_hard.df.index) ] = [ kwargs["NUM_BLOCK"], kwargs["TRIAL"], kwargs["STEP"], kwargs["STEP_TOTAL"], step_soft.posterior_normalized, step_soft.likelihood, "soft", NUM_CLONE ]   # 맨 끝에 하나씩 추가
                    cluster_soft.df.loc [ len(cluster_soft.df.index) ] = [ kwargs["NUM_BLOCK"], kwargs["TRIAL"], kwargs["STEP"], kwargs["STEP_TOTAL"], step_soft.posterior_normalized, step_soft.likelihood, "soft", NUM_CLONE ]   # 맨 끝에 하나씩 추가
                    if ( len ( set (step_soft.membership) ) < kwargs["NUM_CLONE_ITER"] ) :
                        print ("\t\t\t\t→ 빈 clone이 있어서 종료")
                        break
                    step_soft = Mstep.main(input_containpos, df, np_vaf, np_BQ, step_soft, "Soft", **kwargs)     # 새 memberhsip에서 새 mixture구하기
                    print("\t\t\tCLEMENT_bm.py : makeone_index : {}\tfp_index : {}\ttn_index : {}".format( step_soft.makeone_index, step_soft.fp_index, step_soft.tn_index  ))

                    if step_soft.makeone_index == []:   # if failed  (첫 k개는 lenient로 보면서 웬만하면 봐주려고 하지만, 그래도 잘 안될 때)
                        print ("\t\t\t\t→ checkall_lenient == False라서 종료\t{}".format( "\t".join(str(np.round(row, 2)) for row in step_soft.mixture ) ))
                        break
                    if ( miscellaneous.iszerocolumn (step_soft, **kwargs) == True) :
                        print ("\t\t\t\t→ 빈 mixture가 있어서 종료\t{}".format(step_soft.mixture))
                        break
                    if ( len ( set (step_soft.membership) ) < kwargs["NUM_CLONE_ITER"] ) :
                        print ("\t\t\t\t→ 빈 clone이 있어서 종료")
                        break

                    step_soft.acc(step_soft.mixture, step_soft.membership, step_soft.posterior, step_soft.likelihood, step_soft.membership_p, step_soft.membership_p_normalize, step_soft.makeone_index, step_soft.tn_index,  step_soft.fp_index, step_index + 1, step_soft.fp_member_index, step_soft.lowvafcluster_index,  step_soft.includefp, step_soft.fp_involuntary, step_soft.checkall_strict, step_soft.checkall_lenient, step_index, step_index)

                    if (miscellaneous.GoStop(step_soft, **kwargs) == "Stop")  :
                        break


                # i, nnn = step_soft.max_step_index, nnn = step_soft.find_max_likelihood_strictfirst(1, step_soft.stepindex - 2 )   # 합쳐서 무조건 1이 되게 한다면 현실과 안 맞을수도 있음...
                # i = step_soft.max_step_index = step_index - 1     # soft는 그냥 맨 마지막을 제일 좋아해주자
                i = step_soft.max_step_index = step_soft.find_max_likelihood_strictfirst_fromtheend ( 1, kwargs["STEP"] - 1 ) [0]    # 되도록이면 맨 마지막을 신뢰


                # soft clustering에서 아예 답을 못 찾을 경우
                if (kwargs["STEP"] == 1):
                    print ("\n\t1번째 soft step부터 망해서 이번 clone은 망함")
                elif  (step_soft.likelihood_record [i]  <= -9999999) :
                    print ("\n\t모든 step에서 망해서 (-9999999) 이번 clone은 망함")
                else:  # 대부분의경우:  Soft clustering에서 답을 찾은 경우
                    if kwargs["VERBOSE"] >= 1:
                        print ("\n\t✓ (CLEMENT_bm.py) In NUM_CLONE = {}, we chose {}th (= TOTAL {}th)\tstrict = {}, lenient = {}\tstep.likelihood_record [max_step] = {}".format( NUM_CLONE, i , i + cluster_hard.stepindex_record [ kwargs["NUM_CLONE_ITER"] ] - 1 , step_soft.checkall_strict_record [i], step_soft.checkall_lenient_record [i], round ( step_soft.likelihood_record [i] )  ))
                    os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE_ITER"]) + "." + str( kwargs["TRIAL"] ) + "-"  + str(step_soft.max_step_index  + cluster_hard.stepindex_record [ kwargs["NUM_CLONE_ITER"] ] - 1) + "\(soft\)." + kwargs["IMAGE_FORMAT"] + "  " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE_ITER"])  + ".\(soft\)." + kwargs["IMAGE_FORMAT"]  )
                    cluster_soft.acc ( step_soft.mixture_record [i], step_soft.membership_record [i], step_soft.posterior_record [i], step_soft.likelihood_record [i], step_soft.membership_p_record [i], step_soft.membership_p_normalize_record [i], step_soft.stepindex_record[i], cluster_hard.trialindex, step_soft.max_step_index_record[i], step_soft.makeone_index_record[i], step_soft.checkall_strict_record[i], step_soft.checkall_lenient_record[i], step_soft.tn_index_record[i], step_soft.fp_index_record[i], step_soft.includefp_record[i], step_soft.fp_involuntary_record[i], step_soft.fp_member_index_record[i]   ,**kwargs )


            else:   # hard clustering에서 아예 답을 못 찾은 경우
                print ("\tHard clustering에서조차 모두 망해서 이번 clone은 더 돌리기 싫다")


    # subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/candidate  " + kwargs["CLEMENT_DIR"]  ], shell = True)
    # subprocess.run (["cp -r " + kwargs["CLEMENT_DIR"] + "/trial  " + kwargs["CLEMENT_DIR"]  ], shell = True)






    print ("\n\n\n\n==================================== STEP #5.  OPTIMAL K DETERMINATION  =======================================")

    NUM_CLONE_hard , NUM_CLONE_soft = [], []    # Hard clustering에서의 order, Soft clustering에서의 order

    print ("\n\n★★★ Gap Statistics method (Hard clustering)\n")

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        NUM_CLONE_hard = miscellaneous.decision_gapstatistics (cluster_hard, np_vaf, np_BQ, **kwargs)
    print ( f.getvalue() )
    with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_hard.gapstatistics.txt", "w", encoding = "utf8") as gap_CLEMENT:
        print (f.getvalue(), file = gap_CLEMENT)

    subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_hard.gapstatistics.txt  " + kwargs["CLEMENT_DIR"]  + "/trial/CLEMENT_hard.gapstatistics.txt"], shell = True)
    subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_hard.gapstatistics.txt  " + kwargs["CLEMENT_DIR"]  + "/candidate/CLEMENT_hard.gapstatistics.txt"], shell = True)
    subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_hard.gapstatistics.txt  " + kwargs["CLEMENT_DIR"]  + "/candidate/CLEMENT_hard.gapstatistics.txt"], shell = True)

    if kwargs["MODE"] in ["Soft", "Both"]:
        if NUM_BLOCK >= 1:
            print ("\n\n\n★★★ XieBeni index method (2D, 3D Soft clustering)\n")

            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                NUM_CLONE_soft = miscellaneous.decision_XieBeni (cluster_soft, np_vaf, **kwargs)

            print ( f.getvalue() )
            with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.xiebeni.txt", "w", encoding = "utf8") as xiebeni_myEM:
                print (f.getvalue(), file = xiebeni_myEM)
    

        if NUM_BLOCK == 1:
            print ("\n\n\n★★★ Max likelihood method (1D Soft clustering)\n")

            f = io.StringIO()
            with contextlib.redirect_stdout(f):
                NUM_CLONE_soft = miscellaneous.decision_max (cluster_soft, np_vaf, **kwargs)

            print ( f.getvalue() )
            with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.maxlikelihood.txt", "w", encoding = "utf8") as maxlikelihood_myEM:
                print (f.getvalue(), file = maxlikelihood_myEM)
            subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.maxlikelihood.txt  " + kwargs["CLEMENT_DIR"]  + "/CLEMENT_soft.maxlikelihood.txt"], shell = True)
            subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.maxlikelihood.txt  " + kwargs["CLEMENT_DIR"]  + "/candidate/CLEMENT_soft.maxlikelihood.txt"], shell = True)



    print ("\n현재 시각 : {}h:{}m:{}s    (걸린 시간 : {})\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))






    print ("\n\n\n\n=============================================== STEP #6.  SCORING :  EM HARD  =======================================")

    #subprocess.run (["cp -r " +  kwargs["CLEMENT_DIR"]+ "/candidate  " + kwargs["CLEMENT_DIR"] ], shell = True)
    DECISION = "hard_1st"
    max_score_CLEMENT = 0
    moved_col_len = 0
    ARI_CLEMENT_hard, ARI_CLEMENT_soft = 0, 0

    print ("NUM_CLONE_hard (by order) : {}".format(NUM_CLONE_hard))

    if kwargs["MODE"] in ["Hard", "Both"]:
        for i, priority in enumerate(["1st"]):   # "2nd"
            if ( i >= len (NUM_CLONE_hard) ) |  ( cluster_hard.mixture_record [NUM_CLONE_hard[i]] == []) | ( cluster_hard.likelihood_record [NUM_CLONE_hard[i]] == float("-inf")) :
                break

            if cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] == []:       # 그에 해당하는 Soft clustering이 없을 때  BioData면 soft clustering[0] 결과를 신뢰해주자, 나머지면 hard clustering으로 밀고 가자
                print ( cluster_soft.mixture_record [NUM_CLONE_hard[i]] )
                if (priority == "1st") & (kwargs["MODE"] in ["Both"]):
                    DECISION = "hard_1st"
                    print ("DECISION\t{}".format(DECISION))

            else:
                if (priority == "1st") & (kwargs["MODE"] in ["Both"]):
                    DECISION, posterior_sim = miscellaneous.posterior_similarity ( cluster_hard, cluster_soft, NUM_CLONE_hard[i] , **kwargs )
                    if (kwargs["NUM_BLOCK"] >= 2) & ("soft" in DECISION) :     # 2차원 이상에서는  비록 soft를 선택하더라도 NUM_CLONE_hard 기준으로 마련해주자
                        DECISION = "soft(hard)_1st"
                    print ("posterior_sim = {}\tDECISION = {}".format (posterior_sim, DECISION))

            NUM_CLONE_CLEMENT = NUM_CLONE_hard[i]
            NUM_CHILD_CLEMENT = len (cluster_hard.makeone_index_record [ NUM_CLONE_hard[i] ])


            #1. Scoring & Print summary
            if kwargs["SCORING"] == False:
                with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt", "w", encoding = "utf8") as output_file:
                    print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\ndecision\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}\nfp_member_index\t{}".
                            format( NUM_CLONE_CLEMENT , NUM_CHILD_CLEMENT, "hard_" + priority,  round((datetime.datetime.now() - START_TIME).total_seconds()), cluster_hard.includefp_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]] , cluster_hard.tn_index_record [NUM_CLONE_hard[i]] , cluster_hard.makeone_index_record [NUM_CLONE_hard[i]] ,  cluster_hard.fp_member_index_record [NUM_CLONE_hard[i]], ), file = output_file)
        
            elif kwargs["SCORING"] == True:     ################### 정답set과 점수 맞춰보고 색깔 맞추기  (Hard clustering)
                # FP가 있다면 mixture 값을 구해주기
                if NUM_CLONE_hard[i]  in cluster_hard.membership_record [NUM_CLONE_hard[i]] :
                    cluster_hard.mixture_record [NUM_CLONE_hard[i]][:, -1] =  np.mean (  np_vaf[ np.where( cluster_hard.membership_record [NUM_CLONE_hard[i]] == NUM_CLONE_hard[i]  )[0], : ]  * 2, axis = 0 )

                max_score, sample_dict_PtoA, sample_dict_AtoP, score_df = scoring.Scoring ( membership_answer, membership_answer_numerical ,  cluster_hard.membership_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]], 
                                                                                                                                set( list (range(0, NUM_CLONE_hard[i] )) ) - set( cluster_hard.makeone_index_record [NUM_CLONE_hard[i]] ) - set ( [ cluster_hard.fp_index_record [NUM_CLONE_hard[i]]  ] ), **kwargs  )
                if (max_score > max_score_CLEMENT) & ("hard" in DECISION) & (priority == "1st"):
                    max_score_CLEMENT = max_score
                # FP의 mixture 값을 0, 0 으로 원상복구해주기
                cluster_hard.mixture_record [NUM_CLONE_hard[i]][:, -1] =  np.zeros  ( kwargs["NUM_BLOCK"], dtype = "float64") 
                
                # FP 계산 못한다는 걸 감안하고 그냥 ARI 구해보자
                ARI_CLEMENT_hard = result.ARI ( np.array ( membership_answer_numerical)  , 
                                                                        np.array ( cluster_hard.membership_record [NUM_CLONE_hard[i]] ) )
                fARI_CLEMENT_hard = result.fARI ( np.array ( [ membership_answer_numerical [j] for j in membership_answer_numerical_nofp_index  ] ) , 
                                                        np.array ( [ cluster_hard.membership_record [NUM_CLONE_hard[i]] [j] for j in membership_answer_numerical_nofp_index] ),
                                                        kwargs["CLEMENT_DIR"], SCRIPT_DIR,
                                                        transform = True  )
                print ("\nARI_CLEMENT_hard = {}\tfARI_CLEMENT_hard = {}".format ( round (ARI_CLEMENT_hard, 2) , round ( fARI_CLEMENT_hard, 2) ))

                # rescue data가 얼마나 되는지 세보기
                rescue_data = []
                for k in range (kwargs["NUM_MUTATION"]):
                    if ( cluster_hard.membership_record [NUM_CLONE_hard[i]][k] == cluster_hard.fp_index_record [NUM_CLONE_hard[i]] ) :
                        continue
                    if ( cluster_hard.membership_record [NUM_CLONE_hard[i]][k] in cluster_hard.tn_index_record [NUM_CLONE_hard[i]]  ): 
                        continue
                    if np.any( np_vaf[k] == 0 ) == True:
                        rescue_data.append (k)

                print ("\n[■ {} BEST RESULTS]\n\nCLEMENT\t{}/{}\nNUM_CLONE\t{}\nNUM_CHILD\t{}\nARI\t{}\nmakeone_index\t{}\nfp_index\t{}\ntn_index\t{}\nfp_member_index\t{}\nrescue_data\t{}\nlen(rescue_data)\t{}".format(priority, max_score,  kwargs["NUM_MUTATION"], NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, round (ARI_CLEMENT_hard, 2) , cluster_hard.makeone_index_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]], cluster_hard.tn_index_record [NUM_CLONE_hard[i]],  cluster_hard.fp_member_index_record [NUM_CLONE_hard[i]], rescue_data, len(rescue_data)  ))
                print ("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format ( max_score, kwargs["NUM_MUTATION"], sample_dict_AtoP))
                print ("{}".format(score_df))
                with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt", "w", encoding = "utf8") as output_file:
                    print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\ndecision\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}\nfp_member_index\t{}\nrescue_data\t{}\nlen(rescue_data)\t{}".
                        format(NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, "hard_" + priority, max_score, kwargs["NUM_MUTATION"], round(ARI_CLEMENT_hard, 2) ,  round((datetime.datetime.now() - START_TIME).total_seconds()),  cluster_hard.includefp_record [NUM_CLONE_hard[i]], cluster_hard.fp_index_record [NUM_CLONE_hard[i]],  cluster_hard.tn_index_record [NUM_CLONE_hard[i]], cluster_hard.makeone_index_record [NUM_CLONE_hard[i]], cluster_hard.fp_member_index_record [NUM_CLONE_hard[i]], rescue_data, len(rescue_data) ), file = output_file)

            


            #2. Save results
            subprocess.run (["cp " + kwargs["CLEMENT_DIR"] + "/result/CLEMENT_hard_" + priority + ".results.txt  " + kwargs["CLEMENT_DIR"]  + "/CLEMENT_hard_" + priority + ".results.txt"], shell = True)
            for COPYTOPATH in [ kwargs["CLEMENT_DIR"], kwargs["CLEMENT_DIR"], kwargs["CLEMENT_DIR"] +"/result" ]:
                pd.DataFrame(cluster_hard.membership_record [NUM_CLONE_hard[i]]).to_csv ( COPYTOPATH + "/CLEMENT_hard_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )            
                pd.DataFrame(cluster_hard.mixture_record [NUM_CLONE_hard[i]]).to_csv ( COPYTOPATH + "/CLEMENT_hard_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )            
                pd.DataFrame( np.unique( cluster_hard.membership_record [NUM_CLONE_hard[i]], return_counts = True ) ).to_csv ( COPYTOPATH + "/CLEMENT_hard_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
                if kwargs["SCORING"] == True:
                    pd.DataFrame(score_df).to_csv ( COPYTOPATH + "/CLEMENT_hard_" + priority + ".scoredf.tsv", index = False, header= True,  sep = "\t")            
                cluster_hard.df.to_csv ( COPYTOPATH + "/CLEMENT_hard_" + priority + ".df.tsv", index = False, header= False,  sep = "\t" )

            
            #3. Visualization
            subprocess.run (["cp -rf " +  kwargs["CLEMENT_DIR"]+ "/candidate/clone" + str( NUM_CLONE_hard[i] ) + ".\(hard\)." + kwargs["IMAGE_FORMAT"] + " "  + kwargs["CLEMENT_DIR"]+ "/CLEMENT_hard_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)    


            #4. Phylogeny reconstruction
            ISPARENT = True if NUM_CLONE_CLEMENT  > NUM_CHILD_CLEMENT else False

            if ISPARENT == True:
                print ("\n\n\n\n\t▶▶▶  PHYLOGENY RECONSTRUCTION IN HARD CLUSTERING ◀◀◀")
                kwargs["PHYLOGENY_DIR"] = kwargs["CLEMENT_DIR"] + "/CLEMENT_hard_" + priority + ".phylogeny.txt"
                f = io.StringIO()
                with contextlib.redirect_stdout(f):
                    membership_child = set ( cluster_hard.makeone_index_record[ NUM_CLONE_hard[i] ] )            # child의 번호만 뽑아준다 ( e.g.  0, 1, 3)
                    membership_outside = set (range (0, NUM_CLONE_hard [i] )) - membership_child - set ( [cluster_hard.fp_index_record[ NUM_CLONE_hard[i] ] ] )   # outlier의 번호만 뽑아준다 (e.g. 2)

                    g = phylogeny.main(membership_child, membership_outside, cluster_hard.mixture_record[ NUM_CLONE_hard[i] ],  **kwargs)

                print ( f.getvalue() )
                with open ( kwargs["PHYLOGENY_DIR"] , "w", encoding = "utf8") as phylogeny_file:
                    print (f.getvalue(), file = phylogeny_file)
                subprocess.run (["cp " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["CLEMENT_DIR"]  + "/CLEMENT_hard_" + priority + ".phylogeny.txt"], shell = True)
                subprocess.run (["cp " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_hard_" + priority + ".phylogeny.txt"], shell = True)
                        
            print ("Hard {} results printed".format(priority))






    # DECISION = "soft(hard)_1st"


    if kwargs["MODE"] in ["Soft", "Both"]:
        print ("\n\n\n\n=============================================== STEP #7.  SCORING :  EM SOFT  =======================================")

        if DECISION == "soft(hard)_1st":
            print ("\n\nNUM_CLONE_soft  : {}".format(NUM_CLONE_hard [0] ))
            if ( cluster_soft.mixture_record [NUM_CLONE_hard[0]] == [] ):
                DECISION = "soft_1st"
            else:   # 정말로 있을 경우에만 수행
                NUM_CLONE_CLEMENT = NUM_CLONE_hard[0]
                NUM_CHILD_CLEMENT = len (cluster_soft.makeone_index_record [ NUM_CLONE_hard[0] ])

                #1. Scoring & Print summary
                if kwargs["SCORING"] == False:
                    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft(hard)_1st.results.txt", "w", encoding = "utf8") as output_file:
                        try:
                            print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\ndecision\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}\nfp_member_index\t{}".
                                format( NUM_CLONE_CLEMENT , NUM_CHILD_CLEMENT,  "soft(hard)_1st", round((datetime.datetime.now() - START_TIME).total_seconds()),  cluster_soft.includefp_record [NUM_CLONE_CLEMENT ], cluster_soft.fp_index_record [NUM_CLONE_CLEMENT ] , cluster_soft.tn_index_record [NUM_CLONE_CLEMENT ] , cluster_soft.makeone_index_record [NUM_CLONE_CLEMENT ],  cluster_soft.fp_member_index_record [NUM_CLONE_CLEMENT ],   ), file = output_file)
                        except:
                            print ("NUM_CLONE\t-1\nNUM_CHILD\t{}\ndecision\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}", file = output_file)

                elif kwargs["SCORING"] == True:  #정답set과 점수 맞춰보고 색깔 맞추기  (Soft clustering)
                    # FP가 있다면 mixture 값을 구해주기
                    if NUM_CLONE_CLEMENT  in cluster_hard.membership_record [NUM_CLONE_CLEMENT] :
                        cluster_soft.mixture_record [NUM_CLONE_CLEMENT][:, -1] =  np.mean (  np_vaf[ np.where( cluster_soft.membership_record [NUM_CLONE_CLEMENT] == NUM_CLONE_CLEMENT  )[0], : ]  * 2, axis = 0 )

                    try:
                        score_df, score = scoring.mixturebased(mixture_answer, cluster_soft.mixture_record [NUM_CLONE_CLEMENT], membership_answer, cluster_soft.membership_record [NUM_CLONE_CLEMENT], kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"], cluster_soft.includefp_record [NUM_CLONE_CLEMENT], cluster_soft.fp_index_record [NUM_CLONE_CLEMENT], "CLEMENT", **kwargs)        
                    except:
                        score_df = pd.DataFrame(  columns = ["answer", "predicted",  "shared"] )
                    max_score, sample_dict_PtoA, sample_dict_AtoP, score_df = scoring.Scoring ( membership_answer, membership_answer_numerical ,  cluster_soft.membership_record [NUM_CLONE_CLEMENT], cluster_soft.fp_index_record [NUM_CLONE_CLEMENT], 
                                                                                                                                    set( list (range(0, NUM_CLONE_CLEMENT )) ) - set( cluster_soft.makeone_index_record [NUM_CLONE_CLEMENT] ) - set ( [ cluster_soft.fp_index_record [NUM_CLONE_CLEMENT]  ] ), **kwargs  )
                    if (max_score_CLEMENT < max_score) & (DECISION == "soft(hard)_1st") :
                        max_score_CLEMENT = max_score
                    # FP의 mixture 값을 0, 0 으로 원상복구해주기
                    cluster_soft.mixture_record [NUM_CLONE_CLEMENT][:, -1] =  np.zeros  ( kwargs["NUM_BLOCK"], dtype = "float64") 


                    ARI_CLEMENT_soft = result.ARI ( np.array ( [ membership_answer_numerical [j] for j in membership_answer_numerical_nofp_index  ] ) , 
                                                                            np.array ( [ cluster_soft.membership_record [NUM_CLONE_CLEMENT] [j] for j in membership_answer_numerical_nofp_index] ) )
                    fARI_CLEMENT_soft = result.fARI ( np.array ( [ membership_answer_numerical [j] for j in membership_answer_numerical_nofp_index  ] ) , 
                                                                        np.array ( [ cluster_soft.membership_p_record [NUM_CLONE_CLEMENT] [j] for j in membership_answer_numerical_nofp_index] ),
                                                                        kwargs["CLEMENT_DIR"], SCRIPT_DIR,
                                                                        transform = False  )

                    #print ("\nARI_CLEMENT_soft = {}\tfARI_CLEMENT_soft = {}".format ( round (ARI_CLEMENT_soft, 2) , round (fARI_CLEMENT_soft, 2) ))

                    print ("\n[■ {} BEST RESULTS]\n\nCLEMENT\t{}/{}\nNUM_CLONE\t{}\nNUM_CHILD\t{}\nARI\t{}\nmakeone_index\t{}\nfp_index\t{}\ntn_index\t{}\nfp_member_index\t{}".format(priority, max_score,  kwargs["NUM_MUTATION"], NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, round ( ARI_CLEMENT_soft, 2) , cluster_soft.makeone_index_record [ NUM_CLONE_CLEMENT ], cluster_soft.fp_index_record [ NUM_CLONE_CLEMENT ], cluster_soft.tn_index_record [ NUM_CLONE_CLEMENT ],  cluster_soft.fp_member_index_record [ NUM_CLONE_CLEMENT ]  ))
                    print ("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format ( max_score, kwargs["NUM_MUTATION"], sample_dict_AtoP))
                    print ("{}".format(score_df))
                    
                    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft(hard)_1st.results.txt", "w", encoding = "utf8") as output_file:
                        print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\ndecision\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}\nfp_member_index\t{}".
                            format(NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, "soft(hard)_1st", max_score, kwargs["NUM_MUTATION"],  ARI_CLEMENT_soft, round((datetime.datetime.now() - START_TIME).total_seconds()), cluster_soft.includefp_record [ NUM_CLONE_CLEMENT ], cluster_soft.fp_index_record [ NUM_CLONE_CLEMENT ], cluster_soft.tn_index_record [ NUM_CLONE_CLEMENT ], cluster_soft.makeone_index_record [ NUM_CLONE_CLEMENT ],  cluster_soft.fp_member_index_record [ NUM_CLONE_CLEMENT ]   ), file = output_file)
                
                #2. Save results
                #subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft(hard)_1st.results.txt  " + kwargs["CLEMENT_DIR"]  + "/CLEMENT_soft(hard)_1st.results.txt"], shell = True)
                for COPYTOPATH in [ kwargs["CLEMENT_DIR"],  kwargs["CLEMENT_DIR"] +"/result" ]:
                    pd.DataFrame(cluster_soft.membership_record [NUM_CLONE_CLEMENT]).to_csv (COPYTOPATH + "/CLEMENT_soft(hard)_1st.membership.txt", index = False, header= False,  sep = "\t" )
                    pd.DataFrame(cluster_soft.mixture_record [NUM_CLONE_CLEMENT],).to_csv ( COPYTOPATH + "/CLEMENT_soft(hard)_1st.mixture.txt", index = False, header= False,  sep = "\t" )
                    pd.DataFrame( np.unique( cluster_soft.membership_record [NUM_CLONE_CLEMENT], return_counts = True ) ).to_csv ( COPYTOPATH + "/CLEMENT_soft(hard)_1st.membership_count.txt", index = False, header= False,  sep = "\t" )
                    pd.DataFrame( cluster_soft.membership_p_normalize_record [NUM_CLONE_CLEMENT] ).to_csv ( COPYTOPATH + "/CLEMENT_soft(hard)_1st.posterior_allsample_normalize.txt", index = False, header= False,  sep = "\t" )            
                    if kwargs["SCORING"] == True:
                        pd.DataFrame(score_df).to_csv ( COPYTOPATH + "/CLEMENT_soft(hard)_1st.scoredf.tsv", index = False, header= True,  sep = "\t")
                    cluster_soft.df.to_csv ( COPYTOPATH + "/CLEMENT_soft(hard)_1st.df.tsv", index = False, header= False,  sep = "\t" )


                # 3.Visualization
                subprocess.run (["cp " +  kwargs["CLEMENT_DIR"] + "/candidate/clone" + str(NUM_CLONE_CLEMENT ) + ".\(soft\)." + kwargs["IMAGE_FORMAT"] + " "  + kwargs["CLEMENT_DIR"] + "/CLEMENT_soft\(hard\)_1st." + kwargs["IMAGE_FORMAT"]], shell = True)    


                #4. Phylogeny reconstruction
                ISPARENT = True if NUM_CLONE_CLEMENT  > NUM_CHILD_CLEMENT else False
                
                if ISPARENT == True:
                    print ("\n\n\n\n\t▶▶▶  PHYLOGENY RECONSTRUCTION IN SOFT CLUSTERING ◀◀◀")
                    kwargs["PHYLOGENY_DIR"] = kwargs["CLEMENT_DIR"] + "/CLEMENT_soft(hard)_1st.phylogeny.txt"
                    f = io.StringIO()
                    with contextlib.redirect_stdout(f):
                        membership_child = set ( cluster_soft.makeone_index_record[ NUM_CLONE_CLEMENT ] )            # child의 번호만 뽑아준다 ( e.g.  0, 1, 3)
                        membership_outside = set (range (0, NUM_CLONE_CLEMENT )) - membership_child - set ( [cluster_soft.fp_index_record[NUM_CLONE_CLEMENT ] ] )   # outlier의 번호만 뽑아준다 (e.g. 2)

                        g = phylogeny.main(membership_child, membership_outside, cluster_soft.mixture_record[NUM_CLONE_CLEMENT ],  **kwargs)

                    print ( f.getvalue() )
                    with open ( kwargs["PHYLOGENY_DIR"] , "w", encoding = "utf8") as phylogeny_file:
                        print (f.getvalue(), file = phylogeny_file)
                    subprocess.run (["cp -rf " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["CLEMENT_DIR"]  + "/CLEMENT_soft\(hard\)_1st.phylogeny.txt"], shell = True)
                    subprocess.run (["cp -rf " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_soft\(hard\)_1st.phylogeny.txt"], shell = True)
                
                print ("\tSoft {} results printed".format(priority))



        if DECISION == "soft_1st":
            print ("\n\nNUM_CLONE_soft (by order) : {}".format(NUM_CLONE_soft))
            for i, priority in enumerate(["1st"]):        
                if ( i >= len (NUM_CLONE_soft)) | ( cluster_hard.mixture_record [NUM_CLONE_soft[i]] == [] ): 
                    break
                if len (cluster_soft.makeone_index_record [NUM_CLONE_soft[i]]) == 0:
                    print ( "NUM_CHILD = 0" )
                    print ( cluster_soft.mixture_record [NUM_CLONE_soft[i]] )
                    break

                NUM_CLONE_CLEMENT = NUM_CLONE_soft[i]
                NUM_CHILD_CLEMENT = len (cluster_soft.makeone_index_record [ NUM_CLONE_soft[i] ])

                #1. Scoring & Print summary
                if kwargs["SCORING"] == False:
                    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt", "w", encoding = "utf8") as output_file:
                        try:
                            print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\ndecision\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}\nfp_member_index\t{}".
                                format( NUM_CLONE_CLEMENT , NUM_CHILD_CLEMENT,  "soft_" + priority, round((datetime.datetime.now() - START_TIME).total_seconds()),  cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]] , cluster_soft.tn_index_record [NUM_CLONE_soft[i]] , cluster_soft.makeone_index_record [NUM_CLONE_soft[i]],  cluster_soft.fp_member_index_record [NUM_CLONE_soft[i]],   ), file = output_file)
                        except:
                            print ("NUM_CLONE\t-1\nNUM_CHILD\t{}\ndecision\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}", file = output_file)


                elif kwargs["SCORING"] == True:  #정답set과 점수 맞춰보고 색깔 맞추기  (Soft clustering)
                    # FP가 있다면 mixture 값을 구해주기
                    if NUM_CLONE_soft[i]  in cluster_hard.membership_record [NUM_CLONE_soft[i]] :
                        cluster_soft.mixture_record [NUM_CLONE_soft[i]][:, -1] =  np.mean (  np_vaf[ np.where( cluster_soft.membership_record [NUM_CLONE_soft[i]] == NUM_CLONE_soft[i]  )[0], : ]  * 2, axis = 0 )

                    try:
                        score_df, score = scoring.mixturebased(mixture_answer, cluster_soft.mixture_record [NUM_CLONE_soft[i]], membership_answer, cluster_soft.membership_record [NUM_CLONE_soft[i]], kwargs["samplename_dict_CharacterToNum"], kwargs["samplename_dict_NumToCharacter"], cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]], "CLEMENT", **kwargs)        
                    except:
                        score_df = pd.DataFrame(  columns = ["answer", "predicted",  "shared"] )
                    max_score, sample_dict_PtoA, sample_dict_AtoP, score_df = scoring.Scoring ( membership_answer, membership_answer_numerical ,  cluster_soft.membership_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]], 
                                                                                                                                    set( list (range(0, NUM_CLONE_soft[i] )) ) - set( cluster_soft.makeone_index_record [NUM_CLONE_soft[i]] ) - set ( [ cluster_soft.fp_index_record [NUM_CLONE_soft[i]]  ] ), **kwargs  )
                    if (max_score_CLEMENT < max_score) & (DECISION == "soft_1st") & (priority == "1st"):
                        max_score_CLEMENT = max_score
                    # FP의 mixture 값을 0, 0 으로 원상복구해주기
                    cluster_soft.mixture_record [NUM_CLONE_soft[i]][:, -1] =  np.zeros  ( kwargs["NUM_BLOCK"], dtype = "float64") 
        

                    ARI_CLEMENT_soft = result.ARI ( np.array ( [ membership_answer_numerical [j] for j in membership_answer_numerical_nofp_index  ] ) , 
                                                                            np.array ( [ cluster_soft.membership_record [NUM_CLONE_soft[i]] [j] for j in membership_answer_numerical_nofp_index] ) )
                    fARI_CLEMENT_soft = result.fARI ( np.array ( [ membership_answer_numerical [j] for j in membership_answer_numerical_nofp_index  ] ) , 
                                                                        np.array ( [ cluster_soft.membership_p_record [NUM_CLONE_soft[i]] [j] for j in membership_answer_numerical_nofp_index] ),
                                                                        kwargs["CLEMENT_DIR"], SCRIPT_DIR,
                                                                        transform = False  )

                    #print ("\nARI_CLEMENT_soft = {}\tfARI_CLEMENT_soft = {}".format ( round (ARI_CLEMENT_soft, 2) , round (fARI_CLEMENT_soft, 2) ))

                    print ("\n[■ {} BEST RESULTS]\n\nCLEMENT\t{}/{}\nNUM_CLONE\t{}\nNUM_CHILD\t{}\nARI\t{}\nmakeone_index\t{}\nfp_index\t{}\ntn_index\t{}\nfp_member_index\t{}".format(priority, max_score,  kwargs["NUM_MUTATION"], NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, round ( ARI_CLEMENT_soft, 2) , cluster_soft.makeone_index_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]], cluster_soft.tn_index_record [NUM_CLONE_soft[i]],  cluster_soft.fp_member_index_record [NUM_CLONE_soft[i]]  ))
                    print ("(모두 다 돌려본 결과) score : {}점 / {}점,  변환표 (A to P)= {}\n".format ( max_score, kwargs["NUM_MUTATION"], sample_dict_AtoP))
                    print ("{}".format(score_df))
                    
                    with open (kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt", "w", encoding = "utf8") as output_file:
                        print ("NUM_CLONE\t{}\nNUM_CHILD\t{}\ndecision\t{}\nscore\t{}/{}\nARI\t{}\nrunningtime\t{}\nFPexistence\t{}\nFPindex\t{}\nTNindex\t{}\nmakeone_index\t{}\nfp_member_index\t{}".
                            format(NUM_CLONE_CLEMENT, NUM_CHILD_CLEMENT, "soft_" + priority, max_score, kwargs["NUM_MUTATION"],  ARI_CLEMENT_soft, round((datetime.datetime.now() - START_TIME).total_seconds()), cluster_soft.includefp_record [NUM_CLONE_soft[i]], cluster_soft.fp_index_record [NUM_CLONE_soft[i]], cluster_soft.tn_index_record [NUM_CLONE_soft[i]], cluster_soft.makeone_index_record [NUM_CLONE_soft[i]],  cluster_soft.fp_member_index_record [NUM_CLONE_soft[i]],   ), file = output_file)
                
                    

                #2. Save results
                subprocess.run (["cp -rf " + kwargs["CLEMENT_DIR"] + "/result/CLEMENT_soft_" + priority + ".results.txt  " + kwargs["CLEMENT_DIR"]  + "/CLEMENT_soft_" + priority + ".results.txt"], shell = True)
                for COPYTOPATH in [ kwargs["CLEMENT_DIR"], kwargs["CLEMENT_DIR"], kwargs["CLEMENT_DIR"] +"/result" ]:
                    pd.DataFrame(cluster_soft.membership_record [NUM_CLONE_soft[i]]).to_csv (COPYTOPATH + "/CLEMENT_soft_" + priority + ".membership.txt", index = False, header= False,  sep = "\t" )
                    pd.DataFrame(cluster_soft.mixture_record [NUM_CLONE_soft[i]],).to_csv ( COPYTOPATH + "/CLEMENT_soft_" + priority + ".mixture.txt", index = False, header= False,  sep = "\t" )
                    pd.DataFrame( np.unique( cluster_soft.membership_record [NUM_CLONE_soft[i]], return_counts = True ) ).to_csv ( COPYTOPATH + "/CLEMENT_soft_" + priority + ".membership_count.txt", index = False, header= False,  sep = "\t" )
                    pd.DataFrame( cluster_soft.membership_p_normalize_record [NUM_CLONE_soft[i]] ).to_csv ( COPYTOPATH + "/CLEMENT_soft_" + priority + ".posterior_allsample_normalize.txt", index = False, header= False,  sep = "\t" )            
                    if kwargs["SCORING"] == True:
                        pd.DataFrame(score_df).to_csv ( COPYTOPATH + "/CLEMENT_soft_" + priority + ".scoredf.tsv", index = False, header= True,  sep = "\t")
                    cluster_soft.df.to_csv ( COPYTOPATH + "/CLEMENT_soft_" + priority + ".df.tsv", index = False, header= False,  sep = "\t" )


                # 3.Visualization
                subprocess.run (["cp " +  kwargs["CLEMENT_DIR"] + "/candidate/clone" + str( NUM_CLONE_soft[i] ) + ".\(soft\)." + kwargs["IMAGE_FORMAT"] + " "  + kwargs["CLEMENT_DIR"]+ "/CLEMENT_soft_" + priority + "." + kwargs["IMAGE_FORMAT"]], shell = True)    



                #4. Phylogeny reconstruction
                ISPARENT = True if NUM_CLONE_CLEMENT  > NUM_CHILD_CLEMENT else False
                
                if ISPARENT == True:
                    print ("\n\n\n\n▶▶▶  PHYLOGENY ANALYSIS ◀◀◀")
                    kwargs["PHYLOGENY_DIR"] = kwargs["CLEMENT_DIR"] + "/CLEMENT_soft_" + priority + ".phylogeny.txt"
                    f = io.StringIO()
                    with contextlib.redirect_stdout(f):
                        membership_child = set ( cluster_soft.makeone_index_record[ NUM_CLONE_soft[i] ] )            # child의 번호만 뽑아준다 ( e.g.  0, 1, 3)
                        membership_outside = set (range (0, NUM_CLONE_soft [i] )) - membership_child - set ( [cluster_soft.fp_index_record[ NUM_CLONE_soft[i] ] ] )   # outlier의 번호만 뽑아준다 (e.g. 2)

                        g = phylogeny.main(membership_child, membership_outside, cluster_soft.mixture_record[ NUM_CLONE_soft[i] ],  **kwargs)

                    print ( f.getvalue() )
                    with open ( kwargs["PHYLOGENY_DIR"] , "w", encoding = "utf8") as phylogeny_file:
                        print (f.getvalue(), file = phylogeny_file)
                    subprocess.run (["cp -rf " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["CLEMENT_DIR"]  + "/CLEMENT_soft_" + priority + ".phylogeny.txt"], shell = True)
                    subprocess.run (["cp -rf " + kwargs["PHYLOGENY_DIR"] + "  " + kwargs["CLEMENT_DIR"]  + "/result/CLEMENT_soft_" + priority + ".phylogeny.txt"], shell = True)
                
                print ("\tSoft {} results printed".format(priority))







    if kwargs["MODE"] in ["Soft", "Both"]:
        #DECISION : Hard 인지 Soft인지
        print ("\n\n\n======================================= STEP #8. DECISION:  HARD VS SOFT  =======================================")

        if  ( NUM_CLONE_soft == []):     # Soft clustering에서 완전히 망했을 때
            DECISION = "hard_1st"
            soft_std = float("inf")
            print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\tNone\nSoft clustering에서 완전히 망해서 hard를 선택".format( DECISION,  "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] )    ))
            with open ( kwargs["CLEMENT_DIR"]+ "/CLEMENT_decision.evidence.txt"  , "w", encoding = "utf8") as output_file:
                print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\tNone\nSoft clustering에서 완전히 망해서 hard를 선택".format( DECISION, "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] )  ), file = output_file)


        elif (cluster_soft.mixture_record [ NUM_CLONE_hard[0] ] == []):  # hard에 해당되는 soft clustering이 없어서 발생했을 때
            if (kwargs ["MAKEONE_STRICT"] == 3) & (cluster_soft.mixture_record [ NUM_CLONE_soft[0] ] != []) & (cluster_soft.checkall_strict_record  [ NUM_CLONE_soft[0]  ] == True ) :  # BioData이고, best soft가 strict == True라면  soft clustering[0] 결과를 신뢰해주자
                DECISION = "soft_1st"
                hard_std, soft_std = float("inf"), float("inf")
                moved_col_list, moved_col_len = [], 0
                print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\t{}\n해당되는 Soft clustering이 없어서 soft clustering 기준으로 best를 선택".format( DECISION,  "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] ), "\t".join(str(np.round(row, 2)) for row in cluster_soft.mixture_record [NUM_CLONE_soft[0]] )    ))
                with open ( kwargs["CLEMENT_DIR"]+ "/CLEMENT_decision.evidence.txt"  , "w", encoding = "utf8") as output_file:
                    print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\t{}\n해당되는 Soft clustering이 없어서 soft clustering 기준으로 best를 선택".format( DECISION,  "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] ), "\t".join(str(np.round(row, 2)) for row in cluster_soft.mixture_record [NUM_CLONE_soft[0]]) ) , file = output_file)
            else: # 나머지라면 hard clustering 결과를 신뢰해주자
                DECISION = "hard_1st"
                hard_std, soft_std = float("inf"), float("inf")
                print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\tNone\n해당되는 Soft clustering이 없어서 hard를 선택".format( DECISION,  "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] )    ))
                with open ( kwargs["CLEMENT_DIR"]+ "/CLEMENT_decision.evidence.txt"  , "w", encoding = "utf8") as output_file:
                    print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\tNone\n해당되는 Soft clustering이 없어서 hard를 선택".format( DECISION,  "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] )  ), file = output_file)

        else :   # most of the cases
            if  DECISION == "soft(hard)_1st":
                print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\t{}\nsoft 선택기준 > {}\nposterior_sim\t{}".format( DECISION, "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] ), "\t".join(str(np.round(row, 2)) for row in cluster_soft.mixture_record [NUM_CLONE_hard[0]] ), round (kwargs["DECISION_STANDARD"], 2), round( posterior_sim, 2) ))
                with open ( kwargs["CLEMENT_DIR"]+ "/CLEMENT_decision.evidence.txt"  , "w", encoding = "utf8") as output_file:
                    print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\t{}\nsoft 선택기준 > {}\nposterior_sim\t{}".format( DECISION, "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] ), "\t".join(str(np.round(row, 2)) for row in cluster_soft.mixture_record [NUM_CLONE_hard[0]] ), round (kwargs["DECISION_STANDARD"], 2), round ( posterior_sim , 2) ), file = output_file)
            else:
                print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\t{}\nsoft 선택기준 > {}\nposterior_sim\t{}".format( DECISION, "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] ), "\t".join(str(np.round(row, 2)) for row in cluster_soft.mixture_record [NUM_CLONE_soft[0]] ), round (kwargs["DECISION_STANDARD"], 2), round( posterior_sim, 2) ))
                with open ( kwargs["CLEMENT_DIR"]+ "/CLEMENT_decision.evidence.txt"  , "w", encoding = "utf8") as output_file:
                    print ( "DECISION\t{}\nhard_mixture\t{}\nsoft_mixture\t{}\nsoft 선택기준 > {}\nposterior_sim\t{}".format( DECISION, "\t".join(str(np.round(row, 2)) for row in cluster_hard.mixture_record [NUM_CLONE_hard[0]] ), "\t".join(str(np.round(row, 2)) for row in cluster_soft.mixture_record [NUM_CLONE_soft[0]] ), round (kwargs["DECISION_STANDARD"], 2), round ( posterior_sim , 2) ), file = output_file)
                


        # 마지막 result 출력을 위한
        if "soft" in DECISION:
            NUM_CLONE_CLEMENT = NUM_CLONE_soft[i] 
            NUM_CHILD_CLEMENT = len (cluster_soft.makeone_index_record [NUM_CLONE_soft[0]])
            ARI_CLEMENT_DECISION = ARI_CLEMENT_soft
        elif "hard" in DECISION:
            NUM_CLONE_CLEMENT = NUM_CLONE_hard[i] 
            NUM_CHILD_CLEMENT = len (cluster_hard.makeone_index_record [NUM_CLONE_hard[0]])
            ARI_CLEMENT_DECISION = ARI_CLEMENT_hard
            

        # Copy results
        subprocess.run ([ "cp -rf \"" +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_" + DECISION + "." + kwargs["IMAGE_FORMAT"]  + "\" " + kwargs["CLEMENT_DIR"] + "/1.CLEMENT_decision." + kwargs["IMAGE_FORMAT"] ], shell = True)
        for COPYTOPATH in [ kwargs["CLEMENT_DIR"], kwargs["CLEMENT_DIR"] + "/result",  kwargs["CLEMENT_DIR"] +"/candidate" ]:
            subprocess.run ([ "cp -rf \"" +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_decision.evidence.txt\"" + " " + COPYTOPATH + "/CLEMENT_decision.evidence.txt" ], shell = True)
            subprocess.run ([ "cp -rf \"" +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_" + DECISION + "." + kwargs["IMAGE_FORMAT"]  + "\" " + COPYTOPATH + "/CLEMENT_decision." + kwargs["IMAGE_FORMAT"] ], shell = True)
        for COPYTOPATH in [ kwargs["CLEMENT_DIR"], kwargs["CLEMENT_DIR"] + "/result" ]:
            subprocess.run ([ "cp -rf \"" +  kwargs["CLEMENT_DIR"]+ "/result/CLEMENT_" + DECISION + ".results.txt\"" + " " + COPYTOPATH + "/CLEMENT_decision.results.txt" ], shell = True)
            subprocess.run ([ "cp -rf \"" +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_" + DECISION + ".membership.txt\"" + " " + COPYTOPATH + "/CLEMENT_decision.membership.txt" ], shell = True)
            subprocess.run ([ "cp -rf \"" +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_" + DECISION + ".membership_count.txt\"" + " " + COPYTOPATH + "/CLEMENT_decision.membership_count.txt" ], shell = True)
            subprocess.run ([ "cp -rf \"" +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_" + DECISION + ".mixture.txt\"" +  " " + COPYTOPATH + "/CLEMENT_decision.mixture.txt" ], shell = True)


if  __name__ == '__main__':
    main()