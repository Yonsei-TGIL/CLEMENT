import numpy as np
import pandas as pd
import os, subprocess
import datetime
import time
import argparse
import contextlib
import io

import filetype, datapreparation, phylogeny, EMhard, EMsoft, Estep, Mstep, Bunch, miscellaneous
import visualizationsingle, visualizationsinglesoft, visualizationpair, visualizationpairsoft

pd.options.mode.chained_assignment = None

kwargs = {}

parser = argparse.ArgumentParser(description='The below is usage direction.')
parser.add_argument('--INPUT_TSV', type=str, default="/data/project/Alzheimer/EM_cluster/EM_input/MRS_2_sample/M1-3_M1-8_input.txt",
                    help="Input data whether TSV or VCF. The tool automatically detects the number of samples")
parser.add_argument('--MODE', type=str, choices=["Hard", "Soft", "Both"], default="Both")
parser.add_argument('--NPVAF_DIR', default=None, help="Directory where selected datasets are")
parser.add_argument('--RANDOM_PICK', type=int, default=500,  help="The number of mutations that are randomly selected in each trials")
parser.add_argument('--USE_ALL', type=str, default="False",  help="If user want to make use of all the datasets, set True. If user want to downsize the sample, set False and --RANDOM_PICK a integer")
parser.add_argument('--CLEMENT_DIR', default=None, help="Directory where input and output of CLEMENT deposits")
parser.add_argument('--IMAGE_FORMAT', type=str, choices=["jpg", "pdf", "svg"], default="jpg", help="Output image format")


parser.add_argument('--NUM_CLONE_TRIAL_START', type=int, default=5, help="Minimum number of expected cluster_hards (initation of K)")
parser.add_argument('--NUM_CLONE_TRIAL_END', type=int, default=6, choices=range(1, 11), help="Maximum number of expected cluster_hards (termination of K)")
parser.add_argument('--NUM_CLONE_TRIAL_FORCE', type=int, default=4, help="Designate the number of expected cluster_hards by force")
parser.add_argument('--TRIAL_NO', default=3, type=int, choices=range(1, 21), help="Trial number in each candidate cluster_hard number. DO NOT recommend over 15")
parser.add_argument('--KMEANS_CLUSTERNO',  type=int, default=8, choices=range(5, 20), help="Number of initial K-means cluster_hard")
parser.add_argument('--MIN_CLUSTER_SIZE', type=int, default=15, help="The minimum cluster size that is acceptable")

parser.add_argument('--DEPTH_CUTOFF', default=100, type=int, help="The mutation of which depth below this values is abandoned")
parser.add_argument('--VERBOSE', type=int, choices=[0, 1, 2], default=2, help = "Print process →  0: no record,  1: simplified record,  2: verbose record")
parser.add_argument('--RANDOM_SEED', type=int, default=1, help="random_seed for regular random sampling")
parser.add_argument('--GAP', type = str, choices = ["new", "old", "sil"], default="new")

# Simulation dataset
# python3 CLEMENT.py --INPUT_TSV "/data/project/Alzheimer/EM_cluster/EM_input/simulation_2D/clone_4/2D_clone4_0.txt" --NPVAF_DIR "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/simulation_2D/clone_4/0"   --NPVAF_DIR "/data/project/Alzheimer/YSscript/EM_MRS/data/CLEMENT/simulation_2D/clone_4/0"   --NUM_CLONE_TRIAL_START 3 --NUM_CLONE_TRIAL_END 4 --DEPTH_CUTOFF 10 --VERBOSE 2 --TRIAL_NO 3  --RANDOM_SEED 1
# python3 CLEMENT.py --INPUT_TSV "/data/project/Alzheimer/EM_cluster/EM_input/simulation_3D/clone_4/3D_clone4_0.txt" --NPVAF_DIR "/data/project/Alzheimer/YSscript/EM_MRS/data/npvaf/simulation_3D/clone_4/0"  --NUM_CLONE_TRIAL_START 3 --NUM_CLONE_TRIAL_END 7 --DEPTH_CUTOFF 10 --VERBOSE 1 --TRIAL_NO 3  --RANDOM_SEED 1  --MODE Both


args = parser.parse_args()
kwargs["INPUT_TSV"] = args.INPUT_TSV

INPUT_TSV = kwargs["INPUT_TSV"]
INPUT_FILETYPE, NUM_BLOCK = filetype.main(INPUT_TSV)
kwargs["NUM_BLOCK_INPUT"], kwargs["NUM_BLOCK"] = NUM_BLOCK, NUM_BLOCK

SAMPLENAME = INPUT_TSV.split("/")[-1].split(".")[0]     # 'M1-5_M1-8_input'

kwargs["MODE"] = args.MODE
kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"], kwargs["NUM_CLONE_TRIAL_FORCE"] = args.NUM_CLONE_TRIAL_START, args.NUM_CLONE_TRIAL_END, args.NUM_CLONE_TRIAL_FORCE
kwargs["RANDOM_PICK"] = int(args.RANDOM_PICK)
if args.USE_ALL in ["True", "true", True]:
    kwargs["USE_ALL"] = True
elif args.USE_ALL in ["False", "false", False]:
    kwargs["USE_ALL"] = False

kwargs["IMAGE_FORMAT"] = str(args.IMAGE_FORMAT)
kwargs["TRIAL_NO"] = int(args.TRIAL_NO)
kwargs["DEPTH_CUTOFF"] = int(args.DEPTH_CUTOFF)
kwargs["VERBOSE"] = int(args.VERBOSE)
kwargs["MIN_CLUSTER_SIZE"] = int(args.MIN_CLUSTER_SIZE)
kwargs["KMEANS_CLUSTERNO"] = args.KMEANS_CLUSTERNO
kwargs["RANDOM_SEED"] = int(args.RANDOM_SEED)
kwargs["GAP"] = str(args.GAP)

kwargs["NUM_MUTATION"] = kwargs["RANDOM_PICK"]
NUM_MUTATION = kwargs["RANDOM_PICK"]

if args.NPVAF_DIR == None:
    kwargs["NPVAF_DIR"] = "/data/project/Alzheimer/YSscript/CLEMENT/data/npvaf/" + SAMPLENAME
else:
    kwargs["NPVAF_DIR"] = args.NPVAF_DIR
    
if args.CLEMENT_DIR == None:
    kwargs["CLEMENT_DIR"] = "/data/project/Alzheimer/YSscript/CLEMENT/data/CLEMENT/" + SAMPLENAME
else:
    kwargs["CLEMENT_DIR"] = args.CLEMENT_DIR



print ( "python3 EMhybrid.py  --INPUT_TSV {} --NPVAF_DIR {} --CLEMENT_DIR {}  --NUM_CLONE_TRIAL_START {}  --NUM_CLONE_TRIAL_END {}  --DEPTH_CUTOFF {}  --VERBOSE {}   --TRIAL_NO {}   --RANDOM_SEED {}   --MODE {}".
                format (kwargs["INPUT_TSV"], kwargs["NPVAF_DIR"], kwargs["CLEMENT_DIR"], kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"], kwargs["DEPTH_CUTOFF"], kwargs["VERBOSE"], kwargs["TRIAL_NO"], kwargs["RANDOM_SEED"], kwargs["MODE"] ) )

print("\nNOW RUNNING IS STARTED  :  {}h:{}m:{}s\n\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec)))


print("=============== STEP #1.   DATA EXTRACTION FROM THE INPUT  ===============")

os.system("mkdir -p " + kwargs["NPVAF_DIR"])
os.system("mkdir -p " + kwargs["CLEMENT_DIR"])
os.system("rm -rf " + kwargs["CLEMENT_DIR"] + "/trial")
os.system("mkdir -p " + kwargs["CLEMENT_DIR"] + "/trial")
os.system("rm -rf " + kwargs["CLEMENT_DIR"] + "/Kmeans")
os.system("mkdir -p " + kwargs["CLEMENT_DIR"] + "/Kmeans")
os.system("rm -rf " + kwargs["CLEMENT_DIR"] + "/candidate")
os.system("mkdir -p " + kwargs["CLEMENT_DIR"] + "/candidate")

inputdf, df, np_vaf, membership_answer, mixture_answer,  mutation_id, samplename_dict_input = datapreparation.main(**kwargs)
membership_answer_numerical = np.zeros(NUM_MUTATION, dtype="int")
membership_answer_numerical_nofp_index = []
if type(inputdf) != type(False):
    samplename_dict_input_rev = {v: k for k, v in samplename_dict_input.items()}   #  {0: 'FP', 1: 'V2', 2: 'S0', 3: 'V1'}

    print ( "samplename_dict_input : {}\nsamplename_dict_input_rev : {}".format(samplename_dict_input, samplename_dict_input_rev ))

    with open(kwargs["CLEMENT_DIR"] + "/0.input_membership_numerical.txt", "w", encoding="utf8") as result_answer:
        for k in range(NUM_MUTATION):
            membership_answer_numerical[k] = samplename_dict_input[membership_answer[k]]
            if ( membership_answer[k] != "FP"):
                membership_answer_numerical_nofp_index.append ( k )
            print(membership_answer_numerical[k], file=result_answer)
            
    with open(kwargs["CLEMENT_DIR"] + "/0.input_membership_letter.txt", "w", encoding="utf8") as result_answer:
        for k in range(NUM_MUTATION):
            print(membership_answer[k], file=result_answer)
            
    with open(kwargs["CLEMENT_DIR"] + "/0.input_mixture.txt", "w", encoding="utf8") as result_answer:
        print (mixture_answer, file=result_answer)
        for i in range (NUM_BLOCK):
            sum_mixture = 0
            for j in range (mixture_answer.shape[1]):
                if "FP" not in samplename_dict_input:
                    sum_mixture = sum_mixture + mixture_answer[i][j]
                else:
                    if samplename_dict_input ["FP"] != j:
                        sum_mixture = sum_mixture + mixture_answer[i][j]
            print ("Sum of mixture (except FP) in sample {} : {}".format(i, sum_mixture) , file=result_answer  )
        print (samplename_dict_input_rev ,  file=result_answer)
        print (np.unique(membership_answer_numerical, return_counts=True) [1], file=result_answer )
            

    if kwargs["NUM_BLOCK"] == 1:
        x_median = miscellaneous.VAFdensitogram(np_vaf, "INPUT DATA", kwargs["CLEMENT_DIR"] + "/0.inputdata_densitogram." + kwargs["IMAGE_FORMAT"], **kwargs)
        visualizationsingle.drawfigure_1d(membership_answer_numerical, "ANSWER_SET (n={})".format(kwargs["RANDOM_PICK"]), kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict_input_rev, False, -1)
    elif kwargs["NUM_BLOCK"] == 2:
        visualizationsingle.drawfigure_2d(membership_answer, "ANSWER_SET (n={})".format(kwargs["RANDOM_PICK"]), kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict_input, False, -1)
    elif kwargs["NUM_BLOCK"] >= 3:
        visualizationsingle.drawfigure_2d(membership_answer, "ANSWER_SET (n={})".format(kwargs["RANDOM_PICK"]), kwargs["CLEMENT_DIR"] + "/0.inputdata." + kwargs["IMAGE_FORMAT"], np_vaf, samplename_dict_input, False, -1, "SVD")
    
    

START_TIME = datetime.datetime.now()

print("NUM_SAMPLE : {}".format(NUM_BLOCK))
print ("NUM_MUTATION : {}\n\n\n".format( len(np_vaf)))


print ("\n\n--------------------------- [Step2-1. Hard clustering] ----------------------------------------------------")


kwargs["method"] = "gap+normal"
kwargs["adjustment"] = "half"


NUM_BLOCK, kwargs["NUM_BLOCK"]= len(df[0]), len(df[0])
NUM_MUTATION, kwargs["NUM_MUTATION"]  = kwargs["RANDOM_PICK"], kwargs["RANDOM_PICK"]
kwargs["STEP_NO"] = 30

np_vaf = miscellaneous.np_vaf_extract(df)
mixture_kmeans = miscellaneous.initial_kmeans (np_vaf, kwargs["KMEANS_CLUSTERNO"])

cluster_hard = Bunch.Bunch2(**kwargs)
cluster_soft = Bunch.Bunch2(**kwargs)

cluster_hard = EMhard.main (df, np_vaf, mixture_kmeans, **kwargs)





print ("\n\n--------------------------- [Step2-2. Soft clustering] ----------------------------------------------------\n\n")

for NUM_CLONE in range(kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
    print("\n#########################################################################\n\nNUM_CLONE = {0}".format(NUM_CLONE))
    kwargs["NUM_CLONE"], kwargs["NUM_CLONE_NOMINAL"] = NUM_CLONE, NUM_CLONE
    kwargs["OPTION"] = "soft"

    extinction = False

    if cluster_hard.likelihood_record[ NUM_CLONE ] !=  float("-inf"):
        print("\n\n\tSequential Soft clustering (TRIAL_NO = {}, STEP_NO = {})".format ( cluster_hard.trialindex_record[ NUM_CLONE ], cluster_hard.stepindex_record [ NUM_CLONE ] ))
        step_soft = Bunch.Bunch1(NUM_MUTATION , NUM_BLOCK, NUM_CLONE, cluster_hard.stepindex_record [ NUM_CLONE ] + kwargs["STEP_NO"])
        step_soft.copy (cluster_hard, 0, NUM_CLONE)  # Copy cluster_hard to cluster_soft_step 0 

        for step_index in range(1, kwargs["STEP_NO"]):   # Start from #1        ∵ 0 was filled with cluster_hard
            kwargs["STEP"], kwargs["TRIAL"] = step_index, cluster_hard.trialindex_record[ NUM_CLONE ]
            kwargs["STEP_TOTAL"] = step_index + cluster_hard.stepindex_record [ NUM_CLONE ] - 1

            step_soft = Estep.main(df, np_vaf, step_soft, **kwargs)                
            step_soft = Mstep.main(df, np_vaf, step_soft, "Soft", **kwargs) 

            print("\t\t{}th step ( = TOTAL {}th step) :  fp_index {}\tmakeone_index {}\tlikelihood : {}".format(kwargs["STEP"], kwargs["STEP_TOTAL"], step_soft.fp_index, step_soft.makeone_index, round (step_soft.likelihood, 1), step_soft.mixture ))

            step_soft.acc(step_soft.mixture, step_soft.membership, step_soft.likelihood, step_soft.membership_p, step_soft.membership_p_normalize, step_soft.makeone_index, step_soft.fp_index, step_index + 1, step_soft.outlier_index, step_soft.includeoutlier, step_index, step_index)

            if (miscellaneous.GoStop(step_soft, **kwargs) == "Stop")  :
                break
            if ( EMsoft.iszerocolumn (step_soft, **kwargs) == True) :
                print ("\t\t\t\t→ Terminated due to empty mixture\t{}".format(step_soft.mixture))
                break
            if ( len ( set (step_soft.membership) ) < NUM_CLONE ) :
                print ("\t\t\t\t→ Terminated due to empty clone")
                break


        step_soft.max_step_index =  step_soft.find_max_likelihood(1, step_soft.stepindex - 2 )   # Alleviating the condition that sum of mixture should be 1, because it is too stringent especially in 1-D decompoistion
        i = step_soft.max_step_index

        if (i == -1) | (step_soft.likelihood_record [i]  <= -9999999) | (extinction == True):   # If an error occured during the soft clustering iteration
            if i == -1: 
                print ("\t\t\tTerminated because because of inadequate beginning (1st step)")
            elif  (step_soft.likelihood_record [i]  <= -9999999) :
                print ("\t\t\tTerminated because there is no proper answer (satisfying sum of mixture = 1) in spite of itertating all the step")

        else:  # If no error occured & found the proper anse in Soft clustering
            print ("\t\tsoft max : {} → clone{}.{}-{}(soft).jpg to candidate." + kwargs["IMAGE_FORMAT"].format ( i,  NUM_CLONE, kwargs["TRIAL"] , i  + cluster_hard.stepindex_record [ NUM_CLONE ] - 1 ) )
            print ("\t\tmixture = {}".format ( step_soft.mixture_record [i] ))

            os.system ("cp " + kwargs["CLEMENT_DIR"] + "/trial/clone" + str (kwargs["NUM_CLONE"]) + "." + str( kwargs["TRIAL"] ) + "-"  + str(step_soft.max_step_index  + cluster_hard.stepindex_record [ NUM_CLONE ] - 1) + "\(soft\)." + kwargs["IMAGE_FORMAT"] + "  " + kwargs["CLEMENT_DIR"] + "/candidate/clone" + str (kwargs["NUM_CLONE"])  + ".\(soft\)." + kwargs["IMAGE_FORMAT"]  )
            cluster_soft.acc ( step_soft.mixture_record [i], step_soft.membership_record [i], step_soft.likelihood_record [i], step_soft.membership_p_record [i], step_soft.membership_p_normalize_record [i], step_soft.stepindex_record[i], cluster_hard.trialindex, step_soft.max_step_index_record[i], step_soft.makeone_index_record[i], step_soft.fp_index_record[i], step_soft.includeoutlier_record[i], step_soft.outlier_index_record[i]   ,**kwargs )


    else:   # No proper answer satisfying sum of mixture = 1 even in hard clustering
        print ("Terminated even in Hard clustering. Therefore, no further soft clustering would be")
        


print ("\n\nCURRENT TIME : {}h:{}m:{}s    (TIME LAPSED : {})\n\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME ))




NUM_CLONE_hard , NUM_CLONE_soft = [], []


print ("\n\n§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\nGap Statistics method (Hard clustering)")

f = io.StringIO()
with contextlib.redirect_stdout(f):
    NUM_CLONE_hard = miscellaneous.decision_gapstatistics (cluster_hard, np_vaf, **kwargs)
print ( f.getvalue() )
with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_hard.gapstatistics.txt", "w", encoding = "utf8") as gap_myEM:
    print (f.getvalue(), file = gap_myEM)


if NUM_BLOCK >= 2:
    print ("\n\n\n§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\nXieBeni index method (Soft clustering)")

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        NUM_CLONE_soft = miscellaneous.decision_XieBeni (cluster_soft, np_vaf, **kwargs)

    print ( f.getvalue() )
    with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.xiebeni.txt", "w", encoding = "utf8") as xiebeni_myEM:
        print (f.getvalue(), file = xiebeni_myEM)
        

elif NUM_BLOCK == 1:
    print ("\n\n\n§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§\nMax likelihood method (Soft clustering)")

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        NUM_CLONE_soft = miscellaneous.decision_max (cluster_soft, np_vaf, **kwargs)

    print ( f.getvalue() )
    with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_soft.maxlikelihood.txt", "w", encoding = "utf8") as maxlikelihood_myEM:
        print (f.getvalue(), file = maxlikelihood_myEM)







print ("\n\n\n\n=============== STEP #3.  DECISION :  Hard or Fuzzy  ===============")

DECISION = "hard_1st"
cluster_answer = cluster_hard
NUM_CLONE_answer = NUM_CLONE_hard

### Decision : Hard or Fuzzy
for i, priority in enumerate(["1st", "2nd"]):
    if (i >= len (NUM_CLONE_hard)) | ( cluster_hard.mixture_record [NUM_CLONE_hard[i]] == []):
        break

    if (priority == "1st"):
        moved_col_list = miscellaneous.movedcolumn ( cluster_hard.mixture_record [ NUM_CLONE_hard[i] ],  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] )
        hard_std = np.std(  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
        soft_std = np.std(  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ] [ : , moved_col_list ]   )
        if soft_std < hard_std * 0.8:
            DECISION = "soft_1st"
            cluster_answer = cluster_soft
            NUM_CLONE_answer = NUM_CLONE_soft
        
        
        print ( "DECISION : {}\t\thard_std : {}\tsoft_std : {}".format( DECISION, hard_std, soft_std ))
        
        with open (kwargs["CLEMENT_DIR"] + "/CLEMENT_decision.results.txt", "w", encoding = "utf8") as output_hard_vs_fuzzy:
            print ("moved column\t{}".format(moved_col_list), file = output_hard_vs_fuzzy)
            print ( "Hard (n = {}) : std = {}\tstep = {}\nhard_mixture = {}\n".format( cluster_hard.mixture_record [NUM_CLONE_hard[i]].shape[1],  hard_std,  cluster_hard.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_hard.mixture_record [ NUM_CLONE_hard[i] ]   ) , file = output_hard_vs_fuzzy  )
            print ( "Soft (n = {}) : std = {}\tstep = {}\nsoft_mixture = {}".format( cluster_soft.mixture_record [NUM_CLONE_hard[i]].shape[1],  soft_std,  cluster_soft.stepindex_record [ NUM_CLONE_hard[i] ],  cluster_soft.mixture_record [ NUM_CLONE_hard[i] ]   )  , file = output_hard_vs_fuzzy )
            print ("DECISION\t{}".format(DECISION) , file = output_hard_vs_fuzzy )



###### Print results ######

pd.DataFrame(cluster_answer.membership_record [NUM_CLONE_answer[i]]).to_csv (kwargs["CLEMENT_DIR"] + "/CLEMENT_decision.membership.txt", index = False, header= False,  sep = "\t" )
pd.DataFrame(cluster_answer.mixture_record [NUM_CLONE_answer[i]],).to_csv (kwargs["CLEMENT_DIR"] + "/CLEMENT_decision.mixture.txt", index = False, header= False,  sep = "\t" )
subprocess.run (["cp " +  kwargs["CLEMENT_DIR"]+ "/candidate/clone" + str(NUM_CLONE_answer[i]) + ".(" + DECISION.split("_")[0] + ")." + kwargs["IMAGE_FORMAT"] + "  " +  kwargs["CLEMENT_DIR"]+ "/CLEMENT_" + DECISION +"." + kwargs["IMAGE_FORMAT"]], shell = True)

        
        
        

print ("\n\n\n\n=============== STEP #4.  PHYLOGENY RECONSTRUCTION  ===============")

if "hard" in DECISION:
    if len (cluster_hard.makeone_index_record [NUM_CLONE_hard[i]]) + int ( cluster_hard.includeoutlier_record [NUM_CLONE_hard[i]])  < NUM_CLONE_hard[i]:    # FP
        ISPARENT = True
    else:
        ISPARENT = False


    if ISPARENT == True:
        kwargs["PHYLOGENY_DIR"] = kwargs["CLEMENT_DIR"] + "/CLEMENT_hard_" + priority + ".phylogeny.txt"
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            membership_child = set ( cluster_hard.makeone_index_record[ NUM_CLONE_hard[i] ] )            # Pick makeone_index (= child_index,  e.g.,  0, 1, 3)
            membership_outside = set (range (0, NUM_CLONE_hard [i] )) - membership_child - set ( [cluster_hard.fp_index_record[ NUM_CLONE_hard[i] ] ] )   # Pick outlier_index  (= putative parent index, e.g., 2) 

            g = phylogeny.main(membership_child, membership_outside, cluster_hard.mixture_record[ NUM_CLONE_hard[i] ],  **kwargs)

        print ( f.getvalue() )
        with open ( kwargs["PHYLOGENY_DIR"] , "w", encoding = "utf8") as phylogeny_file:
            print (f.getvalue(), file = phylogeny_file)

elif "soft" in DECISION:
    if len (cluster_soft.makeone_index_record [NUM_CLONE_soft[i]]) + int ( cluster_soft.includeoutlier_record [NUM_CLONE_soft[i]])  < NUM_CLONE_soft[i]:    # FP
        ISPARENT = True
    else:
        ISPARENT = False

    if ISPARENT == True:
        kwargs["PHYLOGENY_DIR"] = kwargs["CLEMENT_DIR"] + "/CLEMENT_soft_" + priority + ".phylogeny.txt"
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            membership_child = set ( cluster_soft.makeone_index_record[ NUM_CLONE_soft[i] ] )            # Pick makeone_index (= child_index,  e.g.,  0, 1, 3)
            membership_outside = set (range (0, NUM_CLONE_soft [i] )) - membership_child - set ( [cluster_soft.fp_index_record[ NUM_CLONE_soft[i] ] ] )   # Pick outlier_index  (= putative parent index, e.g., 2) 

            g = phylogeny.main(membership_child, membership_outside, cluster_soft.mixture_record[ NUM_CLONE_soft[i] ],  **kwargs)

        print ( f.getvalue() )
        with open ( kwargs["PHYLOGENY_DIR"] , "w", encoding = "utf8") as phylogeny_file:
            print (f.getvalue(), file = phylogeny_file)



print ("\n\nCURRENT TIME : {}h:{}m:{}s    (TIME LAPSED : {})\n\n".format(time.localtime().tm_hour, time.localtime().tm_min, round(time.localtime().tm_sec), datetime.datetime.now() - START_TIME )) 
