import numpy as np
import scipy
import pandas as pd

def firstapperance  (score_df, i):
    answer_firstappearance, myanswer_firstappearance = False, False
    if np.where ( np.array(score_df["answer"]) == score_df.iloc[i]["answer"]) [0][0] == i:
        answer_firstappearance = True
    if np.where ( np.array(score_df["predicted"]) == score_df.iloc[i]["predicted"]) [0][0] == i:
        myanswer_firstappearance = True
    
    return answer_firstappearance, myanswer_firstappearance


def Yindex (score_df, coefficient = 2):
    distance = np.zeros (score_df.shape[0], dtype = "float")
    weight = np.zeros(score_df.shape[0], dtype = "int")

    for j in range (0, score_df.shape[0]):
        distance[j] =  (score_df.iloc[j]["distance"])

        if firstapperance (score_df, j) == (True, True):
            weight[j] = (score_df.iloc[j]["shared"])
        else:
            if firstapperance (score_df, j) == (False, True):     # 내 clone이 정답지에 없을 때  
                weight[j] =  ( score_df.iloc[j]["n(predicted)"] * coefficient)    
                #print ("{} ({} {}) : 내 clone이 정답지에 없어서 다른 정답과 매칭함".format(j, score_df.iloc[j]["answer"],  score_df.iloc[j]["predicted"]))
            else:                                                 # 내가 정답지에 있는 clone을 못 맞혔을 때
                weight[j] =  ( score_df.iloc[j]["n(answer)"] * coefficient)    
                #print ("{} ({} {}) : 정답지에 있는 clone을 못 찾았음".format(j, score_df.iloc[j]["answer"],  score_df.iloc[j]["predicted"]))

    try:
        Y_index = round(np.average(distance, weights=weight), 3)
    except:
        Y_index = 999

    return Y_index


def ARI (membership_answer_numerical, membership_predicted):
    from sklearn.metrics.cluster import adjusted_rand_score
    return adjusted_rand_score ( membership_answer_numerical, membership_predicted )


def fARI ( membership_answer_numerical, membership_p_predicted, temp_dir, SCRIPT_DIR, transform):
    import subprocess
    import pandas as pd
    import numpy as np

    with open( temp_dir + "/membership_answer_temp.tsv", "w") as file:
        file.write( "\t".join( [str(item) for item in membership_answer_numerical]  )  )

    if transform == False: # Soft clustering의 matrix를 그대로 받았던 경우
        pd.DataFrame (membership_p_predicted).to_csv (temp_dir + "/membership_p_predicted_temp.tsv",  index = True, header= True,  sep = "\t")
    else:  # Hard clustering인데 fuzzy matrix로 만들어줘야 하는 경우
        mat = np.zeros (  ( len(membership_answer_numerical) , np.max (membership_p_predicted) + 1 ),  dtype = "float")  
        # print ( np.unique (membership_p_predicted, return_counts = True ) )
        # print (mat.shape)
        # print ( sorted(set(membership_p_predicted) ) )
        for k in range (len (membership_p_predicted)):
            mat[k, membership_p_predicted [k] ] = 1
        pd.DataFrame (mat).to_csv (temp_dir + "/membership_p_predicted_temp.tsv",  index = True, header= True,  sep = "\t")


    # Run
    subprocess.run (["Rscript " + SCRIPT_DIR + "/fARI.R  " + temp_dir + "/membership_answer_temp.tsv "  + temp_dir + "/membership_p_predicted_temp.tsv " + temp_dir + "/fARI_temp.txt " +  str(transform) ], shell = True)

    # Read ARI
    with open(temp_dir + "/fARI_temp.txt", "r") as file:
        fARI = [float(line.strip()) for line in file]

    return fARI[0]

##############################
def FPmatrix (score_df):
    if "FP" not in list (score_df["answer"]):
        return 0, 0, 0, None, None, None

    FP_df = score_df [score_df["answer"] == "FP"].reset_index().iloc[0]       # 무조건 맨 위를 잡는다



    try:
        sensitivity =  round( int (FP_df["shared"])  / int (FP_df["n(answer)"]) , 2)
    except:
        sensitivity = None
    try:
        PPV = round( int (FP_df["shared"])  / int (FP_df["n(predicted)"]) , 2)
    except:
        PPV  = None
    try:
        F1 = round( 2 * (sensitivity*PPV) / (sensitivity + PPV), 2)
    except:
        F1 = None

    #print ("\nSensitivity : {}\nPPV : {}\nF1 : {}".format(sensitivity, PPV, F1))

    return  int (FP_df["n(answer)"]) - int (FP_df["shared"]),  int (FP_df["shared"]),  int (FP_df["n(predicted)"])  -   int(FP_df["shared"]), sensitivity, PPV, F1