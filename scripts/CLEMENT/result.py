import numpy as np
import scipy
import pandas as pd

def firstapperance  (score_df, i):
    answer_firstappearance, myanswer_firstappearance = False, False
    if np.where ( np.array(score_df["answer"]) == score_df.iloc[i]["answer"]) [0][0] == i:
        answer_firstappearance = True
    if np.where ( np.array(score_df["myanswer"]) == score_df.iloc[i]["myanswer"]) [0][0] == i:
        myanswer_firstappearance = True
    
    return answer_firstappearance, myanswer_firstappearance


def Yindex (score_df, coefficient = 2):
    distance = np.zeros (score_df.shape[0], dtype = "float")
    weight = np.zeros(score_df.shape[0], dtype = "int")

    for j in range (0, score_df.shape[0]):
        distance[j] =  (score_df.iloc[j]["distance"])

        if firstapperance (score_df, j) == (True, True):
            weight[j] = (score_df.iloc[j]["num"])
        else:
            if firstapperance (score_df, j) == (False, True):     # 내 clone이 정답지에 없을 때  
                weight[j] =  ( score_df.iloc[j]["n(myanswer)"] * coefficient)    
                #print ("{} ({} {}) : 내 clone이 정답지에 없어서 다른 정답과 매칭함".format(j, score_df.iloc[j]["answer"],  score_df.iloc[j]["myanswer"]))
            else:                                                 # 내가 정답지에 있는 clone을 못 맞혔을 때
                weight[j] =  ( score_df.iloc[j]["n(answer)"] * coefficient)    
                #print ("{} ({} {}) : 정답지에 있는 clone을 못 찾았음".format(j, score_df.iloc[j]["answer"],  score_df.iloc[j]["myanswer"]))

    try:
        Y_index = round(np.average(distance, weights=weight), 3)
    except:
        Y_index = 999

    return Y_index


def ARI (membership_answer_numerical, membership_tool):
    from sklearn.metrics.cluster import adjusted_rand_score

    return adjusted_rand_score ( membership_answer_numerical, membership_tool )


def FPmatrix (score_df):
    if "FP" not in list (score_df["answer"]):
        return 0, 0, 0, None, None, None

    FP_df = score_df [score_df["answer"] == "FP"].reset_index().iloc[0]       # 무조건 맨 위를 잡는다



    try:
        sensitivity =  round( int (FP_df["num"])  / int (FP_df["n(answer)"]) , 2)
    except:
        sensitivity = None
    try:
        PPV = round( int (FP_df["num"])  / int (FP_df["n(myanswer)"]) , 2)
    except:
        PPV  = None
    try:
        F1 = round( 2 * (sensitivity*PPV) / (sensitivity + PPV), 2)
    except:
        F1 = None

    print ("\nSensitivity : {}\nPPV : {}\nF1 : {}".format(sensitivity, PPV, F1))

    return  int (FP_df["n(answer)"]) - int (FP_df["num"]),  int (FP_df["num"]),  int (FP_df["n(myanswer)"])  -   int(FP_df["num"]), sensitivity, PPV, F1