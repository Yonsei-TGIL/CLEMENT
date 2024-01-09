import numpy as np
import itertools
import pandas as pd
pd.options.mode.chained_assignment = None


def Sorting(membership):

    NUM_MUTATION = len(membership)

    # membership을 앞에 나오는 순서대로 numbering을 다시 해줌
    membership_sort = []
    dict = {}
    shared = 0

    # membership을 정렬해서 보내줌
    for k in range(NUM_MUTATION):
        if membership[k] not in dict:
            membership_sort.append(shared)
            dict[membership[k]] = shared
            shared = shared + 1
        else:
            membership_sort.append(dict[membership[k]])
    membership_sort = np.array(membership_sort)

    return membership_sort


def SoftScoring (membership_answer_character, membership_answer_numerical, cluster_soft, NUM_CLONE):
    ANSWER_NUM_CLONE = int(np.max(membership_answer_numerical) + 1)    # 정답set 의 clone
    NUM_MUTATION = cluster_soft.membership_p_normalize_record [ NUM_CLONE ].shape[0]
    # mixture 는 NUM_CLONE + 1 (FP 포함)

    # Soft clustering에서 나온 Clone No >  ANSWER_NUM_CLONE일 때
    if NUM_CLONE  >= ANSWER_NUM_CLONE:      # 이러면 answer membership을 뒤섞어보면 된다
        max_index_dict, max_score = {}, -1
        
        for tt in list(itertools.permutations( list(range(0, NUM_CLONE)),  ANSWER_NUM_CLONE ) ):
            index_dict = {}    # answer -> CLEMENT
            for i in range(0, np.max(membership_answer_numerical) + 1 ):
                index_dict[i] = tt[i]

            # answer membership을 돌려보기
            temp = []
            for i in range( len(membership_answer_numerical) ):
                temp.append(index_dict[membership_answer_numerical[i]])

            # maxmaxmax_membership과의 일치도 확인
            sum = 0
            for k in range( NUM_MUTATION ):
                sum += cluster_soft.membership_p_normalize_record  [NUM_CLONE] [k] [temp [k]]

            if sum > max_score:
                max_score = sum
                max_temp = temp
                max_index_dict = index_dict
                print ( "index_dict = {}\tmax_score = {}".format (index_dict, max_score))


        sample_dict_PtoA = {}   # 
        sample_dict_AtoP = {}   # 변환표 (R)
        for i in range( len(max_temp)):
            sample_dict_PtoA [max_temp[i]] = membership_answer_character[i]   # {3: 'MRC5_het', 4: 'FP', 0: 'V3_het', 1: 'V5_het', 2: 'V1_het'}
            sample_dict_AtoP [membership_answer_character[i]] = max_temp[i]




def Scoring(membership_answer_character, membership_answer_numerical, my_membership, fp_index, parent_index, **kwargs):
    global ANSWER_NUM_CLONE
    max_score = 0
    max = []

    if ("FP" not in membership_answer_character):
        FP_MATCH = None
    if ("FP" in membership_answer_character) & (fp_index != -1):   # FP를 제대로 designate 했을 경우
        FP_MATCH = True
    elif ("FP" in membership_answer_character) & (fp_index == -1):  # FP가 있는데 제대로 designate 하지 못했을 경우
        FP_MATCH = False


    # Celldata의 경우
    
    PARENT_MATCH = None
    for i in range ( len  (membership_answer_character) ):
        if "," in membership_answer_character[i]:
            PARENT_MATCH = True

    if PARENT_MATCH != None:
        if ( len(parent_index) == 0 ):  # parent가 있는데 제대로 designate 하지 못했을 경우
            PARENT_MATCH = False

    print ("membership_answer_numerical = {}".format ( np.unique ( membership_answer_numerical, return_counts = True  )) )

    ANSWER_NUM_CLONE = len ( set ( membership_answer_numerical)  ) # 정답set 의 clone 개수 (FP 포함)

    if np.max(my_membership)  >= np.max(membership_answer_numerical):      # Predicted cluster no가 더 많을 때에는 membership_answer를 돌린다
        max_index_dict = {}
        
        for tt in list(itertools.permutations( list(range(0, np.max(my_membership) + 1)),  int ( np.max(membership_answer_numerical) + 1) ) ):
            index_dict, index_dict_rev = {}, {}
            for i in range(0, np.max(membership_answer_numerical) + 1 ):
                index_dict[i] = tt[i]
                index_dict_rev[ tt[i] ] = i

            # answer membership을 돌려보기
            temp = []
            for i in range(len(membership_answer_numerical)):
                temp.append(index_dict[membership_answer_numerical[i]])

            # answer가 FP인데 designate 하지 못한 경우, 일부러 99로 만들어서 어긋나게 만들어준다
            if FP_MATCH == False:
                for k in range(len(membership_answer_numerical)):
                    if membership_answer_character[k] == "FP":
                        index_dict [temp[k]] = 99
                        index_dict_rev [ 99 ] = temp[k]
                        temp[k] = 99
                        kwargs["samplename_dict_NumToCharacter"][99] = "FP"
            if PARENT_MATCH == False:
                for k in range(len(membership_answer_numerical)):
                    if "," in membership_answer_character[k]:
                        index_dict [temp[k]] = 999
                        index_dict_rev [ 999 ] = temp[k]
                        temp[k] = 999
                        kwargs["samplename_dict_NumToCharacter"][999] = "FP"

            # maxmaxmax_membership과의 일치도 확인
            if np.sum(np.equal(temp, my_membership)) > max_score:
                max_score = np.sum(np.equal(temp, my_membership))
                max = temp
                max_index_dict = index_dict
                kkk =  np.unique (   [ temp[k] for k in np.array ( np.where ( np.equal ( temp , my_membership ) ) [0] )  ]   , return_counts = True)    # 공통인 것들의 index, count
                # print ("temp = {}".format ( np.unique (temp , return_counts = True )))
                # print ("index_dict = {}\tindex_dict_rev = {}".format (index_dict, index_dict_rev))
                # print ("\tmy_membership = {}".format ( np.unique ( my_membership, return_counts = True  )) )
                # print ( "\tkkk = {}".format ( kkk) ) 
                score_df = pd.DataFrame ( {"answer" : [  kwargs["samplename_dict_NumToCharacter"][ index_dict_rev[k] ] for k in kkk[0] ] ,
                                                                    "predicted" : [ k for k in kkk[0] ] , 
                                                                    "shared" : kkk[1],
                                                                    "n(answer)" : [ np.unique ( temp, return_counts = True ) [1][k] for k in  [ index_dict_rev[k] for k in kkk[0] ] ],      # 돌린다
                                                                    "n(predicted)" : [ np.unique ( my_membership, return_counts = True ) [1][k] for k in kkk[0]  ]    } )   # 가만히 있음


        sample_dict_PtoA = {}   #
        sample_dict_AtoP = {}   # 변환표 (R)
        for i in range(len(max)):
            sample_dict_PtoA [max[i]] = membership_answer_character[i]   # {3: 'MRC5_het', 4: 'FP', 0: 'V3_het', 1: 'V5_het', 2: 'V1_het'}
            sample_dict_AtoP [membership_answer_character[i]] = max[i]



    else:   # Answer cluster no가 더 많을 때에는 my_membership를 돌린다
        max_index_dict = {}
        for tt in list(itertools.permutations( list(range(0, np.max(membership_answer_numerical) + 1)),  int ( np.max(my_membership) + 1) ) ):
            index_dict, index_dict_rev = {}, {}
            for i in range(0, np.max(my_membership) + 1):    # predicted -> answer
                index_dict[i] = tt[i]
                index_dict_rev[ tt[i] ] = i


            # predicted membership을 돌려보기
            temp = []
            for k in range(len(my_membership)):
                temp.append(index_dict[my_membership[k]])

            # answer가 FP인데 designate 하지 못한 경우, 일부러 99로 만들어서 어긋나게 만들어준다
            if FP_MATCH == False:
                for k in range(len(membership_answer_numerical)):
                    if membership_answer_character[k] == "FP":
                        temp[k] = 99
                        index_dict [temp[k]] = 99
                        index_dict_rev [ 99 ] = temp[k]
                        kwargs["samplename_dict_NumToCharacter"][99] = "FP"
            if PARENT_MATCH == False:
                for k in range(len(membership_answer_numerical)):
                    if "," in membership_answer_character[k]:
                        temp[k] = 999
                        index_dict [temp[k]] = 999
                        index_dict_rev [ 999 ] = temp[k]
                        kwargs["samplename_dict_NumToCharacter"][999] = "FP"
            
            # maxmaxmax_membership과의 일치도 확인
            if np.sum(np.equal(temp, membership_answer_numerical)) > max_score:
                max_score = np.sum(np.equal(temp, membership_answer_numerical))
                max = temp
                max_index_dict = index_dict
                kkk =  np.unique (   [ temp[k] for k in np.array ( np.where ( np.equal ( temp , membership_answer_numerical ) ) [0] )  ]   , return_counts = True)    # 공통인 것들의 index, count

                score_df = pd.DataFrame ( {"answer" : [ kwargs["samplename_dict_NumToCharacter"] [k] for k in kkk[0] ] ,
                                                    "predicted" : [ index_dict_rev[k]  for k in kkk[0] ] , 
                                                    "shared" : kkk[1],
                                                    "n(answer)" : [ np.unique ( membership_answer_numerical, return_counts = True ) [1][k] for k in kkk[0] ],
                                                    "n(predicted)" : [ np.unique ( temp, return_counts = True ) [1][k] for k in  [ index_dict_rev[k] for k in kkk[0] ]  ]    } ) 

        sample_dict_PtoA = {}          # predicted -> answer  {0: 'V1', 1: 'V2', 2: 'FP', 3: 'S0', 4: 'S0,V2'}
        sample_dict_AtoP = {}                 # answer -> predicted  {'V1': 0, 'V2': 1, 'FP': 2, 'S0': 3, 'S0,V1': 4, 'S0,V2': 1}
        for key, value in max_index_dict.items() :
            answer_char = ""
            for k in range ( len( membership_answer_character) ):
                if membership_answer_numerical [k] == value:   # 해당 answer numerical value가 character로는 뭔지 궁금
                    answer_char = membership_answer_character [k]
                    break
            sample_dict_PtoA[ key  ] = answer_char
            sample_dict_AtoP[ answer_char ] = key
            
    #print ("sample_dict_PtoA = {}\nsample_dict_AtoP = {}".format(sample_dict_PtoA, sample_dict_AtoP))

    return max_score, sample_dict_PtoA, sample_dict_AtoP, score_df


def find_mindistance(j, mixture_my, mixture_target):
    import scipy
    distance_temp = []
    mixture_my[:, j][np.isnan(mixture_my[:, j])] = 0   # FP가 nan으로 뜰 수가 있는데 그럴 떄에는 0으로 바꿔줘야 한다
    for k in range(mixture_target.shape[1]):
        try:
            distance_temp.append(scipy.spatial.distance.euclidean( mixture_my[:, j], mixture_target[:, k]))
        except:
            print ( mixture_my[:, j] )
            print ( mixture_target[:, k] )

    #print (distance_temp)
    mindist, minarg = np.min(distance_temp), np.argmin(distance_temp)
    return (mindist, minarg)


def mixturebased(mixture_left, mixture_right, membership_left, membership_right, samplename_dict_input, samplename_dict_input_rev, includeoutlier, fp_index, tool, **kwargs):
    import pandas as pd
    import scipy

    membership_left, membership_right = list(membership_left), list(membership_right)
    intersect_list = []

    # 우선 Greedy하게 가장 죽이 잘 맞는 화살표를 만들어준다

    # FP가 answer set에 있고,  내 tool도 FP를 불렀으면 강제 매칭
    FP_MATCH = False  # 기본적으로 FP match가 없다고 생각
    # Outlier일 경우 FP와 강제로 matching시켜준다
    if (includeoutlier == True) & ("FP" in samplename_dict_input.keys()):
        j_index = np.where(np.array(membership_right) == fp_index)[0]
        k_index = np.where(np.array(membership_left) == "FP")[0]
        intersect_list.append(["FP", fp_index, len(set(j_index) & set(k_index))])
        FP_MATCH = True

    #print ("includeoutlier = {}\tfp_index = {}\tFP_MATCH = {}".format (includeoutlier, fp_index, FP_MATCH))

    for j in range(mixture_right.shape[1]):
        if (FP_MATCH == True) & (j == fp_index):  # 위에서 했으니 넘어가자  (FP 우선matching)
            continue
        else:       # 대부분의 경우에
            intersect_temp = []
            # 주어진 answer (S0, V1 등)에 대해 가장 잘 맞는 predicted(0, 3, 4, 2,1)를 찾기
            for k in range(mixture_left.shape[1]):
                j_index = np.where(np.array(membership_right) == j)[0]
                k_index = np.where(np.array(membership_left) == samplename_dict_input_rev[k])[
                    0]  # {0: 'FP', 1: 'V2', 2: 'S0', 3: 'V1'}
                intersect_temp.append(len(set(j_index) & set(k_index)))
            maxintersect, maxintersectarg = np.max( intersect_temp), np.argmax(intersect_temp)
            intersect_list.append(  [samplename_dict_input_rev[maxintersectarg], j, maxintersect])

        # if (includeoutlier == "Yes") & (j == mixture_right.shape[1] - 1) & ("FP" in samplename_dict_input.keys()) :  #Outlier일 경우 FP와 강제로 matching시켜준다
            # j_index = np.where ( np.array(membership_right) == j) [0]
            # k_index = np.where ( np.array(membership_left) == "FP" ) [0]
            # intersect_list.append ( [ "FP", j, len( set(j_index) & set(k_index)) ]  )
        # else:
        #     intersect_list.append( [ samplename_dict_input_rev [maxintersectarg]  , j , maxintersect] )

    if includeoutlier == True:
        score_df = pd.concat([pd.DataFrame([intersect_list[0]], columns=["answer", "predicted",  "shared"]), pd.DataFrame(intersect_list[1:], columns=[
                             "answer", "predicted",  "shared"]).sort_values(by="shared", axis=0, ascending=False).reset_index(drop=True)], ignore_index=True)
    else:
        score_df = pd.DataFrame(intersect_list, columns=["answer", "predicted",  "shared"]).sort_values(
            by="shared", axis=0, ascending=False).reset_index(drop=True)
    score_df["distance"] = 0.0
    score_df = score_df.astype({'distance': 'float'})

    p = []
    for i in range(score_df.shape[0]):
        if score_df.iloc[i]["answer"] not in p:      # 일대일 매칭이 최우선으로 된 경우
            score_df["distance"].iloc[i] = float(scipy.spatial.distance.euclidean( mixture_right[:, score_df.iloc[i]["predicted"]], mixture_left[:,  samplename_dict_input[score_df.iloc[i]["answer"]]]))
        else:      # 이미 answer(left)가 다른 짝이 있어서 혼자 남은 불쌍한 predicted(right)
            # 1. 남은 애들 중 그 다음 짝을 찾아준다
            j = score_df.iloc[i]["predicted"]
            if kwargs["VERBOSE"] >= 3:
                print("{}는 혼자 남았습니다 남은 것 \t\t{}".format(j,  set(samplename_dict_input.keys()) - set(score_df.iloc[:i]["answer"])))
            remarry_list, intersect_temp = [], []
            for unmet_answer in set(samplename_dict_input.keys()) - set(score_df.iloc[:i]["answer"]):
                j_index = np.where(np.array(membership_right) == j)[0]
                k_index = np.where(np.array(membership_left) == unmet_answer)[0]
                intersect_temp.append(len(set(j_index) & set(k_index)))
                remarry_list.append(unmet_answer)
                if kwargs["VERBOSE"] >= 3:
                    print("{} - {} : 겹치는 것 {}".format(j, unmet_answer,  len(set(j_index) & set(k_index))))
            try:
                maxintersect, maxintersectarg = np.max(
                    intersect_temp), np.argmax(intersect_temp)
            except:
                maxintersect = 0

            if maxintersect > 0:  # 겹치는 게 뭐라도 있을 경우 재혼이라도 시켜준다
                score_df["answer"].iloc[i] = remarry_list[maxintersectarg]
                if kwargs["VERBOSE"] >= 3:
                    print("불쌍한 {}와 불쌍한 {}가 만나 다시 재혼했습니다.  (겹치는 것 : {})".format( j,  remarry_list[maxintersectarg], maxintersect))
                score_df["shared"].iloc[i] = maxintersect
                score_df["distance"].iloc[i] = float(scipy.spatial.distance.euclidean(mixture_right[:, score_df.iloc[i]["predicted"]],
                                                     mixture_left[:, samplename_dict_input[score_df["answer"].iloc[i]]]))     # 혼자 남은 불쌍한 predicted(right)  - 버려졌다가 구햊딘 answer(left)

            # 전혀 겹치는게 없으면 그냥 포기하고 가장 mixture 값이 가까운 애로 찾아준다  (걔가 이미 짝이 있어도 어쩔 수 없음)
            elif maxintersect == 0:
                #print ("전혀 겹치는 게 없다 : j = {}\tmixture_right = {}\tmixture_left = {}".format (score_df.iloc[i]["predicted"], mixture_right, mixture_left ))
                mindist, minarg = find_mindistance( score_df.iloc[i]["predicted"],  mixture_right, mixture_left)
                j = score_df.iloc[i]["predicted"]
                k = samplename_dict_input_rev[minarg]

                if kwargs["VERBOSE"] >= 3:
                    print("불쌍한 {}는 결국 짝을 찾지 못하고 {}만 바라보면서 거리를 계산합니다.".format(j, k))

                j_index = np.where(np.array(membership_right) == j)[0]
                k_index = np.where(np.array(membership_left) == k)[0]
                #겹치는 게 없으니 0을 줄 수 밖에...
                score_df["answer"].iloc[i], score_df["shared"].iloc[i], score_df["distance"].iloc[i] = k,  0, mindist
                #score_df["answer"].iloc[i], score_df["shared"].iloc[i], score_df["distance"].iloc[i] = 0,  0, 0

        p.append(score_df.iloc[i]["answer"])


# 선택받지 못한 answer(left)를 매칭시켜준다
    unmet_list = []
    for unmet_answer in set(samplename_dict_input.keys()) - set(score_df["answer"]):
        if kwargs["VERBOSE"] >= 3:
            print("아직 선택받지 못한 {}".format(unmet_answer))
        i = samplename_dict_input[unmet_answer]
        if i >= mixture_left.shape[1]:
            continue
        if mixture_right.shape[1] == 0:
            continue

        mindist, minarg = find_mindistance(i,  mixture_left, mixture_right)
        unmet_list.append([unmet_answer, minarg,  0, mindist])

    score_df = score_df.append(pd.DataFrame(unmet_list, columns=["answer", "predicted", "shared", "distance"])).reset_index(drop=True)

    # membership 개수를 추가시켜주고 마무리
    score_df["n(answer)"] = 0
    score_df["n(predicted)"] = 0

    score = 0

    for i in range(score_df.shape[0]):
        score_df["n(answer)"].iloc[i] = membership_left.count( score_df.iloc[i]["answer"])
        score_df["n(predicted)"].iloc[i] = membership_right.count( score_df.iloc[i]["predicted"])
        if score_df.iloc[i]["shared"] == 0:
            if i < mixture_right.shape[1]:
                score_df["n(answer)"].iloc[i] = 0    # 매칭받지 못했던 애에게는 0을 줘야 한다
            else:
                # 매칭받지 못했던 애에게는 0을 줘야 한다
                score_df["n(predicted)"].iloc[i] = 0

        score = score + score_df.iloc[i]["shared"]
        if (fp_index == -1):
            # FP를 제대로 designate 못한 경우에는 score를 강제로 박탈시켜준다
            if (score_df["answer"].iloc[i] == "FP"):
                score = score - score_df.iloc[i]["shared"]

    score_df["distance"] = score_df["distance"].round(3)

    return score_df, score


def main(membership_answer, maxmaxmax_membership, **kwargs):
    global NUM_BLOCK_INPUT, NUM_BLOCK, RANDOM_PICK, NUM_MUTATION, FP_RATIO, INPUT_DIR, OUTPUT_DIR

    NUM_BLOCK_INPUT = kwargs["NUM_BLOCK_INPUT"]
    NUM_BLOCK = kwargs["NUM_BLOCK"]
    RANDOM_PICK = kwargs["RANDOM_PICK"]
    NUM_MUTATION = RANDOM_PICK
    FP_RATIO = kwargs["FP_RATIO"]
    INPUT_DIR = "/data/project/Alzheimer/EM_cluster/pilot/04.EM_input/"
    OUTPUT_DIR = "./output/"

    membership_answer_old = membership_answer
    # 숫자로 바꿔주고 앞에서부터 numbering
    membership_answer_new = Sorting(membership_answer_old)
    membership_answer_max, score, sample_dict_rev = Scoring(  membership_answer_old, membership_answer_new, maxmaxmax_membership)

    return membership_answer_max, score, sample_dict_rev          # 100점 만점 score 반환
