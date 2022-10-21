import scipy.stats
import numpy as np
import comb, math

def greaterall (a, b, boole):
    if boole == "<":
        for i in range(len(a)):
            if a[i] > b[i]:
                return False
        return True
    if boole == ">":
        for i in range(len(a)):
            if a[i] < b[i]:
                return False
        return True

def checkall (sum_mixture):
    for i in range ( len(sum_mixture) ):
        if ( sum_mixture[i] < 0.77 ) | ( sum_mixture[i] > 1.25 ):       # If at least one sample dissatisfies this rule (sum of mixture = 1), than reject  (return False)
            return False
    return True


def makeone (step, **kwargs):
    membership = step.membership
    mixture = step.mixture

    global subset_list_acc, subset_mixture_acc, sum_mixture_acc

    
    if step.fp_index != -1:       # Let's not consider putative FP as a child clone candidate
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(list(set(membership) - set( [step.fp_index ])), mixture)   # All possible combinations
    elif step.fp_index == -1:
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(list(set(membership))[: ], mixture)   # All possible combinations


    if kwargs["NUM_CLONE_NOMINAL"] == 1:   #  It is weird that comb.comball does not work in NUM_CLONE == 1
        return [0], [[1, 0]], -1          # makeone_index , p_list, fp_index

    
    p_max = float("-inf")
    p_list, j2_list = [], []
    fp_possible = []
    fp_possible_clone = np.zeros (len(subset_mixture_acc), dtype = "int")

    for j2 in range(len(subset_mixture_acc)):      #All possible combinations
        subset_list, subset_mixture, sum_mixture = subset_list_acc[j2], subset_mixture_acc[j2], sum_mixture_acc[j2]

        if checkall(sum_mixture) == False:         # 한 If at least one sample dissatisfies this rule (sum of mixture = 1), than reject  (return False)
            continue

        p = 0
        for i in range (kwargs["NUM_BLOCK"]):
            depth = 1000
            a = int(sum_mixture[i] * 1000 / 2) 
            b = depth - a

            target_a = 500
            try:
                p = p + math.log10(scipy.stats.betabinom.pmf(target_a, depth, a + 1, b+1))
            except:
                p = p - 400
        if p > -400:
            ISSMALLER_cnt, SMALLER_ind = 0, []
            ISPARENT = False

            for j3 in (set(range(0, mixture.shape[1]  )) - set(subset_list)) :    # Remaining clone
                if (j3 == step.fp_index) | ( len (set (membership)) <=  j3 )  :
                    continue
                ISLARGER_cnt = 0
                possible_parent_cnt = 0
                for j4 in subset_list :    #선정된 boundary 후보          나머지 clone중에 child clone보다 작으면 안됨.  그러나 딱 1개만 있고 FP clone이면 용서해준다
                    if greaterall ( mixture[:,j3] , mixture[:,j4], "<" ) == True:
                        ISLARGER_cnt = ISLARGER_cnt + 1
                    if greaterall ( mixture[:,j3] , mixture[:,j4], ">" ) == True:   # Parent clone이 있는지 알아보기
                        #print (mixture[:,j3],  mixture[:,j4])
                        possible_parent_cnt = possible_parent_cnt + 1
                
                if ISLARGER_cnt >= 1:       # 적어도 1개의 boundry 후보 (child 후보)들 보다 작다는게 말이 안됨.  다만 FP clone이 있을 경우에는 다르다
                    ISSMALLER_cnt = ISSMALLER_cnt + 1
                    SMALLER_ind.append (j3)

                if possible_parent_cnt == len (subset_list):
                    #print ("{} is larger than all the child clone : {}".format (j3, mixture))
                    ISPARENT = True
                    break
            
            if ( ( kwargs ["NUM_BLOCK"] == 1 ) & ( ISPARENT == True ) ):   # 1D: Ignore the possibility that parent clone exist
                continue
                
                
            if (step.includeoutlier == False) &  (ISSMALLER_cnt == 1):       # 그동안 FP clone 없었는데 유일하게 안쪽에서 발견할 경우
                #print ("0 :  {} 번째 clone은 {}보다 안쪽에 있는 유일한 clone이라 살리는걸 고려.  개수는 {}개. ".format(SMALLER_ind[0], subset_list,   int(np.unique(membership, return_counts=True)[1] [SMALLER_ind[0]] )))
                check = 0

                for j3 in set(range(0, mixture.shape[1] )) - set(subset_list) - set([SMALLER_ind[0]]):   # 나머지 clone (outlier, FP도 빼고) 을 돈다
                    if j3 == step.fp_index:
                        continue

                    tt = []
                    for j4 in subset_list:   # boundary clone (putative child clone)을 돈다
                        if greaterall ( mixture[:,j3] , mixture[:,j4], ">" ) == True:
                            tt.append (mixture[:,j4])

                    if len(tt) < 2:        # 나머지 clone은 2개 이상의 child 조합으로 만들어지는 parent여야 한다. 그게 만족 안하면 이 조합이 오류임을 알 수 있다
                        #print ("{} ( {} ) 는 2개 이상의 child clone의 합으로 만들어지지지 않아서 이 조합 (FP 포함, {}) 은 아예 기각".format(j3 , mixture[:, j3], subset_mixture))
                        check = 1
                        break

                if check == 0:
                    if SMALLER_ind[0] >= np.max(membership) + 1:
                        print ("mixture : {}".format(list(mixture)))
                        print ("SMALLER_ind : {}".format(SMALLER_ind))
                        print ("membership : {}".format(np.unique(membership, return_counts=True)[1]))

                    try:
                        if ( int(np.unique(membership, return_counts=True)[1] [SMALLER_ind[0]] )  > kwargs["MIN_CLUSTER_SIZE"]  ):
                            fp_possible.append (j2)
                            fp_possible_clone[j2] =  SMALLER_ind[0]
                            p_list.append ( [ p, j2] )
                    except:
                        print ("Something error happend in isparent.py")
                        print ("\tFP Clone number : {}".format(SMALLER_ind[0]))
                        print ( int(np.unique(membership, return_counts=True) ) )



            if ISSMALLER_cnt == 0:      # No FP clone included & Pure child clone
                check = 0
                for j3 in set(range(0, mixture.shape[1] )) - set(subset_list):   # Iterating remaining clone
                    if j3 == step.fp_index:
                        continue
                    tt = []
                    for j4 in subset_list:   # Iterating boundary clone (putative child clone)
                        if greaterall ( mixture[:,j3] , mixture[:,j4], ">" ) == True:
                            tt.append (mixture[:,j4])

                    if len(tt) < 2:        # 나머지 clone은 2개 이상의 child 조합으로 만들어지는 parent여야 한다. 그게 만족 안하면 이 조합이 오류임을 알 수 있다
                        #print ("{} ( {} ) 는 2개 이상의 child clone의 합으로 만들어지지지 않아서 이 조합은 (FP 미포함, {})은 아예 기각".format(j3 , mixture[:, j3], subset_mixture ))
                        check = 1
                        break
                if check == 0:
                    p_list.append ( [ p, j2] )        # Put boundary clone (child clone) into p_list

            #print ( subset_list, "ISLARGER_cnt = {}\tISSMALLER_cnt = {}\tmixture = {}".format (ISLARGER_cnt, ISSMALLER_cnt, mixture) )


    if p_list == []:
        return [], [], -1

    p_list = np.array(p_list).transpose()
    p_list = p_list[ :, p_list[0].argsort()]
    p_list = np.flip(p_list , 1)

    if kwargs["VERBOSE"] >= 2:
        if step.fp_index != -1:
            for i in range(len(p_list[1])):
                print ("step={}\tclone = {}, sum = {} : p = {}".format (kwargs["STEP_TOTAL"], subset_list_acc[ int(p_list[1][i]) ] ,  sum_mixture_acc[ int(p_list[1][i]) ],  p_list[0][i]))
            print ("")

    best_j2 = int( p_list[1,0] )

    if best_j2 in fp_possible:
        if kwargs["VERBOSE"] >= 2:
            print ("\t\t{}th clone is the only clone inside {}, so survived. At the same time, this combination ({}) became the highest probable makeone".format(fp_possible_clone[best_j2], subset_list_acc[best_j2], subset_list_acc[best_j2]))
        return subset_list_acc[best_j2], p_list, fp_possible_clone[best_j2]

        # try:
        #     print ("NUM_CLONE : {}\t{}-{}\tmakeone_index\t{}\tsum\t{}".format(mixture.shape[1], kwargs["TRIAL"], kwargs["STEP"], subset_list_acc[best_j2], sum_mixture_acc[best_j2] ) )
        # except:
        #     print ("NUM_CLONE : {}\tmakeone_index\t{}\tsum\t{}".format(mixture.shape[1], subset_list_acc[best_j2], sum_mixture_acc[best_j2] ) )

    return subset_list_acc[best_j2], p_list, step.fp_index

