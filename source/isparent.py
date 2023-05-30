import scipy.stats
import numpy as np
import itertools
import comb
import math
import Estep


def greaterall(a, b, boole):
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


def checkall(sum_mixture, **kwargs):      # 여는 좀  느슨하게 잡는다
    NUM_BLOCK = len(sum_mixture) 

    if kwargs["MAKEONE_STRICT"] == 1:
        makeone_standard = np.array ( [ [0.84, 1.15], [0.84, 1.15] ],dtype = float)
    elif kwargs["MAKEONE_STRICT"] == 2:
        makeone_standard = np.array ( [ [0.77, 1.25], [0.77, 1.25] ],dtype = float)
    
    if kwargs["NUM_BLOCK"] == 1:      # 1D
        if (sum_mixture[0] < makeone_standard[0][0]) | (sum_mixture[0] > makeone_standard[0][1]):
            return False
        else:
            return True
    else:
        for i in range( kwargs["NUM_BLOCK"] ):  # 한 block이라도 1에 걸맞지 않은 게 있으면 False를 return함
            if (sum_mixture[i] < makeone_standard[1][0]) | (sum_mixture[i] > makeone_standard[1][1]): 
                return False
        return True

############################################################################################################
class PhylogenyObject:
    def __init__(self):
        self.appropriate = True
        self.p = float ("-inf")
        self.child_list = []
        self.sum_mixture_child = []
        self.j3 = -1
        self.j3_mixture = []
        
class PhylogenyObjectAcc:
    def __init__(self):
        self.appropriate_record = []
        self.p_record = []
        self.child_list_record = []
        self.sum_mixture_child_record = []
        self.j3_record = []
        self.j3_mixture_record = []
    def acc(self, Phy):
        self.appropriate_record.append (Phy.appropriate)
        self.p_record.append (Phy.p)
        self.child_list_record.append (Phy.child_list)
        self.sum_mixture_child_record.append (Phy.sum_mixture_child)
        self.j3_record.append (Phy.j3)
        self.j3_mixture_record.append (Phy.j3_mixture)
############################################################################################################

def Appropriate_Phylogeny_Verification (PhyAcc, subset_list, j2, j3, step, **kwargs):
    import itertools

    Phy = PhylogenyObject ()
    Phy.j3_mixture =  step.mixture[:,j3]
    Phy.j3 = j3

    NUM_BLOCK, NUM_CLONE = step.mixture.shape[0], step.mixture.shape[1]

    combi =  itertools.chain(*map(lambda x: itertools.combinations(subset_list, x), range(2, len(subset_list)+1)))   # 2 ~ n개씩 짝지은 조합
    

    for child_list in combi:
        if child_list in PhyAcc.child_list_record:   # 다른 parent가 같은 구성으로 이루어졌다면, 피하자
            #print ("j3 = {}, child_list = {} → PhyAcc에 있다".format( j3, child_list))
            continue
        
        # 가능한 child clone끼리 합친 sum_mixture_child를 구하기
        sum_mixture_child = np.zeros (  NUM_BLOCK, dtype = "float")
        for i in range (NUM_BLOCK):
            for j4 in list ( child_list ) : 
                sum_mixture_child [i] += step.mixture [ i ][ j4 ]

        # Beta binomial을 바탕으로 계산하기
        p = 0
        for i in range (kwargs["NUM_BLOCK"]):
            depth = 1000
            a = int(sum_mixture_child[i] * 1000 / 2) 
            b = depth - a
            target_a = int (step.mixture[i][j3] * 1000/ 2)
            try:
                p = p + math.log10(scipy.stats.betabinom.pmf(target_a, depth, a + 1, b+1))
            except:
                p = p - 400

        if p > Phy.p:
            Phy.p = round (p, 2)
            Phy.child_list = child_list
            Phy.sum_mixture_child = sum_mixture_child
            #print ("\t\t\t\t\tp = {}, child_list = {}, sum_mixture_child = {}".format ( round(p, 2), child_list, sum_mixture_child))

    if kwargs["VERBOSE"] >= 3:
        print ("\t\t\t\t→ j3 = {}, Phy.child_list = {}, Phy.sum_mixture_child = {}, Phy.j3_mixture[:,j3] = {}, Phy.p = {}".format ( j3, Phy.child_list, Phy.sum_mixture_child, step.mixture[:,j3], round (Phy.p, 2) ) )
    
                
    Threshold = -3 - NUM_BLOCK
    if kwargs ["STEP"] <= 4:
        Threshold = Threshold - (4 - kwargs["STEP"]) / 2
                
    if Phy.p < ( Threshold ):           # 이 기준이 애매하긴 한데.... 나중에 distribution을 보고 판단할 필요가 있다
        Phy.appropriate = False
    
    return Phy






def makeone(df, np_vaf,  np_BQ, step, **kwargs):
    membership = step.membership
    mixture = step.mixture
    NUM_BLOCK = step.mixture.shape[0]
    NUM_CLONE = step.mixture.shape[1]

    global subset_list_acc, subset_mixture_acc, sum_mixture_acc

    if step.fp_index != -1:       # Outlier가 있으면 그건 빼자
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(   list(set(membership) - set([step.fp_index])), mixture)   # 모든 덧셈 조합을 구하기
    elif step.fp_index == -1:
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(  list(set(membership))[:], mixture)   # 모든 덧셈 조합을 구하기
    
    #subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(  list(set(membership))[:], mixture)   # 모든 덧셈 조합을 구하기

    if kwargs["NUM_CLONE_NOMINAL"] == 1:   # 이상하게 NUM_CLONE == 1일 때에는 comb.comball이 작동하지 않는다
        return [0], [[1, 0]], -1          # makeone_index , p_list, fp_index

    p_max = float("-inf")
    p_list, j2_list = [], []
    second_circumstance_j2 = []
    second_circumstance_fp_index = np.zeros(len(subset_mixture_acc), dtype="int")
    fourth_circumstance_j2 = []
    
        
    # 여러 조합을 돈다 
    for j2 in range(len(subset_mixture_acc)):
        subset_list, subset_mixture, sum_mixture = subset_list_acc[ j2 ], subset_mixture_acc[ j2 ], sum_mixture_acc[ j2 ]
        PhyAcc = PhylogenyObjectAcc()


        if checkall(sum_mixture, **kwargs) == False:         # 한 block이라도 1에 걸맞지 않은 게 있으면 False를 return함
            #print ("\t\t\t\t→ boundary clone의 sum of mixture의 이상으로 기각 ({})".format(sum_mixture))
            continue

        p = 0
        for i in range(kwargs["NUM_BLOCK"]):
            depth = 1000
            a = int(sum_mixture[i] * 1000 / 2)
            b = depth - a

            target_a = 500
            try:
                p = p + math.log10(scipy.stats.betabinom.pmf( target_a, depth, a + 1, b+1))
            except:
                p = p - 400
                
        
        if p > -400:
            ISSMALLER_cnt, SMALLER_ind = 0, []

            # 나머지 clone   (전체 clone -  child 후보 clone)
            Phy =  PhylogenyObject ()
            
            for j3 in (set(range(0, mixture.shape[1])) - set(subset_list) - set ( [step.fp_index]) ):
                if (len(set(membership)) <= j3):
                    continue
                j3_lt_j4_cnt = 0
                j3gtj4_cnt = 0
                for j4 in subset_list:  # 선정된 boundary 후보          나머지 clone중에 child clone보다 작으면 안됨.  그러나 딱 1개만 있고 그게 FP clone이면 용서해준다
                    if greaterall(mixture[:, j3], mixture[:, j4], "<") == True:
                        j3_lt_j4_cnt = j3_lt_j4_cnt + 1
                    if greaterall(mixture[:, j3], mixture[:, j4], ">") == True:                        # Parent clone이 있는지 알아보기
                        j3gtj4_cnt = j3gtj4_cnt + 1

                # 명색이 parent 후보인데, 1개라도 boundry 후보 (child 후보)들 보다 작다는게 말이 안됨.  다만 FP clone이 있을 경우에는 다르다
                if j3_lt_j4_cnt >= 1:     # 1개라도 boundary 후보들보다 작은 경우 -> 이런 j3 (SAMLLER_ind)가 딱 1개이고 그게 FP clone이면 용서해준다
                    ISSMALLER_cnt = ISSMALLER_cnt + 1
                    SMALLER_ind.append(j3)

                if j3gtj4_cnt >= 2:       # j3가 2개 이상의  j4 (child clone) 보다 클 때 -> parent의 가능성이 있다
                    if (kwargs["VERBOSE"] >= 3) & ( len(PhyAcc.p_record) == 0):
                        print ("\t\t\tj2 = {}, subset_list = {}".format(j2, subset_list))
                    Phy = Appropriate_Phylogeny_Verification (PhyAcc, subset_list, j2, j3, step, **kwargs)
                    if Phy.appropriate == False:
                        if kwargs["VERBOSE"] >= 4:
                            print ("\t\t\t\t\t→ Phylogeny 확률이 너무 안좋다. 따라서 이 j2는 withdrawl".format(j3))
                        break
                    else:
                        PhyAcc.acc (Phy)

           

            
            if (Phy.appropriate == False): # Phylogeny가 영 믿음이 안가서 이 j2는 기각하고 다음으로 넘어간다
                continue
        
            if len( PhyAcc.j3_record )  > kwargs["MAXIMUM_NUM_PARENT"]:      # len( PhyAcc.j3_record ) = parent의 숫자
                if kwargs["VERBOSE"] >= 3:
                    print ( "\t\t\t\t→ possible_parent_num ( {} ) > MAXIMUM_NUM_PARENT ( {})".format ( len( PhyAcc.j3_record ), kwargs["MAXIMUM_NUM_PARENT"]))
                continue
            
            # 1D일 경우에는 parent clone의 존재를 인정하지 말자
            if ((kwargs["NUM_BLOCK"] == 1) & (len( PhyAcc.j3_record )  > 0 )  ):
                #print ("1D + parent clone mode OFF 인데 parent clone을 잡았으므로 pass    ->  parent : {}\t child : {}".format (j3, subset_list ))
                continue






            # 1. 그동안 FP clone 있었는데, 안쪽에서 또 발견될 경우 → 살려줄 가치가 없다
            # 2. 그동안 FP clone 없었는데, 유일하게 안쪽에서 발견될 경우  → 별다른 문제 없으면 살려줄 수 있다
            # 3. 그동안 FP clone 없었는데, 나머지들로 makeone 잘 할 경우 → 별다른 문제 없으면 살려줄 수 있다
            # 4. 그동안 FP clone 있었는데, 나머지들로 makeone 잘 할 경우 → 별다른 문제 없으면 살려줄 수 있다

            # 2. 그동안 FP clone 없었는데 유일하게 안쪽에서 발견할 경우
            if (step.includefp == False) & (ISSMALLER_cnt == 1):   # 1개라도 boundary 후보들보다 작은 경우 -> 이런 j3가 딱 1개이고 그게 FP clone이면 용서해준다
                if kwargs["VERBOSE"] >= 4:
                    print ("\t\t\t\t{} 번째 clone은 {}보다 안쪽에 있는 유일한 clone이라 살리는걸 고려.  개수는 {}개. ".format(SMALLER_ind[0], subset_list,   int(np.unique(membership, return_counts=True)[1] [SMALLER_ind[0]] )))
                check = 0

                # 나머지 clone (putative parent clone) 을 돈다
                for j3 in set(range(0, mixture.shape[1])) - set(subset_list) - set([SMALLER_ind[0]]):
                    if j3 == step.fp_index:
                        continue
                    tt = []
                    # boundary clone (FP clone 제외한 putative child clone)을 돈다
                    for j4 in set (subset_list) - set([SMALLER_ind[0]]) :
                        if greaterall(mixture[:, j3], mixture[:, j4], ">") == True:
                            tt.append(mixture[:, j4])

                    # 나머지 clone은 2개 이상의 child 조합으로 만들어지는 parent여야 한다. 그게 만족 안하면 이 조합이 오류임을 알 수 있다
                    if len(tt) < 2:
                        if kwargs["VERBOSE"] >= 4:
                            print ("\t\t\t\t→ {} 가 2개 이상의 child clone의 합으로 만들어지지지 않아서 이 j2는 아예 기각".format(j3  ))
                        check = 1
                        break

                if check == 0:     # 2번 상황
                    if ( int(np.unique(membership, return_counts=True)[1][SMALLER_ind[0]]) > kwargs["MIN_CLUSTER_SIZE"]):
                        second_circumstance_j2.append(j2)
                        second_circumstance_fp_index [j2] = SMALLER_ind[0]
                        p_list.append([p, j2, 2 ])     # makeone_p, j2 (subset_list), 2번상황, p_Phylogeny
                        if kwargs["VERBOSE"] >= 3:
                            if len( PhyAcc.j3_record ) != 0:
                                print ("\t\t\t\t→ parent : {}, PhyAcc.p_record = {}, sum_mixture = {}, p = {}".format( PhyAcc.j3_record, PhyAcc.p_record, sum_mixture , round (p, 2) ))
                            else:
                                print ("\t\t\t\t→ parent는 없음, p = {}".format(round (p, 2) ))
                    else:
                        if kwargs["VERBOSE"] >= 3:
                            print ("\t\t\t\t→ {} <= MIN_CLUSTER_SIZE ({})라서 이 j2는 기각".format( int(np.unique(membership, return_counts=True)[1][SMALLER_ind[0]]) ,  kwargs["MIN_CLUSTER_SIZE"] ))
                        continue
                        

            if ISSMALLER_cnt == 0:      # 3번, 4번 상황
                check = 0
                for j3 in set(range(0, mixture.shape[1])) - set(subset_list): # j3: 나머지 clone (outlier를 빼고) 을 돈다
                    if j3 == step.fp_index:
                        continue
                    tt = []
                    for j4 in subset_list:   # j4: boundary clone (putative child clone)을 돈다
                        if greaterall(mixture[:, j3], mixture[:, j4], ">") == True:
                            tt.append(mixture[:, j4])

                    if len(tt) < 2:         # 나머지 clone은 2개 이상의 child 조합으로 만들어지는 parent여야 한다. 그게 만족 안하면 이 조합이 오류임을 알 수 있다
                        if kwargs["VERBOSE"] >= 3:
                            print ("\t\t\t\t→ {} 가 2개 이상의 child clone의 합으로 만들어지지지 않아서 이 j2는 아예 기각".format(j3  ))
                        check = 1
                        break
                    
                if check == 0:   # 별 문제 없으면
                    if step.fp_index != -1:   # 기존 fp가 있었다면 4번상황
                        #print  ("기존 fp_index ( {} ) 는 유지하되, 나머지들로 온전한 makeone을 만듬 (subset_list = {})".format (step.fp_index, subset_list) )
                        p_list.append([p, j2, 4 ])   # makeone_p, j2 (subset_list), 4번상황, p_Phylogeny
                    else:    # 기존 fp가 없었다면 3번상황
                        p_list.append([p, j2, 3 ])    # makeone_p, j2 (subset_list), 3번상황, p_Phylogeny

                    if kwargs["VERBOSE"] >= 3:
                        if len( PhyAcc.j3_record ) != 0:
                            print ("\t\t\t\t→ parent : {}, PhyAcc.p_record = {}, sum_mixture = {}, p = {}".format( PhyAcc.j3_record, PhyAcc.p_record, sum_mixture , round (p, 2) ))
                        else:
                            print ("\t\t\t\t→ parent는 없음, p = {}".format(round (p, 2) ))

    if p_list == []:
        return [], [], -1      # 그 어떤 makone도 만들지 못해서 그냥 넘긴다

    

    p_list = np.array(p_list).transpose()
    p_list = p_list[:, p_list[0].argsort()]  # p_list[0]  (확률)을 기준으로  p_list (나머지 row) 를 다 sort 해버리자
    p_list = np.flip(p_list, 1)
    

    if kwargs["VERBOSE"] >= 2:
        for i in range (0, p_list.shape[1]):
            j2 = int(p_list[1, i])
            print ("\t\t\t∇ {}위 : subset_list_acc [j2] = {}\tsum_mixture_acc [j2] = {}\t{}번상황\tp = {}".format ( i + 1 , subset_list_acc [ j2  ], sum_mixture_acc [ j2  ], int (p_list[2, i]) , round( p_list[0, i], 2)  ))

    best_j2 = int(p_list[1, 0])    # 1 : subset_list의 index  0 : 제일 잘한것 (0등)
    optimal, optimal_j2 = 0, best_j2
    if p_list.shape[1] >= 2: # 2위가 있어야 아래 상황도 고려
        if p_list [2, 0] in [2, 4]:   # 2: 몇번 상황인지  0 : 제일 잘한 것 (0등)  # 제일 잘한 subset에 FP가 포함되어 있냐?
            if p_list [2, 1] in [3]:    # 3번상황 : 그동안 FP clone 없었는데, 나머지들로 makeone 잘 할 경우
                print ("\t\t\t∇∇ FP를 새로 만들어준 2번 상황이 최적의 makeone이긴한데, FP 가 없는 3번상황이 바로 다음순위. 승부를 겨뤄보자")
                # for i in range (0, p_list.shape[1]):
                #     j2 = int(p_list[1, i])
                #     print ("\t\t\t∇ {}위 : j2 = {}, subset_list_acc [j2] = {}\tsum_mixture_acc [j2] = {}\t{}번상황\tp = {}".format ( i + 1 , j2, subset_list_acc [ j2  ], sum_mixture_acc [ j2  ], int (p_list[2, i]) , round( p_list[0, i], 2)  ))
                step.makeone_index = subset_list_acc[ best_j2 ]
                step.fp_index = second_circumstance_fp_index [best_j2]                        
                step_best = Estep.main (df, np_vaf, np_BQ, step, **kwargs )
                print ("\t\t\t\t∇ first j2 = {}\tstep.makeone_index = {}\tstep.fp_index = {}\tstep_best.likelihood = {}".format ( best_j2, step.makeone_index, step.fp_index, round(step_best.likelihood )  ))
                step_best_likelihood = step_best.likelihood
                
                secondbest_j2 = int(p_list[1, 1])
                step.makeone_index = subset_list_acc[ secondbest_j2 ]
                step.fp_index = -1                  
                step_secondbest = Estep.main (df, np_vaf, np_BQ, step, **kwargs )
                step_secondbest_likelihood = step_secondbest.likelihood
                print ("\t\t\t\t∇ secondbest j2 = {}\tstep.makeone_index = {}\tstep.fp_index = {}\tstep_secondbest.likelihood = {}".format ( secondbest_j2, step.makeone_index, step.fp_index, round (step_secondbest.likelihood) ))
            
    
                if step_secondbest_likelihood > step_best_likelihood: # 2위더라도 1위보다 더 좋은 
                    print ("\t\t\t\t→ secondbest를 선택한다")
                    optimal = 1      # 2위를 배정
                    optimal_j2 = secondbest_j2
            
    




    #print ("second_circumstance_j2 = {}".format(second_circumstance_j2))

    if optimal_j2 in second_circumstance_j2:      # 2번상황 (그동안 FP clone 없었는데, 유일하게 안쪽에서 발견될 경우)    
        new_fp_index = second_circumstance_fp_index[optimal_j2]
        if kwargs["VERBOSE"] >= 3:
            print ("\t\t\t2번상황 발생 (그동안 FP clone 없었는데, 안쪽에서 발견된걸 강제로 FP로 지정) → old_fp = {}\tnew_fp = {}".format (step.fp_index, new_fp_index))
        step.fp_involuntary = True
        return subset_list_acc[optimal_j2], p_list, new_fp_index
    
    else:
        if ( p_list[2, optimal] == 3):     # optimal : 0 or 1   ← 웬만하면 best인 0을 고르겠지만...
            if kwargs["VERBOSE"] >= 3:
                print ("\t\t\t3번상황 발생 (그동안 FP clone 없었는데, 나머지들로 makeone 잘 할 경우) → old_fp = {}\tnew_fp = {}".format (step.fp_index, step.fp_index) )
        elif ( p_list[2, optimal] == 4):
            if kwargs["VERBOSE"] >= 3:
                print ("\t\t\t4번상황 발생 (그동안 FP clone 있었는데, 나머지들로 makeone 잘 할 경우) → old_fp = {}\tnew_fp = {}".format (step.fp_index, step.fp_index) )
        step.fp_involuntary = False
        return subset_list_acc[optimal_j2], p_list, step.fp_index
