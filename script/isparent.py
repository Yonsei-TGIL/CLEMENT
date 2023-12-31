import scipy.stats
import numpy as np
import itertools, miscellaneous
import comb
import math
import Estep


def greaterall(a, b, boole, **kwargs):
    if boole == "<":
        for i in range(len(a)):
            if a[i] > b[i]:
                return False
        return True
    if boole == ">":
        for i in range(len(a)):
            if a[i] < b[i]:
                if kwargs ["NUM_BLOCK"] == 3:
                    if a[i] != 0:           # 231011 ( for 3D )
                        return False
                else:
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
        #print ("child_list = {}".format(child_list))
        for i in range (kwargs["NUM_BLOCK"]):
            if kwargs ["NUM_BLOCK"] == 3:
                if step.mixture [i][j3] == 0:
                    continue

            depth = 1000
            a = int(sum_mixture_child[i] * 1000 / 2) 
            b = depth - a
            target_a = int (step.mixture[i][j3] * 1000/ 2)
            try:
                p = p + math.log10(scipy.stats.betabinom.pmf(target_a, depth, a + 1, b+1))
        
                #print ( "\tNUM_BLOCK = {}\ttarget_a = {}\ta = {}\tdepth = {}\tp = {}".format (i, target_a, a, depth,  p))
            except:
                p = p - 400

        if p > Phy.p:
            Phy.p = round (p, 2)
            Phy.child_list = child_list
            Phy.sum_mixture_child = sum_mixture_child
            #print ("\t\t\t\t\tp = {}, child_list = {}, sum_mixture_child = {}".format ( round(p, 2), child_list, sum_mixture_child))

    if kwargs["VERBOSE"] >= 3:
        print ("\t\t\t\t\t\t→ j3 = {}, Phy.child_list = {}, Phy.sum_mixture_child = {}, Phy.j3_mixture[:,j3] = {}, Phy.p = {}".format ( j3, Phy.child_list, Phy.sum_mixture_child, step.mixture[:,j3], round (Phy.p, 2) ) )
    
                
    Threshold = -3 - NUM_BLOCK
    if kwargs ["STEP"] < kwargs["COMPULSORY_NORMALIZATION"] :
        Threshold = Threshold - (4 - kwargs["STEP"]) / 2
                
    if Phy.p < ( Threshold ):           # 이 기준이 애매하긴 한데.... 나중에 distribution을 보고 판단할 필요가 있다
        Phy.appropriate = False
    
    return Phy






def makeone(input_containpos, df, np_vaf,  np_BQ, step, **kwargs):
    import miscellaneous, copy

    membership = step.membership
    mixture = step.mixture
    NUM_BLOCK = step.mixture.shape[0]
    NUM_CLONE = step.mixture.shape[1]
    kwargs["CHECKALLFROM"] = "isparent.py"
    

    global subset_list_acc, subset_mixture_acc, sum_mixture_acc

    step.makeone_experience = False

    if step.includefp == True:       # Outlier가 있으면 그건 빼자
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(   list(set(membership) - set([step.fp_index])), mixture)   # 모든 덧셈 조합을 구하기
    elif step.includefp == False:
        subset_list_acc, subset_mixture_acc, sum_mixture_acc = comb.comball(  list(set(membership))[:], mixture)   # 모든 덧셈 조합을 구하기

    if len ( set ( step.membership) - set ( [kwargs["NUM_CLONE_ITER"]] )) == 1:   # 실제 clone이 1개일 때에는 comb.comball이 작동하지 않는다
        if kwargs["MAKEONE_STRICT"] in  [2, 3] :    # BioData에서는 극단적으로 lenient하게 잡아줘야 하니까
            checkallt_pass, sum_mixture, p, kwargs  = miscellaneous.checkall ( step, "lenient", np_vaf, **kwargs)
        else:
            checkallt_pass, sum_mixture, p, kwargs = miscellaneous.checkallttest ( step, "lenient", np_vaf, **kwargs ) 

        if checkallt_pass == False :
            p_list = []
            return p_list, step, "lenient", kwargs      # 그 어떤 makeone도 만들지 못해서 그냥 넘긴다
        else:  # 성공한 경우
            step.makeone_index = [0]
            p_list = [[1, 0]]
            return p_list, step, "lenient", kwargs          # makeone_index , p_list, fp_index

    p_max = float("-inf")
    p_list, j2_list = [], []
    
        
    step_subset = copy.deepcopy ( step )
    if step.hard_or_soft == "hard":
        condition = "lenient" if ( kwargs ["STEP"] < kwargs["COMPULSORY_NORMALIZATION"]  ) else "strict"      # hard는 0, 1, 2를 봐줌
    else:
        condition = "lenient" if ( kwargs ["STEP"] <= kwargs["SOFT_LENIENT"]  ) else "strict"   # Soft는 1, 2, 3을 봐줌
        #condition = "lenient"    # Soft는 그냥 다 lenient로 봐줌

    # 여러 조합을 돈다 
    for j2 in range(len(subset_mixture_acc)):        # j2 : 번호, subset_list = [1, 3, 4]
        subset_list, subset_mixture, sum_mixture = subset_list_acc[ j2 ], subset_mixture_acc[ j2 ], sum_mixture_acc[ j2 ]
        step_subset.makeone_index = subset_list

        PhyAcc = PhylogenyObjectAcc()

        if condition == "lenient":
            checkallt_pass, sum_mixture, p, kwargs = miscellaneous.checkall ( step_subset, condition, np_vaf, **kwargs)
        elif condition == "strict":
            if kwargs["MAKEONE_STRICT"] in  [2, 3] : #  SimData, CellData Low depth, BioData
                checkallt_pass, sum_mixture, p, kwargs = miscellaneous.checkall ( step_subset, condition, np_vaf, **kwargs)
            else:
                checkallt_pass, sum_mixture, p, kwargs = miscellaneous.checkallttest ( step_subset, condition, np_vaf, **kwargs)
            
        if checkallt_pass == False:   
            if kwargs["VERBOSE"] >= 4:
                print ("\t\t\t\t(isparent.py)  makeone clone의 sum of mixture의 이상으로 기각 ({})".format( np.array(sum_mixture).flatten() ))
            continue

        step.makeone_experience = True

        # Beta binomial을 바탕으로 p를 계산
        if p == 0:
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
            
            possible_parent_index =  list ( set(range(0, mixture.shape[1])) - set(subset_list) - set ( [step.fp_index]) )
            for j3 in possible_parent_index:
                if (len(set(membership)) <= j3):
                    continue
                j3_lt_j4_cnt = 0
                j3gtj4_cnt = 0
                for j4 in subset_list:  # 선정된 makeone 후보          나머지 clone중에 child clone보다 작으면 안됨.  그러나 딱 1개만 있고 그게 FP clone이면 용서해준다
                    if greaterall(mixture[:, j3], mixture[:, j4], "<", **kwargs) == True:
                        j3_lt_j4_cnt = j3_lt_j4_cnt + 1
                    if greaterall(mixture[:, j3], mixture[:, j4], ">", **kwargs) == True:                        # Parent clone이 있는지 알아보기
                        j3gtj4_cnt = j3gtj4_cnt + 1

                #print ("j3 = {}\tj3_lt_j4_cnt = {}\tj3gtj4_cnt = {}".format (mixture[:, j3], j3_lt_j4_cnt, j3gtj4_cnt))

                # 명색이 parent 후보인데, 1개라도 makeone 후보 (child 후보)들 보다 작다는게 말이 안됨.  다만 FP clone이 있을 경우에는 다르다
                if j3_lt_j4_cnt >= 1:     # 1개라도 makeone 후보들보다 모든 차원에서 작은 경우 -> 이런 j3 (SAMLLER_ind)가 딱 1개이고 그게 FP clone이면 용서해준다
                    ISSMALLER_cnt = ISSMALLER_cnt + 1
                    SMALLER_ind.append(j3)

                if j3gtj4_cnt >= 2:       # j3가 2개 이상의  j4 (child clone) 보다 클 때 -> parent의 가능성이 있다
                    if (kwargs["VERBOSE"] >= 3) & ( len(PhyAcc.p_record) == 0):
                        print ("\t\t\t\t\tj2 = {}, subset_list = {}".format(j2, subset_list))
                    Phy = Appropriate_Phylogeny_Verification (PhyAcc, subset_list, j2, j3, step, **kwargs)
                    if Phy.appropriate == False:
                        if kwargs["VERBOSE"] >= 3:
                            print ("\t\t\t\t\t\t\t→ Phylogeny 확률이 너무 안좋다. 따라서 이 j2는 withdrawl".format(j3))
                        break
                    else:
                        PhyAcc.acc (Phy)

           

            
            if (Phy.appropriate == False): # Phylogeny가 영 믿음이 안가서 이 j2는 기각하고 다음으로 넘어간다
                continue

            if len( PhyAcc.j3_record )  > kwargs["MAXIMUM_NUM_PARENT"]:      # len( PhyAcc.j3_record ) = parent의 숫자
                if kwargs["VERBOSE"] >= 1:
                    print ( "\t\t\t\t\t\t→ possible_parent_num ( {} ) > MAXIMUM_NUM_PARENT ( {})".format ( len( PhyAcc.j3_record ), kwargs["MAXIMUM_NUM_PARENT"]))
                continue
            
            # 1D일 경우에는 parent clone의 존재를 인정하지 말자
            # if ((kwargs["NUM_BLOCK"] == 1) & (len( PhyAcc.j3_record )  > 0 )  ):
            #     if kwargs["VERBOSE"] >= 1:
            #         print ("\t\t\t1D이므로 parent clone을 잡는걸 인정하지 않는다    ->  parent : {}\t child : {}".format (j3, subset_list ))
            #     continue
            
            if kwargs["VERBOSE"] >= 4:
                print ( "\t\t\t\t\t\tisparent.py : subset_list = {}\tpossible_parent_index = {}\tISSMALLER_CNT = {}\tstep.includefp = {}".format (subset_list, possible_parent_index, ISSMALLER_cnt, step.includefp) )





            # 1. 그동안 FP clone 있었는데, 안쪽에서 또 발견될 경우 → 살려줄 가치가 없다
            # 2. (이젠 X) 그동안 FP clone 없었는데, 유일하게 안쪽에서 발견될 경우  → 별다른 문제 없으면 살려주고, 걔를 fp로 지명해줌
            # 3. 그동안 FP clone 없었는데, 나머지들로 makeone 잘 할 경우 → 그래도 한 번 정도는 살려줄 수 있다
            # 4. 그동안 FP clone 있었는데, 나머지들로 makeone 잘 할 경우 → 당연히 살려줄 수 있음

            if (step.includefp == True) & (ISSMALLER_cnt > 0):      # 1번 상황
                if kwargs["VERBOSE"] >= 4:
                    print ("1번 상황이라 제낀다")
                continue

            if ISSMALLER_cnt == 0:      # 3번, 4번 상황
                check = 0
                # Parent가 정당한지 검사
                for j3 in possible_parent_index: # j3: 나머지 clone (fp cluster를 빼고) 을 돈다
                    tt = []
                    for j4 in subset_list:   # j4: makeone clone (putative child clone)을 돈다
                        if greaterall(mixture[:, j3], mixture[:, j4], ">", **kwargs ) == True:
                            tt.append(mixture[:, j4])

                    if len(tt) < 2:         # 나머지 clone은 2개 이상의 child 조합으로 만들어지는 parent여야 한다. 그게 만족 안하면 이 조합이 오류임을 알 수 있다
                        if kwargs["VERBOSE"] >= 4:
                            print ("\t\t\t\t\t\t→ {} 가 2개 이상의 child clone의 합으로 만들어지지지 않아서 이 j2는 아예 기각".format(j3  ))
                        check = 1
                        break


                if check == 0:   # Parent 검사에서 통과하면
                    if step.includefp == True:   # 기존 fp가 있었다면 4번상황
                        p_list.append([p, j2, 4 ])   # makeone_p, j2 (subset_list), 4번상황, p_Phylogeny
                    else:    # 기존 fp가 없었다면 3번상황
                        p_list.append([p, j2, 3 ])    # makeone_p, j2 (subset_list), 3번상황, p_Phylogeny

                    if kwargs["VERBOSE"] >= 3:
                        if len( PhyAcc.j3_record ) != 0:
                            print ("\t\t\t\t→ parent : {}, PhyAcc.p_record = {}, sum_mixture = {}, p = {}".format( PhyAcc.j3_record, PhyAcc.p_record, sum_mixture , round (p, 2) ))
                        else:
                            print ("\t\t\t\t→ parent는 없음, p = {}".format(round (p, 2) ))


    if p_list == []: # 실패한 경우
        return p_list, step, condition, kwargs     # 그 어떤 makeone도 만들지 못해서 그냥 넘긴다

    else:  # 성공한 경우
        p_list = np.array(p_list).transpose()
        p_list = p_list[:, p_list[0].argsort()]  # p_list[0]  (확률)을 기준으로  p_list (나머지 row) 를 다 sort 해버리자
        p_list = np.flip(p_list, 1)
        

        if kwargs["VERBOSE"] >= 3:
            for i in range (0, p_list.shape[1]):
                j2 = int(p_list[1, i])
                print ("\t\t\t∇ {}nd place : subset_list_acc [j2] = {}\tsum_mixture_acc [j2] = {}\t{}th cirumstance\tp = {}".format ( i + 1 , subset_list_acc [ j2  ], sum_mixture_acc [ j2  ], int (p_list[2, i]) , round( p_list[0, i], 2)  ))

        best_j2 = int(p_list[1, 0])    # 1 : subset_list의 index  0 : 제일 잘한것 (0등)
        optimal, optimal_j2 = 0, best_j2

        
        if ( p_list[2, optimal] == 3):     # optimal : 0 or 1   ← 웬만하면 best인 0을 고르겠지만...
            if kwargs["VERBOSE"] >= 3:
                print ("\t\t\t3번상황 발생 (그동안 FP clone 없었는데, 나머지들로 makeone 잘 할 경우) → fp_index = {}".format ( step.fp_index) )
        elif ( p_list[2, optimal] == 4):
            if kwargs["VERBOSE"] >= 3:
                print ("\t\t\t4번상황 발생 (그동안 FP clone 있었는데, 나머지들로 makeone 잘 할 경우) → fp_index = {}".format (step.fp_index) )
        step.fp_involuntary = False
        step.makeone_index = subset_list_acc[optimal_j2]
        
        if condition == "lenient":
            step.checkall_lenient = True
        elif condition == "strict":
            step.checkall_strict= True
        
        return p_list, step, condition, kwargs
    


