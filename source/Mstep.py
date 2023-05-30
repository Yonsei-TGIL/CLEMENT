def main(df, np_vaf, np_BQ, step, option, **kwargs):  # 새로운 MIXTURE 정하는 과정
    import math
    import isparent, EMhard
    import visualizationsinglesoft,  visualizationeachstep
    import numpy as np

    NUM_BLOCK, NUM_CLONE, NUM_MUTATIN = kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"], kwargs["NUM_MUTATION"]
    NUM_MUTATION = kwargs["RANDOM_PICK"]

    kwargs["OPTION"] = option

    if option in ["Hard", "hard"]:
        ############################### HARD CLUSTERING ##############################
        for j in range(NUM_CLONE):
            ind_list = np.where(step.membership == j)[0]   # membership == j 인 index를 구하기
            for i in range(NUM_BLOCK):
                sum_depth, sum_alt = 0, 0
                for ind in ind_list:       # depth, alt를 다 더하기
                    if df[ind][i]["alt"] != 0:
                        sum_depth = sum_depth + df[ind][i]["depth"]
                        sum_alt = sum_alt + df[ind][i]["alt"]
                step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0 # j번째 clone만을 생각한 이상적인 분율을 일단 assign

        
        
        previous_fp_index = step.fp_index
        step.makeone_index, p_list, step.fp_index = isparent.makeone(df, np_vaf, np_BQ, step, **kwargs)


        if step.fp_index != -1:  # FP clone이 있다면
            step.includefp = True
            step.fp_member_index = list(np.where(step.membership == step.fp_index)[0])
        else:   # FP clone이 없다면
            step.includefp = False
            step.fp_member_index = []



        # if NUM_CLONE == 1:    # 이란 무조건 1로 꽂고 본다. 근데 이 과정이 꼭 필요한가?
        #     step.mixture = np.array([[1.0]] * kwargs["NUM_BLOCK"])

        if step.makeone_index == []:
            step.likelihood = -9999999
            step.makeone_prenormalization = False
            print ("\t\t\tMstep | 합쳐서 1을 못 만든다", end = "\t"  )
            print( ",".join( str(row) for row in step.mixture ))

        if kwargs["STEP"] <= 4:
            if step.makeone_index != []:
                step.makeone_prenormalization, sum_mixture =  EMhard.checkall (step, **kwargs)  # prenormalization mixture를 가지고 0.9 ~ 1.1 사이에 확실히 들어오는지를 판단
                print ("\t\t\tMstep | step.makeone_prenormalization = {}".format(step.makeone_prenormalization), end = "\t") 
                print(" ".join(str(row) for row in sum_mixture ))

        
                for i in range(NUM_BLOCK):
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:                # boundary clone들만 붙잡고 normalization
                            sum = sum + step.mixture[i][j]
                    step.mixture[i] = np.round( step.mixture[i] / sum, 4) if sum != 0 else 0   # 만약 sum이 0이면 분모에 0이 들어갈 수 있으니까...
                #print ("\t\tafter-normalization : {}".format (step.mixture))
        
    if (kwargs["NUM_BLOCK"] == 1):
        visualizationeachstep.drawfigure_1d_hard(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard).pdf", **kwargs)
    if (kwargs["NUM_BLOCK"] >= 2):
        visualizationeachstep.drawfigure_2d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard).pdf", **kwargs)
    ###############################################################################

    ################################ SOFT CLUSTERING ##############################

    if option in ["Soft", "soft"]:
        #print ("\t\ta. Mixture (before soft clustering) : {}". format(list(step.mixture)))

        makeone_index_i = []
        for k in range(NUM_MUTATION):
            if step.membership[k] in step.makeone_index:
                makeone_index_i.append(k)

        for j in range(NUM_CLONE):
            # child가 아닌 애들 (parent, fp)은 mixutre도 그냥 hard처럼 mixture를 계산해준다
            if j not in step.makeone_index:
                for i in range(NUM_BLOCK):
                    sum_depth, sum_alt = 0, 0
                    # depth, alt를 다 더하기
                    for ind in np.where(step.membership == j)[0]:
                        if df[ind][i]["alt"] != 0:
                            sum_depth = sum_depth + df[ind][i]["depth"]
                            sum_alt = sum_alt + df[ind][i]["alt"]
                    # j 번째 clone만을 생각한 이상적인 분율을 일단 assign
                    step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0

            elif j in step.makeone_index:  # child clone들만 soft clustering 한다
                for i in range(NUM_BLOCK):   # 각 block에 대해 weighted mean을 계산
                    vaf, weight = np.zeros(NUM_MUTATION, dtype="float"), np.zeros(NUM_MUTATION, dtype="float")
                    for k in range(NUM_MUTATION):
                        vaf[k] = int(df[k][i]["alt"]) / int(df[k][i]["depth"])
                        # weight는 block과는 상관없다
                        # if (len(step.membership_p[k]) == NUM_CLONE + 1) & (step.membership_p[k][-1] == 0):   # FP 처리한 mutation의 경우 완전 뺴기 위해서
                        #     weight[k] = 0

                        if step.membership[k] in step.makeone_index:
                            weight[k] = math.pow(10, step.membership_p[k][j])
                        #print ("{} : (clone {})   weight = {},  vaf = {}".format(k, step.membership[k], weight[k], vaf[k]))

                    step.mixture[i][j] = round(np.average(vaf[makeone_index_i], weights=weight[makeone_index_i]), 4) * 2

        #print ("\t\tb. Mixture (after soft clustering) : {}". format(list(step.mixture)))

        # Block당 Mixture의 합이 1이 되도록 재조정
        if kwargs["adjustment"] in ["True", "true", True]:
            for i in range(NUM_BLOCK):
                sum = np.sum(step.mixture[i])
                # 만약 sum이 0이면 분모에 0이 들어갈 수 있으니까...
                step.mixture[i] = np.round(step.mixture[i] / sum, 4) if sum != 0 else 0

        elif kwargs["adjustment"] in ["Half", "half"]:
            # step.makeone_index, p_list, step.fp_index = isparent.makeone (step, **kwargs)

            if step.fp_index != -1:  # FP clone이 있다면
                step.includefp = True
                step.fp_member_index = list(np.where(step.membership == step.fp_index)[0])
            else:   # FP clone이 없다면
                step.includefp = False
                step.fp_member_index = []

            if NUM_CLONE == 1:
                step.mixture = np.array([[1.0]] * kwargs["NUM_BLOCK"])

            if step.makeone_index == []:
                step.likelihood = -9999999
                step.makeone_prenormalization = False
                # print ("여러 조건 (parent, child) 때문에 합쳐서 1을 못 만든다\t", list(mixture))

        if kwargs["STEP"] <= 5:
            if step.makeone_index != []:
                for i in range(NUM_BLOCK):
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:                # 1과 가장 가까운 clone들만 붙잡아줌
                            sum = sum + step.mixture[i][j]
                    # 만약 sum이 0이면 분모에 0이 들어갈 수 있으니까...

                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:                # 1과 가장 가까운 clone들만 붙잡아줌
                            step.mixture[i][j] = np.round(step.mixture[i][j] / sum, 4) if sum != 0 else 0

        #print ("\t\tc. Mixture (after normalization) : {}". format(list(step.mixture)))

        
        step.membership_p_normalize = np.zeros((NUM_MUTATION, step.membership_p.shape[1]), dtype="float64")
        for k in range(NUM_MUTATION):
            if k in step.fp_member_index:
                step.membership_p_normalize[k] = np.zeros(step.membership_p_normalize.shape[1], dtype="float64")  # 1 (FP_index) 0 0 0 0 0 으로 만들어준다
                # 1 (FP_index) 0 0 0 0 0 으로 만들어준다
                step.membership_p_normalize[k][step.fp_index] = 1
            else:
                step.membership_p_normalize[k] = np.round(np.power(10, step.membership_p[k])/np.power(10, step.membership_p[k]).sum(axis=0, keepdims=1), 2)   # 로그를 취해으나 다시 지수를 취해준다
                if step.fp_index != -1:     # fp가 있을 때에만...  fp가 없는데도 -1 column이 0으로 되버리는 참사
                    # 0  (FP_index) 0.2 0.7 0.1 0 0 으로 만들어준다
                    step.membership_p_normalize[k][step.fp_index] = 0

        if (kwargs["NUM_BLOCK"] == 1):
            visualizationeachstep.drawfigure_1d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft).pdf", **kwargs)
        if (kwargs["NUM_BLOCK"] == 2):
            visualizationeachstep.drawfigure_2d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft).pdf", **kwargs)

    #############################################################################

    return step
