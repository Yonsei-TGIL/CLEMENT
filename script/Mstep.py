def lowvaf_and_peak ( step, i, j, sum_alt, sum_depth, ind_list, vaf_list, mode, **kwargs ):
    import numpy as np
    from scipy.stats import kde
    from scipy.stats import betabinom

    if len(ind_list) == 0:
        step.cluster_cdf_list[i].append ( 0 )
        step.mixture[i][j] = 0
        return step
    
    if step.mixture[i][j] == 0: # 2D에서 축에 있는 TN은 이런 계산에서 빼야지
        step.cluster_cdf_list[i].append ( 0 )
        return step
    

    cdf = betabinom.cdf( max( kwargs ["ALT_THRESHOLD"] * 2, int ( kwargs["MEAN_DEPTH"] * 0.05 * 2 ) ), int (sum_depth / len (ind_list) )  , int( (sum_depth / len (ind_list) ) * step.mixture[i][j]) + 1,  int((sum_depth / len (ind_list) ) * ( 1- step.mixture[i][j]) ) + 1 )
    #print ( "{}th block, {}th clone ({})의 alt <= {}일 cdf : {}".format (i, j, step.mixture[i][j], max ( kwargs ["ALT_THRESHOLD"], int ( kwargs["MEAN_DEPTH"] * 0.04  )  ) , round(cdf, 2) ) )
    step.cluster_cdf_list[i].append ( round(cdf, 2) )
    if ( (kwargs ["NUM_BLOCK"] == 1) & (cdf >= 0.01) ) | ( (kwargs ["NUM_BLOCK"] == 1) & (kwargs["MAKEONE_STRICT"] == 3) &  (cdf >= 0.005) )  | ( (kwargs ["NUM_BLOCK"] >= 2) & (cdf >= 0.05) ):
        step.lowvafcluster_index[i].append ( j )
        step.lowvafcluster_cdf_list[i].append ( round(cdf, 2) )
        if ( mode == "hard"):
            if ( len ( set(vaf_list) ) <= 1 ) | ( ( kwargs["MEAN_DEPTH"] >= 70 ) & ( kwargs["NUM_BLOCK"] >= 2) ) :                        
                step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 5) if sum_depth != 0 else 0   # Ideal centroid allocation
            else:     # peak, median 중 작은 것을 고른다
                x = np.linspace(0, 1.5, 601)
                try:
                    kde_np_vaf_new = kde.gaussian_kde( vaf_list )                
                    y = kde_np_vaf_new(x)     
                    x_max = x [np.argmax(y)]
                    step.mixture[i][j] = min ( x_max * 2,  round((sum_alt * 2) / sum_depth, 5) )   # peak, median 중 작은 것을 고른다
                    #print ( "j  = {}, i = {}\tcdf : {}\tcentroid(mean) : {}\tpeak : {}\tdecision : {}".format (j, i , round(cdf, 3), round((sum_alt * 2) / sum_depth, 5),   x_max * 2, step.mixture[i][j] ))
                except:     # 0만 가득차서 kde를 못 구할 떄
                    step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 5) if sum_depth != 0 else 0   # Ideal centroid allocation
            
    return step




def main(input_containpos, df, np_vaf, np_BQ, step, option, **kwargs):  
    import math, copy, math, subprocess
    import isparent, miscellaneous
    import visualizationeachstep
    import numpy as np
    from scipy.stats import kde

    NUM_BLOCK, NUM_CLONE, NUM_MUTATION = kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"], kwargs["NUM_MUTATION"]
    NUM_MUTATION = kwargs["RANDOM_PICK"]

    kwargs["OPTION"] = option

    if option in ["Hard", "hard"]:
        ############################### HARD CLUSTERING ##############################
        step.lowvafcluster_index =  [[] for _ in range(NUM_BLOCK)]
        step.lowvafcluster_cdf_list =  [[] for _ in range(NUM_BLOCK)]
        step.cluster_cdf_list =  [[] for _ in range(NUM_BLOCK)]
        for j in range(NUM_CLONE):
            if j == step.fp_index:  
                step.mixture[i][j] = 0   # FP cluster는 무조건 원점으로 박아 넣는다
            else:
                ind_list = np.where(step.membership == j)[0]   # Find the index  where membership == j                
                for i in range(NUM_BLOCK):
                    sum_depth, sum_alt = 0, 0
                    vaf_list = []
                    for ind in ind_list:       # Summing depth and alt
                        # if j in step.makeone_index:  
                        #     if df[ind][i]["alt"] == 0:  # TP cluster 에서는 FN을 평균치 계산에 넣지 않는다
                        #         continue
                        if j in step.tn_index:
                            zero_dimension = np.where( np.round ( step.mixture [:, j], 2 ) == 0 )[0]   # 0인 축 (sample)
                            if ( i in zero_dimension )  &  (df[ind][i]["alt"] != 0):    # TN : 0이어야 하는 축에서 0이 아니면 넘기자
                                continue
                        
                        #Most of the cases
                        sum_depth = sum_depth + df[ind][i]["depth"]
                        sum_alt = sum_alt + df[ind][i]["alt"]
                        if df[ind][i]["alt"] != 0:      # peak을 고를 때 vaf = 0은 제외
                            vaf_list.append ( df[ind][i]["alt"] / df[ind][i]["depth"])

                    #if kwargs ["MAKEONE_STRICT"] == 1:
                    step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 5) if sum_depth != 0 else 0   # Ideal centroid allocation

                    # lowvaf cluster인지 보고, 맞다면 step.mixture를 peak 기준으로 수정
                    if len ( [ind for ind in ind_list if df[ind][i]["alt"] == 0] ) == 0:       # 2D에서, alt = 0 인게 있다면 lowvaf cluster에 넣어줄 필요가 없음. 이미 mean을 취함으로써 편의를 봐주고 있기 때문
                        step = lowvaf_and_peak ( step, i, j, sum_alt, sum_depth, [ind for ind in ind_list if df[ind][i]["alt"] != 0], vaf_list, "hard", **kwargs )                    
    
        #print ( step.lowvafcluster_index )
        p_list, step, condition, kwargs = isparent.makeone(input_containpos, df, np_vaf, np_BQ, step, **kwargs)
        

        # checkall_strict, checkall_lenient 판단
        kwargs["CHECKALLFROM"] = "Mstep.py"
        if kwargs["MAKEONE_STRICT"] in  [2, 3]:       
            step.checkall_lenient, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "lenient", np_vaf, **kwargs) 
            step.checkall_strict, sum_mixture, p, kwargs     =  miscellaneous.checkall (step, "strict", np_vaf, **kwargs) 
        else:  # strict의 기준을 더 빡세게 잡아도 됨
            step.checkall_lenient, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "lenient", np_vaf, **kwargs) 
            step.checkall_strict, sum_mixture, p, kwargs =  miscellaneous.checkallttest (step, "strict", np_vaf, **kwargs) 


        if p_list == []:   # if failed
            #step.likelihood = step.likelihood_record[kwargs["STEP"]] = -9999999
            #step.checkall_lenient =  step.checkall_strict = False
            sum_mixture = np.sum (step.mixture[:, step.makeone_index], axis = 1)
            #step.makeone_index = []
            if kwargs["VERBOSE"] >= 1:
                if step.makeone_experience == False:
                    print ("\t\t\tMstep.py : condition = {}\tcheckall_strict = {}, checkall_lenient = {}\tUnable to make {} (sum = {}, standard = {}) ".format ( condition, step.checkall_strict, step.checkall_lenient, kwargs["MAKEWHAT"], " ".join(str( np.round(row, 2) ) for row in sum_mixture ), ",".join(str( row ) for row in kwargs["MAKEONE_STANDARD"] )  ) )
                else:
                    step.checkall_lenient =  step.checkall_strict = False
                    print ("\t\t\tMstep.py : condition = {}\tcheckall_strict = {}, checkall_lenient = {}\tUnable to make appropriate superclone-clone structure ({}) ".format ( condition, step.checkall_strict, step.checkall_lenient, ",".join(str(np.round(row, 2)) for row in step.mixture )  ) )

        else:    # most of the cases
            pre_mixture = copy.deepcopy (step.mixture)
            if kwargs["STEP"] < kwargs["COMPULSORY_NORMALIZATION"]:
                if kwargs["VERBOSE"] >= 1:
                    #print ("\t\t\tMstep.py : condition = lenient, checkall_strict = {}, checkall_lenient = {}\t(sum = {})".format(step.checkall_strict, step.checkall_lenient, " ".join(str( np.round(row, 2) ) for row in sum_mixture ) ) ) 
                    print ("\t\t\tMstep.py : condition = {}, checkall_strict = {}, checkall_lenient = {}\t(pre_mixture = {})".format(condition, step.checkall_strict, step.checkall_lenient, " ".join(str( np.round(row, 2) ) for row in pre_mixture ) ) ) 
                for i in range(NUM_BLOCK):     # Normalization 
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:   
                            sum = sum + step.mixture[i][j]
                    step.mixture[i] = np.round( step.mixture[i] / sum, 2) if sum != 0 else 0   # If sum = 0, let mixture = 0
                step.mixture = step.mixture * kwargs["MAKEWHAT"]    # normalization은 MAKEWHAT에 따라 달렸다 (30x 특히)

            else:
                if kwargs["VERBOSE"] >= 1:
                    #print ("\t\t\tMstep.py : condition = strict, checkall_strict = {}, checkall_lenient = {}\t(sum = {})".format(step.checkall_strict, step.checkall_lenient, " ".join(str( np.round(row, 2) ) for row in sum_mixture )  ) ) 
                    print ("\t\t\tMstep.py : condition = {}, checkall_strict = {}, checkall_lenient = {}\t(pre_mixture = {}, standard = {})".format(condition, step.checkall_strict, step.checkall_lenient, " ".join(str( np.round(row, 2) ) for row in pre_mixture ), ",".join(str( row ) for row in kwargs["MAKEONE_STANDARD"] )  ) ) 

                                
        if  kwargs["VISUALIZATION"] == True:
            if (kwargs["NUM_BLOCK"] == 1):
                visualizationeachstep.drawfigure_1d_hard(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture,**kwargs)
            elif (kwargs["NUM_BLOCK"] == 2):
                visualizationeachstep.drawfigure_2d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
            elif (kwargs["NUM_BLOCK"] >= 3):
                #visualizationeachstep.drawfigure_3d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
                visualizationeachstep.drawfigure_3d_SVD(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
            #subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "\(hard\)." + kwargs["IMAGE_FORMAT"] + " " + kwargs["COMBINED_OUTPUT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "\(hard\)." + kwargs["IMAGE_FORMAT"] ], shell = True)
    ###############################################################################

    ################################ SOFT CLUSTERING ##############################

    if option in ["Soft", "soft"]:
        #print ("\t\ta. Mixture (before soft clustering) : {}". format(list(step.mixture)))

        makeone_index_i = []
        for k in range(NUM_MUTATION):
            if step.membership[k] in step.makeone_index:
                makeone_index_i.append(k)

        step.lowvafcluster_index =  [[] for _ in range(NUM_BLOCK)]
        step.lowvafcluster_cdf_list =  [[] for _ in range(NUM_BLOCK)]
        step.cluster_cdf_list =  [[] for _ in range(NUM_BLOCK)]
        for j in range(NUM_CLONE):
            if j == step.fp_index:  
                step.mixture[i][j] = 0   # FP cluster는 무조건 원점으로 박아 넣는다
            else:
                if j not in step.makeone_index:                        
                    for i in range(NUM_BLOCK):
                        sum_depth, sum_alt, vaf_list, ind_list = 0, 0, [], np.where(step.membership == j)[0]
                        for ind in ind_list:
                            if df[ind][i]["alt"] != 0:  # TP cluster 라면 FN을 제거하기 위해 이래야 하고,  TN cluster라면 이럴 필요가 없다
                                # if (kwargs["SEX"] == "M") & ( bool(re.search(r'X|Y', input_containpos.iloc[ind]["pos"]))  == True  ) :
                                #     sum_depth = sum_depth + df[ind][i]["depth"] * 2
                                #     sum_alt = sum_alt + df[ind][i]["alt"]
                                # else: #Most of the cases
                                sum_depth += df[ind][i]["depth"]
                                sum_alt += df[ind][i]["alt"]
                                vaf_list.append ( df[ind][i]["alt"] / df[ind][i]["depth"])
                        step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0 

                        # lowvaf cluster인지 보고, 맞다면 step.mixture를 peak 기준으로 수정
                        step = lowvaf_and_peak ( step, i, j, sum_alt, sum_depth, [ind for ind in ind_list if df[ind][i]["alt"] != 0], vaf_list, "soft", **kwargs )                    

                else:     # 실질적인 soft clustering
                    for i in range(NUM_BLOCK):   # Calculate the weighted mean
                        sum_depth, sum_alt, vaf_list, ind_list = 0, 0, [], np.where(step.membership == j)[0]
                        for ind in ind_list:
                            if df[ind][i]["alt"] != 0:  # TP cluster 라면 FN을 제거하기 위해 이래야 하고,  TN cluster라면 이럴 필요가 없다
                                sum_depth += df[ind][i]["depth"]
                                sum_alt += df[ind][i]["alt"]
                                vaf_list.append ( df[ind][i]["alt"] / df[ind][i]["depth"])

                        vaf, weight = np.zeros(NUM_MUTATION, dtype="float"), np.zeros(NUM_MUTATION, dtype="float")
                        for k in range(NUM_MUTATION):
                            vaf [k] = int(df[k][i]["alt"]) / int(df[k][i]["depth"])
                            #weight [k] = math.pow(10, step.membership_p[k][j])
                            weight [k] = step.membership_p_normalize[k][j]

                        step.mixture[i][j] = round( np.average( vaf[ makeone_index_i ], weights = weight [ makeone_index_i ]), 5) * 2

                        #lowvaf cluster인지 보고, 맞다면 step.mixture를 peak 기준으로 수정
                        step = lowvaf_and_peak ( step, i, j, sum_alt, sum_depth, [ind for ind in ind_list if df[ind][i]["alt"] != 0], vaf_list, "soft", **kwargs )                    
        

 
        if NUM_CLONE == 1:
            step.mixture = np.array([[1.0]] * kwargs["NUM_BLOCK"])
        sum_mixture = np.sum (step.mixture, axis = 1)      # prenormalization
        pre_mixture = copy.deepcopy (step.mixture)

        p_list, step, condition, kwargs = isparent.makeone(input_containpos, df, np_vaf, np_BQ, step, **kwargs)     # 첫 3개는 lenient로, 그 다음은 strict로 거른다


        # checkall_strict, checkall_lenient 판단
        kwargs["CHECKALLFROM"] = "Mstep.py"
        if kwargs["MAKEONE_STRICT"] in  [2, 3]:       
            step.checkall_lenient, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "lenient", np_vaf, **kwargs) 
            step.checkall_strict, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "strict", np_vaf, **kwargs) 
        else:  # strict의 기준을 더 빡세게 잡아도 됨
            step.checkall_lenient, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "lenient", np_vaf, **kwargs) 
            step.checkall_strict, sum_mixture, p, kwargs =  miscellaneous.checkallttest (step, "strict", np_vaf, **kwargs) 


        if p_list == []:   # if failed
            step.likelihood = step.likelihood_record[kwargs["STEP"]] = -9999999
            step.checkall_lenient = step.checkall_strict = False
            step.makeone_index = []
            if kwargs["VERBOSE"] >= 1:
                if step.makeone_experience == False:
                    print ("\t\t\tMstep.py : condition = {}\tcheckall_strict = {}, checkall_lenient = {}\tUnable to make {}  ({}) ".format ( condition, step.checkall_strict, step.checkall_lenient, kwargs["MAKEWHAT"], "\t".join(str(np.round(row, 2)) for row in step.mixture )  ) )
                else:
                    print ("\t\t\tMstep.py : condition = {}\tcheckall_strict = {}, checkall_lenient = {}\tUnable to make appropriate superclone-clone structure ({}) ".format ( condition, step.checkall_strict, step.checkall_lenient, "\t".join(str(np.round(row, 2)) for row in step.mixture )  ) )

        else:    # if succeedeed
            # if kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] ):     
            #     for i in range(NUM_BLOCK):
            #         sum = 0
            #         for j in range(NUM_CLONE):
            #             if j in step.makeone_index:            
            #                 sum = sum + step.mixture[i][j]           

            #         for j in range(NUM_CLONE):
            #             if j in step.makeone_index:      
            #                 step.mixture[i][j] = np.round(step.mixture[i][j] / sum, 2) if sum != 0 else 0
            # step.mixture = step.mixture * kwargs["MAKEWHAT"]


            kwargs["CHECKALLFROM"] = "Mstep.py"
            if kwargs["MAKEONE_STRICT"] in  [2, 3]:       # BioData에서는 극단적으로 lenient하게 잡아줘야 하니까
                step.checkall_lenient, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "lenient", np_vaf, **kwargs) 
                step.checkall_strict, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "strict", np_vaf, **kwargs) 
            else:
                step.checkall_lenient, sum_mixture, p, kwargs =  miscellaneous.checkall (step, "lenient", np_vaf, **kwargs) 
                step.checkall_strict, sum_mixture, p, kwargs =  miscellaneous.checkallttest (step, "strict", np_vaf, **kwargs) 


            if kwargs["VERBOSE"] >= 1:
                # if kwargs["STEP"] <= (kwargs["SOFT_LENIENT"] ) :     
                #     print ("\t\t\tMstep.py : condition = {}\tcheckall_strict = {}\tcheckall_lenient = {}\traw = ({}), normalized ({})".format(condition, step.checkall_strict,  step.checkall_lenient, pre_mixture.flatten (), step.mixture.flatten () )) 
                # else:
                print ("\t\t\tMstep.py : condition = {}, checkall_strict = {}, checkall_lenient = {}\t(step.mixture = {})".format(condition, step.checkall_strict, step.checkall_lenient, " ".join(str( np.round(row, 2) ) for row in pre_mixture )  ) ) 
                
        
            # step.membership_p_normalize = np.zeros((NUM_MUTATION, step.membership_p.shape[1]), dtype="float64")
            # for k in range(NUM_MUTATION):
            #     if k in step.fp_member_index:
            #         step.membership_p_normalize[k] = np.zeros(step.membership_p_normalize.shape[1], dtype="float64")  # Set  1 (FP_index) 0 0 0 0 0    
            #         step.membership_p_normalize[k][step.fp_index] = 1
            #     else:
            #         step.membership_p_normalize[k] = np.round(np.power(10, step.membership_p[k]) / np.power(10, step.membership_p[k]).sum(axis=0, keepdims=1), 2)     # 분율을 따진다 ( = posterior allsample과 같음)
            #         if step.fp_index != -1: 
            #             step.membership_p_normalize[k][step.fp_index] = 0



        if  kwargs["VISUALIZATION"] == True:
            if (kwargs["NUM_BLOCK"] == 1):
                visualizationeachstep.drawfigure_1d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
            if (kwargs["NUM_BLOCK"] == 2):
                visualizationeachstep.drawfigure_2d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
            if (kwargs["NUM_BLOCK"] == 3):
                #visualizationeachstep.drawfigure_3d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
                visualizationeachstep.drawfigure_3d_SVD(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], sum_mixture, **kwargs)
            #subprocess.run ([ "cp -rf " +  kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "\(soft\)." + kwargs["IMAGE_FORMAT"] + " " + kwargs["COMBINED_OUTPUT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_ITER"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "\(soft\)." + kwargs["IMAGE_FORMAT"] ], shell = True)

    #############################################################################

    return step
