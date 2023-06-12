def main(df, np_vaf, np_BQ, step, option, **kwargs):  
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
            ind_list = np.where(step.membership == j)[0]   # Find the index  where membership == j
            for i in range(NUM_BLOCK):
                sum_depth, sum_alt = 0, 0
                for ind in ind_list:       # Summing depth and alt
                    if df[ind][i]["alt"] != 0:
                        sum_depth = sum_depth + df[ind][i]["depth"]
                        sum_alt = sum_alt + df[ind][i]["alt"]
                step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0   # Ideal centroid allocation

        
        
        previous_fp_index = step.fp_index
        step.makeone_index, p_list, step.fp_index = isparent.makeone(df, np_vaf, np_BQ, step, **kwargs)


        if step.fp_index != -1:  # If FP is present
            step.includefp = True
            step.fp_member_index = list(np.where(step.membership == step.fp_index)[0])
        else:   # If FP is absent
            step.includefp = False
            step.fp_member_index = []




        if step.makeone_index == []:
            step.likelihood = -9999999
            step.makeone_prenormalization = False
            print ("\t\t\tUnavailable to make 1 (Mstep.py)", end = "\t"  )
            print( ",".join( str(row) for row in step.mixture ))

        if kwargs["STEP"] <= 4:
            if step.makeone_index != []:
                step.makeone_prenormalization, sum_mixture =  EMhard.checkall (step, **kwargs) 
                print ("\t\t\tstep.makeone_prenormalization = {} (Mstep.py)".format(step.makeone_prenormalization), end = "\t") 
                print(" ".join(str(row) for row in sum_mixture ))

        
                for i in range(NUM_BLOCK):
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:   
                            sum = sum + step.mixture[i][j]
                    step.mixture[i] = np.round( step.mixture[i] / sum, 4) if sum != 0 else 0   # If sum = 0, let mixture = 0
                
        
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
            if j not in step.makeone_index:
                for i in range(NUM_BLOCK):
                    # Summing all depth and alt
                    sum_depth, sum_alt = 0, 0
                    for ind in np.where(step.membership == j)[0]:
                        if df[ind][i]["alt"] != 0:
                            sum_depth = sum_depth + df[ind][i]["depth"]
                            sum_alt = sum_alt + df[ind][i]["alt"]
                    
                    step.mixture[i][j] = round((sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0

            elif j in step.makeone_index:  
                for i in range(NUM_BLOCK):   # Calculate the weighted mean
                    vaf, weight = np.zeros(NUM_MUTATION, dtype="float"), np.zeros(NUM_MUTATION, dtype="float")
                    for k in range(NUM_MUTATION):
                        vaf[k] = int(df[k][i]["alt"]) / int(df[k][i]["depth"])

                        if step.membership[k] in step.makeone_index:
                            weight[k] = math.pow(10, step.membership_p[k][j])
                        #print ("{} : (clone {})   weight = {},  vaf = {}".format(k, step.membership[k], weight[k], vaf[k]))

                    step.mixture[i][j] = round(np.average(vaf[makeone_index_i], weights=weight[makeone_index_i]), 4) * 2

        
        # Normalize to make 1
        if kwargs["adjustment"] in ["True", "true", True]:
            for i in range(NUM_BLOCK):
                sum = np.sum(step.mixture[i])
                
                step.mixture[i] = np.round(step.mixture[i] / sum, 4) if sum != 0 else 0

        elif kwargs["adjustment"] in ["Half", "half"]:
            # step.makeone_index, p_list, step.fp_index = isparent.makeone (step, **kwargs)

            if step.fp_index != -1:  # If FP is present
                step.includefp = True
                step.fp_member_index = list(np.where(step.membership == step.fp_index)[0])
            else:   # If FP is absent
                step.includefp = False
                step.fp_member_index = []

            if NUM_CLONE == 1:
                step.mixture = np.array([[1.0]] * kwargs["NUM_BLOCK"])

            if step.makeone_index == []:
                step.likelihood = -9999999
                step.makeone_prenormalization = False
                

        if kwargs["STEP"] <= 5:
            if step.makeone_index != []:
                for i in range(NUM_BLOCK):
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:            
                            sum = sum + step.mixture[i][j]
                    

                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:      
                            step.mixture[i][j] = np.round(step.mixture[i][j] / sum, 4) if sum != 0 else 0

        #print ("\t\tc. Mixture (after normalization) : {}". format(list(step.mixture)))

        
        step.membership_p_normalize = np.zeros((NUM_MUTATION, step.membership_p.shape[1]), dtype="float64")
        for k in range(NUM_MUTATION):
            if k in step.fp_member_index:
                step.membership_p_normalize[k] = np.zeros(step.membership_p_normalize.shape[1], dtype="float64")  # Set  1 (FP_index) 0 0 0 0 0 
    
                step.membership_p_normalize[k][step.fp_index] = 1
            else:
                step.membership_p_normalize[k] = np.round(np.power(10, step.membership_p[k])/np.power(10, step.membership_p[k]).sum(axis=0, keepdims=1), 2)  
                if step.fp_index != -1: 
                    step.membership_p_normalize[k][step.fp_index] = 0

        if (kwargs["NUM_BLOCK"] == 1):
            visualizationeachstep.drawfigure_1d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft).pdf", **kwargs)
        if (kwargs["NUM_BLOCK"] == 2):
            visualizationeachstep.drawfigure_2d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft).pdf", **kwargs)

    #############################################################################

    return step
