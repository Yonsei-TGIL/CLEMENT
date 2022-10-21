from turtle import color


def np_vaf_extract(df):
    import numpy as np

    NUM_BLOCK = len(df[0])
    NUM_MUTATION = len(df)

    np_vaf = np.zeros((NUM_MUTATION, NUM_BLOCK), dtype = 'float')
    for row in range (NUM_MUTATION):
        for col in range (NUM_BLOCK):
            if df[row][col]["depth"] == 0:
                np_vaf[row][col] = 0
            else:    
                np_vaf[row][col] = round (df[row][col]["alt"] / df[row][col]["depth"] , 3)

    return np_vaf



def initial_kmeans (np_vaf, KMEANS_CLUSTERNO):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.cluster import KMeans
    import palettable

    NUM_BLOCK = np_vaf.shape[1]
    NUM_MUTATION = np_vaf.shape[0]

    tabl = palettable.tableau.Tableau_20.mpl_colors
    colorlist = [i for i in tabl]

    kmeans = KMeans(n_clusters=KMEANS_CLUSTERNO, init='k-means++', max_iter=100,random_state=0)  # model 생성
    kmeans.fit(np_vaf)  # model에 집어넣어줄 dataframe 지정
    membership_kmeans = kmeans.labels_     # K means가 predict 한것

    #plt.figure(figsize=(6, 6))
    #plt.scatter (np_vaf[:,0] * 2, np_vaf[:,1] * 2, alpha = 1 , color = [colorlist[k] for k in membership_kmeans])

    mixture_kmeans = np.zeros (( NUM_BLOCK, KMEANS_CLUSTERNO), dtype = "float")
    for j in range(KMEANS_CLUSTERNO):
        for i in range (NUM_BLOCK):
            mixture_kmeans[i][j] = round(np.mean(np_vaf[membership_kmeans == (j)][:,i] * 2), 3)
        #plt.scatter (mixture_kmeans[0][j], mixture_kmeans[1][j], marker = '*', color = colorlist[j], edgecolor = "black", s = 200, label = "clone " + str(j))

    #plt.legend()
    
    return mixture_kmeans


def kmeans_best (np_vaf, **kwargs ):
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy, math
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_samples, silhouette_score
    import palettable
    
    tabl = palettable.tableau.Tableau_20.mpl_colors
    colorlist = [i for i in tabl]

    KMEANS_CLUSTERNO = kwargs ["KMEANS_CLUSTERNO"]
    NUM_BLOCK = np_vaf.shape[1]
    NUM_MUTATION = np_vaf.shape[0]
    Input_B = 5
    Gap_list, Std_list, S_list = np.zeros ( KMEANS_CLUSTERNO + 1, dtype ="float"), np.zeros (KMEANS_CLUSTERNO + 1, dtype ="float"), np.zeros (KMEANS_CLUSTERNO + 1, dtype ="float")

    for n_clusters in range (1, KMEANS_CLUSTERNO + 1):
        kmeans = KMeans(n_clusters=n_clusters, init='k-means++', max_iter=100,random_state=0)  # model 생성
        kmeans.fit(np_vaf)  # model에 집어넣어줄 dataframe 지정
        
        membership = kmeans.labels_     # K means가 predict 한것
        mixture = np.zeros (( NUM_BLOCK, n_clusters ), dtype = "float")
        for j in range( n_clusters ):
            for i in range (NUM_BLOCK):
                mixture[i][j] = round(np.mean(np_vaf[membership == (j)][:,i] * 2), 3)
                
        plt.figure(figsize=(6, 6))
        plt.scatter (np_vaf[:,0] * 2, np_vaf[:,1] * 2, alpha = 1 , color = [colorlist[k] for k in membership])
        #plt.savefig ( kwargs["CLEMENT_DIR"] + "/Kmeans/Kmeans.clone"  +   str (n_clusters) + ".jpg" )
           
        #1. Intra cluster variation (Wk)
        Wk = 0
        for k in range(NUM_MUTATION):
            j = membership [k]
            Wk = Wk + math.pow (  scipy.spatial.distance.euclidean(np_vaf[k] * 2, mixture[:, j]),  2  )   # Sum of square
        Wk = round(math.log10(Wk), 3)
        print ("n_clusters={}\tMy Clustering\tWk  : {}" .format(n_clusters, Wk))

        #2. Random generation & ICC (Wkb)
        Wkb_list = []
        for b in range (Input_B):
            reference_uniform = np.random.random_sample ( size = (NUM_MUTATION, NUM_BLOCK ) )    # Outlier라면 계산에서 뺴주자
            # reference_f = f_distribution ( NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE]) , NUM_BLOCK   )
            # reference_multiply = multiply_npvaf ( NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE]) , NUM_BLOCK , np_vaf, 
            #                                                                         sorted (set( range(0, NUM_MUTATION) ) - set(cluster.outlier_index_record [NUM_CLONE])) )
            reference_np = reference_uniform

            #kmeans = KMeans(n_clusters=NUM_CLONE - int (cluster.includeoutlier_record [NUM_CLONE]), init='k-means++', max_iter=10, random_state=0)  # model 생성,  outlier (fp) 가 있다면 n_cluster를 하나 줄이는게 맞다
            #kmeans = KMeans(n_clusters=n_clusters,  init='k-means++', max_iter=100,random_state=0)  # model 생성
            kmeans.fit(reference_np)  # model에 집어넣어줄 nparray 지정
            Wkb_list.append ( round (math.log10(kmeans.inertia_), 3) )
            #drawfigure (reference_np, kmeans.labels_,  kwargs["CLEMENT_DIR"] + "/Kmeans/Kmeans.clone"  +   str (n_clusters) + "." + str(b) + ".jpg")

        Gap_list [n_clusters] = round ( np.mean(Wkb_list) - Wk, 3)
        Std_list [n_clusters] = round ( np.std (Wkb_list), 3)
        S_list [n_clusters] = round (Std_list[n_clusters] * math.sqrt(1 + Input_B) , 3 )
        
        print ("\tRandom noise (B = {}) : \tmean Wkb = {}\tsdk = {}\tsk (sdk * sqrt ({})) = {}\n\tGap (Wkb - Wk) = {}\n".format (Input_B, round( np.mean(Wkb_list), 3) , Std_list[n_clusters], Input_B + 1, S_list[n_clusters]  ,  Gap_list[n_clusters] ))


    # if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
    #     print ("Gap list : {}\nS list : {}\n".format(Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1], S_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1] ))

    # # Max Gap Number
    # maxmaxmax_NUM_CLONE  = np.argmax ( Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1]  ) + kwargs["NUM_CLONE_TRIAL_START"]
    
    # Gap_list_df = pd.DataFrame ( Gap_list ).sort_values ( by = 0, ascending = False)
    # Gap_list_df = Gap_list_df [ Gap_list_df[0] != 0]
    # Gap_list_index = []

    # for i in range( len(Gap_list_df) ):
    #     print ("Gap statistics method (max Gap): {}th optimal NUM_CLONE = {}".format(i  + 1, Gap_list_df.index[i]  ))
    #     Gap_list_index.append ( Gap_list_df.index[i]  )

    # return Gap_list_index           # Gap 큰 순서대로 NUM_CLONE의 list를 보내주자


    maxmaxmax_NUM_CLONE  = 1
    for NUM_CLONE in range ( 0, KMEANS_CLUSTERNO + 1 ):
        if Gap_list[NUM_CLONE] >= Gap_list[NUM_CLONE + 1] - S_list [NUM_CLONE +1]:
            maxmaxmax_NUM_CLONE  = NUM_CLONE
            print ("Gap statistics method : optimal NUM_CLONE = {}".format(NUM_CLONE))
            break

    print ("Gap statistics method : optimal NUM_CLONE = {}".format( maxmaxmax_NUM_CLONE ))
    

    

def set_initial_parameter(df, mixture_kmeans, **kwargs):
    import random
    import numpy as np
    
    initial_mixture = np.zeros((kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"]), dtype="float")
    initial_mixture = mixture_kmeans [:, random.sample (range( kwargs["KMEANS_CLUSTERNO"]), kwargs["NUM_CLONE"])]

    return initial_mixture


def movedcolumn (mixture_hard, mixture_soft):
    import scipy.spatial
    import math
    
    if mixture_hard.shape[1] != mixture_soft.shape[1]:
        return []
    
    NUM_CLONE = mixture_hard.shape[1]
    NUM_BLOCK = mixture_hard.shape[0]
    threshold_dist =  (math.sqrt ( NUM_BLOCK ) * 0.03)

    col_list = []
    for i in range ( NUM_CLONE ):
        #print ( "{} ; {}".format (i,  float(scipy.spatial.distance.euclidean( mixture_hard[ : , i ] , mixture_soft [ : , i ]  )) ) )
        if float(scipy.spatial.distance.euclidean( mixture_hard[ : , i ] , mixture_soft [ : , i ]  )) > threshold_dist:
            col_list.append(i)

    return col_list



def GoStop(step, **kwargs):
    import numpy as np
    previous_No = 5

    # 첫 5번 정도는 판단 유예
    if kwargs["STEP"] > previous_No:
        for p in range( kwargs["STEP_TOTAL"] - previous_No,  kwargs["STEP_TOTAL"] ):       # 지난 5개의 step 기록을 보고 비교해본다
            if np.sum(np.equal(step.membership, step.membership_record[p])) >= int(kwargs["NUM_MUTATION"] * 0.995):           # membership이 99% 이상 겹치면 stop 해주자
                #print ( "지난 {0} 번째와 membership 변화가 없어서 stop\n".format(p) )
                return "Stop"

            # mixture가 완전히 똑같으면 stop 해주자
            if (np.round(step.mixture, 2) == np.round( step.mixture_record[p] , 2)).all() == True:
                #print ("지난 {0} 번째와 mixture 변화가 없어서 stop\n".format(p))
                return "Stop"

        # likelihood가 더 나아질 기미가 없으면 stop 해주자
        i =  step.find_max_likelihood(0, kwargs["STEP"]) 
        if (step.likelihood) < 0.99 * step.likelihood_record [ i ]:
            #print ("Likelihood 기준으로 stop\n")
            return "Stop"

    return "Go"


def decision_elbow (cluster, **kwargs):
    import numpy as np
    from kneed import KneeLocator
    import matplotlib.pyplot as plt

    if kwargs["NUM_CLONE_TRIAL_START"] == kwargs["NUM_CLONE_TRIAL_END"]:       # 1개밖에 없으면  kneedle이 오류난다
        threshold_x = kwargs["NUM_CLONE_TRIAL_END"]
    else:            
        x = np.arange(start = kwargs["NUM_CLONE_TRIAL_START"], stop = kwargs["NUM_CLONE_TRIAL_END"] + 1)
        y = np.array( cluster.likelihood_record  [ kwargs["NUM_CLONE_TRIAL_START"] :  kwargs["NUM_CLONE_TRIAL_END"] + 1 ]  )

        for S in list(np.arange(0.1, 3.0, 0.5)):
            kneedle = KneeLocator(x, y, S=S, curve="concave", direction="increasing")
            if (kneedle.knee != None) & (kwargs["VERBOSE"] >= 0):
                kneedle.plot_knee()
                plt.show()
                break
        if kneedle.knee != None:             # 가장 정상적인 경우
            threshold_x = round(kneedle.knee)
        else:                                             # elbow 법에서 제대로 못 찾았을 경우  단순히 가장 높은 값으로 한다
            threshold_x = cluster.find_max_likelihood( kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] )
            print ("최적의 clone 값을 elbow법으로 찾지 못함 : {}".format(threshold_x))

    return [threshold_x]   # maxmaxmax_NUM_CLONE


def VAFdensitogram (np_vaf, output_suptitle, output_filename,  **kwargs):
    import matplotlib
    import seaborn as sns
    from scipy.stats import kde
    import numpy as np
    import palettable

    tabl = palettable.tableau.Tableau_20.mpl_colors
    colorlist = [i for i in tabl]

    fig, ax = matplotlib.pyplot.subplots (figsize = (6, 6))
    matplotlib.pyplot.rcParams["font.family"] = 'arial'
    matplotlib.pyplot.suptitle ("{}".format(output_suptitle), fontsize = 20)
    ax.set_xlabel("VAF", fontdict = {"fontsize" : 14})
    ax.set_ylabel("Density", fontdict = {"fontsize" : 14})
    max_y = 0

    for i in range(kwargs["NUM_BLOCK"]):
        x = np.linspace(0,1,101)
        kde_np_vaf_new = kde.gaussian_kde(np_vaf[:,i])
        y = kde_np_vaf_new(x)
        if max_y < np.max(y):
            max_y = np.max(y)
        ax.plot(x, y, color = colorlist[i], label = "sample {0}".format(i))
        
        x_median = np.median (np_vaf[:, i])
        y_median = y [ round (x_median * 100) ]
        if x_median > 0.4:
            ax.text( x_median, max_y * 1.15, "{0}\n→ Monoclonal".format(x_median), verticalalignment='top', ha = "center", fontsize = 15, color = "g")
        elif x_median > 0.25:
            ax.text( x_median, max_y * 1.15, "{0}\n→ Biclonal".format(x_median), verticalalignment='top', ha = "center", fontsize = 15, color = "g")
        else:
            ax.text( x_median, max_y * 1.15, "{}\n→ Polyclonal".format(x_median), verticalalignment='top', ha = "center", fontsize = 15, color = "g")
        ax.vlines( x = x_median, ymin = 0, ymax = y_median , color="g", linestyle = "--", label = "median VAF")

    ax.axis ([0,  1.01,  0,  max_y * 1.3])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth (3)
    ax.legend()



    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)
    
    return x_median



def drawfigure (np_vaf, membership, output_filename):
    import palettable
    import matplotlib
    import numpy as np
    import seaborn as sns

    tabl = palettable.tableau.Tableau_20.mpl_colors
    colorlist = [i for i in tabl]

    if np_vaf.shape[1] == 2:
        matplotlib.pyplot.figure(figsize=(6, 6))
        matplotlib.pyplot.axis([0,  1,  0,  1])
        matplotlib.pyplot.scatter (x = np_vaf [: ,0], y = np_vaf [:, 1] , color=[colorlist[k] for k in membership]  )
        matplotlib.pyplot.savefig(output_filename)
    elif np_vaf.shape[1] == 1:
        matplotlib.pyplot.figure(figsize=(6, 2))
        matplotlib.pyplot.xlim(0,  1)
        matplotlib.pyplot.yticks ( [] )
        matplotlib.pyplot.scatter (x = np_vaf [: ,0] , y = [0.5] * np_vaf.shape[0],  color=[colorlist[k] for k in membership]  )
        matplotlib.pyplot.savefig(output_filename)



def f_distribution (NUM_MUTATION_local, NUM_BLOCK_local):
    import numpy as np
    import random

    a = np.random.f (3, 10, (NUM_MUTATION_local * 5, NUM_BLOCK_local))

    f_max = np.max (a)
    a_list = []
    for i in range ( len(a) ):
        if np.all ( a [ i ,  : ] < ( f_max / 4 ) ) == True:
            a_list.append(i) 

    a = a [np.array (random.sample (a_list, NUM_MUTATION_local))]
    f_max = np.max(a)
    a = a / f_max

    return a


def multiply_npvaf ( NUM_MUTATION_local, NUM_BLOCK_local, np_vaf, nonoutlier_index ): 
    import numpy as np
    import random

    rr = 4
    rrr = [80, 100, 100, 120]

    model_np_vaf = np.zeros ( (NUM_MUTATION_local * rr , NUM_BLOCK_local),  dtype = "float")

    for k in range (NUM_MUTATION_local):
        for i in range (NUM_BLOCK_local):
            for r_index in range (rr):
                while True:
                    r =  random.randint(70, rrr[r_index]) / 100
                    model_np_vaf [NUM_MUTATION_local * r_index +  k][i] = np_vaf[ nonoutlier_index[k] ][i] * r * 2
                    if np_vaf[ nonoutlier_index[k] ][i] * r * 2 <= 2:          #  (논란) 곱하고 나누고 늘려줘도 1 이하인 것들만 통과시켜준다
                        break
            # if np_vaf[k][0] < 0.04:
            #     print (model_np_vaf [NUM_MUTATION_local * 1 +  k], model_np_vaf [NUM_MUTATION_local * 2 +  k])

    model_sel_np_vaf = model_np_vaf [ random.sample( range(len(model_np_vaf)), NUM_MUTATION_local ) , : ]   # NUM_MUTATION 개만 골라서 return 해준다


    #print (np_vaf[0:20], "\n\n", model_sel_np_vaf [0:20])
    return model_sel_np_vaf


def decision_gapstatistics (cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.cluster import KMeans

    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]

    Input_B = 20
    Gap_list, Std_list, S_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float"), np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float"), np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")
    maxmaxmax_NUM_CLONE  = 0

    
    for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
            print ("NUM_CLONE == {}".format(NUM_CLONE))
            print ("\tIncludeOutlier = {} (fp_index = {})\t\t-> {}".format( cluster.includeoutlier_record [NUM_CLONE], cluster.fp_index_record [NUM_CLONE],  NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE])   ))

        if (cluster.stepindex_record[NUM_CLONE] == 0) | (cluster.likelihood_record[NUM_CLONE] < -9999990):  # 그 어떤 답도 찾을 수 없어서 stepindex = 0인 경우
            Gap_list [NUM_CLONE] = float("-inf")

        else:
            if maxmaxmax_NUM_CLONE == 0:
                maxmaxmax_NUM_CLONE = NUM_CLONE
                
            membership = cluster.membership_record[NUM_CLONE]
            mixture = cluster.mixture_record[NUM_CLONE]

            #1. Intra cluster variation (Wk)
            Wk = 0
            for k in range(NUM_MUTATION):
                j = membership [k]
                if (cluster.includeoutlier_record [NUM_CLONE] == True)  &  (k in cluster.outlier_index_record [NUM_CLONE] ):   # FP라면 계산에서 뺴주자
                    continue
                Wk = Wk + math.pow (  scipy.spatial.distance.euclidean(np_vaf[k] * 2, mixture[:, j]),  2)   # Sum of square
            Wk = round(math.log10(Wk), 3)
            if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
                print ("\tMy Clustering\tWk  : {}" .format(Wk))


            #2. Random generation & ICC (Wkb)
            Wkb_list = []
            for b in range (Input_B):
                reference_uniform = np.random.random_sample(size=(NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE]) , NUM_BLOCK ))    # Outlier라면 계산에서 뺴주자
                # reference_f = f_distribution ( NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE]) , NUM_BLOCK   )
                reference_multiply = multiply_npvaf ( NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE]) , NUM_BLOCK , np_vaf, 
                                                                                        sorted (set( range(0, NUM_MUTATION) ) - set(cluster.outlier_index_record [NUM_CLONE])) )
                reference_np = reference_multiply
                # if NUM_BLOCK == 0:
                #     reference_np = reference_uniform

                #kmeans = KMeans(n_clusters=NUM_CLONE - int (cluster.includeoutlier_record [NUM_CLONE]), init='k-means++', max_iter=10, random_state=0)  # model 생성,  outlier (fp) 가 있다면 n_cluster를 하나 줄이는게 맞다
                if cluster.includeoutlier_record [NUM_CLONE] == False:
                    kmeans = KMeans(n_clusters=NUM_CLONE, init = cluster.mixture_record [NUM_CLONE].transpose() , max_iter=3, random_state=0)  # model 생성,  outlier (fp) 가 있다면 n_cluster를 하나 줄이는게 맞다
                else:   # FP 있을 경우 빼고 돌린다
                    kmeans = KMeans(n_clusters=NUM_CLONE - 1, init = np.delete(cluster.mixture_record [NUM_CLONE].transpose(), cluster.fp_index_record [NUM_CLONE], axis = 0), max_iter=3, random_state=0)  # model 생성,  outlier (fp) 가 있다면 n_cluster를 하나 줄이는게 맞다

                kmeans.fit(reference_np)  # model에 집어넣어줄 nparray 지정
                Wkb_list.append ( round (math.log10(kmeans.inertia_), 3) )
                drawfigure (reference_np, kmeans.labels_,  kwargs["CLEMENT_DIR"] + "/Kmeans/Kmeans.clone"  +   str (NUM_CLONE) + "." + str(b) + ".jpg")

            Gap_list [NUM_CLONE] = round ( np.mean(Wkb_list) - Wk, 3)
            Std_list [NUM_CLONE] = round ( np.std (Wkb_list), 3)
            S_list [NUM_CLONE] = round (Std_list[NUM_CLONE] * math.sqrt(1 + Input_B) , 3 )
            

            if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
                print ("\tRandom noise (B = {}) : \tmean Wkb = {}\tsdk = {}\tsk (sdk * sqrt ({})) = {}\n\tGap (Wkb - Wk) = {}\n\tPosterior = {}".format (Input_B, round( np.mean(Wkb_list), 3) , Std_list[NUM_CLONE], Input_B + 1, S_list[NUM_CLONE]  ,  Gap_list[NUM_CLONE], round (cluster.likelihood_record[NUM_CLONE]) ))

    if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
        print ("Gap list : {}\nS list : {}\n".format(Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1], S_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1] ))


    # Max Gap Number
    if kwargs ["GAP"] in  ["new", "sil"] :
        Gap_list_index = []
        maxmaxmax_NUM_CLONE  = np.argmax ( Gap_list [kwargs["NUM_CLONE_TRIAL_START"] : kwargs["NUM_CLONE_TRIAL_END"] + 1]  ) + kwargs["NUM_CLONE_TRIAL_START"]
        
        Gap_list_df = pd.DataFrame ( Gap_list ).sort_values ( by = 0, ascending = False)
        Gap_list_df = Gap_list_df [ Gap_list_df[0] != 0]

        for i in range( len(Gap_list_df) ):
            print ("Gap statistics method (max Gap): {}th optimal NUM_CLONE = {}".format(i  + 1, Gap_list_df.index[i]  ))
            Gap_list_index.append ( Gap_list_df.index[i]  )

        return Gap_list_index           # Gap 큰 순서대로 NUM_CLONE의 list를 보내주자


    # Original Gap Number 
    elif kwargs ["GAP"] == "old":
        for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"]):
            if Gap_list [NUM_CLONE] != float("-inf"):
                if Gap_list[NUM_CLONE] >= Gap_list[NUM_CLONE + 1] - S_list [NUM_CLONE +1]:
                    maxmaxmax_NUM_CLONE  = NUM_CLONE
                    if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
                        print ("Original gap statistics method : optimal NUM_CLONE = {}".format(NUM_CLONE))
                    return [maxmaxmax_NUM_CLONE]

        print ("Original Gap statistics method : optimal NUM_CLONE = {}".format( maxmaxmax_NUM_CLONE  ))
        return [ maxmaxmax_NUM_CLONE  ]            # 맨 마지막은 그 다음과 비교할 수 없으니...



def decision_silhouette (cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.metrics import silhouette_samples, silhouette_score
    
    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]
    Silhouette_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")
    
    print ("\n\n------------------- Silhouette Method -------------------\n")
    
    for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 0):
            print ("NUM_CLONE == {}".format(NUM_CLONE))
            print ("\tIncludeOutlier = {} (fp_index = {})\t\t-> {}".format( cluster.includeoutlier_record [NUM_CLONE], cluster.fp_index_record [NUM_CLONE],  NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE])   ))

        if (cluster.stepindex_record[NUM_CLONE] == 0) | (cluster.likelihood_record[NUM_CLONE] < -9999990):  # 그 어떤 답도 찾을 수 없어서 stepindex = 0인 경우
            Silhouette_list [NUM_CLONE] = 0
            
        else:
            membership = cluster.membership_record[NUM_CLONE]
            membership_fpexclude = np.array ([membership[i]   for i in set (range (0, len(membership)) ) - set( cluster.outlier_index_record [NUM_CLONE] ) ])
            np_vaf_fpexclude = np.array ([np_vaf[i]   for i in set (range (0, len(np_vaf)) ) - set( cluster.outlier_index_record [NUM_CLONE] ) ])
                
            silhouette_score_alldata = silhouette_samples(np_vaf_fpexclude , membership_fpexclude)   # (500,1 )  500
            Silhouette_list [NUM_CLONE] = np.mean (silhouette_score_alldata)


    arg_list = [ list(Silhouette_list).index(i) for i in sorted(Silhouette_list, reverse=True)][:2]
    print ("Silhouette_list : {}\narg_list : {}".format  (Silhouette_list, arg_list) )
    
    return arg_list



def decision_BIC (cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.cluster import KMeans

    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]

    BIC_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")


    for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
            print ("NUM_CLONE == {}".format(NUM_CLONE))
            print ("\tIncludeOutlier = {} (fp_index = {})\t\t-> {}".format( cluster.includeoutlier_record [NUM_CLONE], cluster.fp_index_record [NUM_CLONE],  NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE])   ))

        if (cluster.stepindex_record[NUM_CLONE] == 0) | (cluster.likelihood_record[NUM_CLONE] < -9999990):  # 그 어떤 답도 찾을 수 없어서 stepindex = 0인 경우
            BIC_list [NUM_CLONE] = float("-inf")

        else:
            # if cluster.includeoutlier_record [NUM_CLONE] == False:
            #     BIC_list [NUM_CLONE] = (-2) * cluster.likelihood_record [NUM_CLONE] + (  NUM_CLONE  * math.log (NUM_MUTATION))
            # else :
            #     BIC_list [NUM_CLONE] = (-2) * cluster.likelihood_record [NUM_CLONE] + (  (NUM_CLONE - 1)  * math.log (NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE]) ) )


            membership = cluster.membership_record[NUM_CLONE]
            mixture = cluster.mixture_record[NUM_CLONE]

            Nc = np.unique ( membership, return_counts = True) [1]
            #print ("Nc : {}  (shape = {})".format (Nc, Nc.shape))
            Vc = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = "float")     # 각 cluster 별로 variance를 계산한 것
            V = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = "float")       # 모든 mutation 별로 variance를 계산한 것
            LL = np.zeros ( NUM_CLONE, dtype = "float" )       # 모든 mutation 별로 variance를 계산한 것


            for j in range (NUM_CLONE):
                if j == cluster.fp_index_record [NUM_CLONE]:
                    continue

                NUM_CLONE_index = np.where (membership == j)[0]
                for i in range (NUM_BLOCK):
                    m = np.mean (np_vaf [NUM_CLONE_index][i] ) * 2
                    inertia = 0
                    for k in NUM_CLONE_index:
                        inertia = inertia + math.pow ( (np_vaf [k][i] * 2) - m, 2 )
                    Vc [i][j] = inertia  /  ( len(NUM_CLONE_index) - 1 )

                    m = np.mean (np_vaf [:, i] ) * 2
                    inertia = 0 
                    for k in range (NUM_MUTATION):
                        inertia = inertia + math.pow ( (np_vaf [k][i] * 2) - m, 2 )
                    V [i][j] = inertia / (NUM_MUTATION - 1)

                    #print ("({},{}) VC[i][j] ={}\tV[i][j] = {}".format (i , j , Vc[i][j], V[i][j]))

            for j in range (NUM_CLONE):
                if j == cluster.fp_index_record [NUM_CLONE]:
                    continue

                csum = 0
                for i in range (NUM_BLOCK):
                    csum = csum + math.log (  ( Vc[i][j] + V[i][j] )  / 2   )
                LL[j] = (-1) * Nc [j]  * csum
            #print ("LL = {}".format (LL))

            if cluster.includeoutlier_record [NUM_CLONE] == False:
                BIC_list [NUM_CLONE] = (-2) * np.sum(LL) + ( 2 * NUM_CLONE * NUM_BLOCK * math.log (NUM_MUTATION))
            else:
                BIC_list [NUM_CLONE] = (-2) * np.sum(LL) + ( 2 * (NUM_CLONE - 1) * NUM_BLOCK * math.log (NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE]) ))

    print (BIC_list)

    BIC_list_df = pd.DataFrame ( BIC_list ).sort_values ( by = 0, ascending = False)
    BIC_list_df = BIC_list_df [ BIC_list_df[0] != 0]
    BIC_list_index = []

    for i in range( len(BIC_list_df) ):
        print ("BIC statistics method (max BIC): {}th optimal NUM_CLONE = {}".format(i  + 1, BIC_list_df.index[i]  ))
        BIC_list_index.append ( BIC_list_df.index[i]  )
   

def XieBeni_calcualtor (cluster, np_vaf, **kwargs):
    import math, scipy
    from sklearn.cluster import KMeans

    if kwargs["NUM_CLONE"] == 1:
        return 1111111111111


    min = float("inf")
    for j1 in range (kwargs["NUM_CLONE"]):   # fp index는 빼고 centroid 돌아야지
        for j2 in range (j1 + 1, kwargs["NUM_CLONE"]) :
            if (j1 == cluster.fp_index_record [ kwargs["NUM_CLONE"] ]) | (j2 == cluster.fp_index_record [ kwargs["NUM_CLONE"] ]) :
                continue

            if min > math.pow (scipy.spatial.distance.euclidean(cluster.mixture_record[kwargs["NUM_CLONE"]][:, j1], cluster.mixture_record[kwargs["NUM_CLONE"]][:, j2]), 2):
                min = math.pow (scipy.spatial.distance.euclidean(cluster.mixture_record[kwargs["NUM_CLONE"]][:, j1], cluster.mixture_record[kwargs["NUM_CLONE"]][:, j2]), 2)

    v = 0
    for j in range (kwargs["NUM_CLONE"])  :   # fp index는 빼고 centroid 돌아야지
        if (j == cluster.fp_index_record [ kwargs["NUM_CLONE"] ]) :
            continue

        for k in range (kwargs["NUM_MUTATION"]):
            if k in cluster.outlier_index_record[kwargs["NUM_CLONE"]]:
                continue
            v1 = math.pow( cluster.membership_p_record[kwargs["NUM_CLONE"]][k, j], 2)
            v2 = math.pow ( (scipy.spatial.distance.euclidean(cluster.mixture_record[kwargs["NUM_CLONE"]][:, j], np_vaf[k] * 2)) , 2)
            v = v + (v1 * v2)

    if min == 0:       #분모가 0이 되는 경우
        return 1111111111111 - kwargs["NUM_CLONE"]
    else:
        return (v / (min * ( kwargs["NUM_MUTATION"] - len(cluster.outlier_index_record[kwargs["NUM_CLONE"]]) )))



def decision_XieBeni(cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.cluster import KMeans


    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]
 
    XieBeni_list = np.zeros (kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype ="float")

    for NUM_CLONE in range (kwargs["NUM_CLONE_TRIAL_START"], kwargs["NUM_CLONE_TRIAL_END"] + 1):
        if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
            print ("NUM_CLONE == {}".format(NUM_CLONE))
            print ("\tIncludeOutlier = {} (fp_index = {})\t\t-> {}".format( cluster.includeoutlier_record [NUM_CLONE], cluster.fp_index_record [NUM_CLONE],  NUM_MUTATION - len(cluster.outlier_index_record [NUM_CLONE])   ))

        if (cluster.stepindex_record[NUM_CLONE] == 0) | (cluster.likelihood_record[NUM_CLONE] < -9999990):  # 그 어떤 답도 찾을 수 없어서 stepindex = 0인 경우
            XieBeni_list [NUM_CLONE] = float("inf")
        else:
            kwargs ["NUM_CLONE"] = NUM_CLONE
            XieBeni_list[NUM_CLONE] = XieBeni_calcualtor (cluster, np_vaf, **kwargs)

        if (kwargs["VERBOSE"] in ["True", "T"]) | (kwargs["VERBOSE"] >= 1):
            print ("\tXie-Beni index : {}".format (XieBeni_list[NUM_CLONE]))

    # Max XieBeni Number
        
    XieBeni_list_df = pd.DataFrame ( XieBeni_list ).sort_values ( by = 0, ascending = True)
    XieBeni_list_df = XieBeni_list_df [ XieBeni_list_df[0] != 0]
    XieBeni_list_index = []

    for i in range( len(XieBeni_list_df) ):
        print ("XieBeni method (min): {}th optimal NUM_CLONE = {}".format(i  + 1, XieBeni_list_df.index[i]  ))
        XieBeni_list_index.append ( XieBeni_list_df.index[i]  )

    return XieBeni_list_index 



def decision_max(cluster, np_vaf, **kwargs):
    import pandas as pd
    import numpy as np
    import math, scipy
    from sklearn.cluster import KMeans


    NUM_MUTATION, NUM_BLOCK = kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]

    likelihood_list = np.array ( cluster.likelihood_record )
    num_child_list = np.zeros ( kwargs["NUM_CLONE_TRIAL_END"] + 1, dtype = "int" )
    for j in range (0, kwargs["NUM_CLONE_TRIAL_END"] + 1):
        num_child_list [j] = len (cluster.makeone_index_record [j])

    candidate_df = pd.DataFrame ( np.concatenate ( ([likelihood_list] ,[num_child_list]) , axis = 0) ).transpose().sort_values ( by = 0, ascending = False)
    candidate_df.columns = ["likelihood", "NUM_CHILD"]
    candidate_df = candidate_df.astype ({"NUM_CHILD" : "int"})
    candidate_df ["NUM_CLONE"] = candidate_df.index
    candidate_df = candidate_df.reset_index(drop = True)


    candidate_df_after = pd.DataFrame (columns = candidate_df.columns)
    check = []
    check_index = []

    for k in range ( candidate_df.shape [0] ) :
        if candidate_df.loc[k, "NUM_CHILD"] not in set(check):
            if candidate_df.loc[k, "likelihood"] == float("-inf"):
                candidate_df_after = pd.concat( [ candidate_df_after, candidate_df.loc[ list ( set(list(range( candidate_df.shape[0] ))) - set(check_index) ) ] ])    
                candidate_df_after = candidate_df_after.reset_index(drop = True)
                break

            temp_df = candidate_df [ candidate_df["NUM_CHILD"] ==candidate_df.loc[k, "NUM_CHILD"] ].sort_values ( by = "NUM_CLONE", ascending = True) 
            check_index.append(temp_df.index [0])
            temp_df = temp_df.reset_index(drop = True)
            
            candidate_df_after = pd.concat( [ candidate_df_after, temp_df.loc[ [0] ] ])
            check.append  (candidate_df.loc[k, "NUM_CHILD"])
        
    candidate_df_after = candidate_df_after [ candidate_df_after["NUM_CLONE"] != 0 ]   # NUM_CLONE =0 은 지워준다
    candidate_df_after = candidate_df_after.astype ({"NUM_CHILD" : "int",  "NUM_CLONE" : "int"})
    candidate_df_after = candidate_df_after.reset_index(drop = True)


    print ("")
    for i in range( len(candidate_df_after) ):
        print ("max likelihood method : {}th optimal NUM_CLONE = {}\tNUM_CHILD = {}".format(i  + 1, candidate_df_after["NUM_CLONE"].iloc [i],  candidate_df_after["NUM_CHILD"].iloc [i]  ))

    return list (candidate_df_after ["NUM_CLONE"] )