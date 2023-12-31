import palettable
import matplotlib
import extract
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import kde
from sklearn.decomposition import TruncatedSVD, PCA


def drawfigure_1d_hard(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    
    # if (bunch.makeone_index == []) | (bunch.makeone_index == None) | (bunch.checkall_strict == False):  
    #     matplotlib.pyplot.style.use("Solarize_Light2")
    # else:
    matplotlib.pyplot.style.use("seaborn-white")


    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includefp == True) & (bunch.fp_index != -1):
        colorlist[bunch.fp_index] = Gr_10[8]


    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    fig.subplots_adjust ( bottom = 0.15, top = 0.85)

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    max_y = 0

    x = np.linspace(0, 2, 200)

    for j in sorted(list(set(bunch.membership))):     
        np_vaf_new_index, np_vaf_new = extract.npvaf(np_vaf, bunch.membership, j)           # extract membership1  == clone_num
        mean_x = np.mean (np_vaf_new [:, 0] * 2)

        if (len ( set(np_vaf_new [:, 0]) ) <= 1):
            continue

        kde_np_vaf_new = kde.gaussian_kde(np_vaf_new[:, 0] * 2)
        weight = len(np_vaf_new) / kwargs["NUM_MUTATION"]  
        y = kde_np_vaf_new(x) * weight

        if max_y < np.max(y):
            max_y = np.max(y)


        if j == bunch.fp_index:  # FP를 따로 그려줌
            ax.plot(x, y, color = Gr_10[10], linewidth=5, label="FP ({}, n = {})".format( np.round( bunch.mixture[0, j], 2), np.bincount(bunch.membership)[j]))
            ax.text(0, np.max(y) * 1.1, "{}\n(n={})".format( round(mean_x, 2),  np.bincount(bunch.membership)[j]), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})
            #ax.axvline( x = 0,  ymin = 0, ymax = np.max(y) / (ax.get_ylim()[1] - ax.get_ylim()[0]) , linestyle = '--', linewidth = 2, color = Gr_10[10] )
        elif j in bunch.makeone_index:
            if j not in bunch.lowvafcluster_index [0]:
                ax.plot(x, y, color = colorlist[j], linewidth=5, label="clone {} ({}, cdf = {}, n = {})".format(j, np.round( bunch.mixture[0, j], 2), bunch.cluster_cdf_list[0][j], np.bincount(bunch.membership)[j]) )
                ax.fill_between(x, y, color = colorlist[j], alpha = 0.5, label = None )
            else:   # most of the cases
                ax.plot(x, y, color = colorlist[j], linewidth=5, label="(lowvaf) clone {} ({}, cdf = {}, n = {})".format(j, np.round( bunch.mixture[0, j], 2), bunch.cluster_cdf_list[0][j], np.bincount(bunch.membership)[j]) )
            ax.text(bunch.mixture[0, j], np.max(y) * 1.1, "{}\n(n={})".format( np.round ( bunch.mixture[0,  j], 2),  np.bincount(bunch.membership)[j]), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})
            #ax.axvline( x = bunch.mixture[0, j],  ymin = 0, ymax = np.max(y) / (ax.get_ylim()[1] - ax.get_ylim()[0]) * weight, linestyle = '-', linewidth = 2, color = colorlist[j] )
        else:
            ax.plot(x, y, color=colorlist[j], linewidth=2, label="(parent) clone {} ({}, cdf = {}, n = {})".format(j, np.round( bunch.mixture[0, j], 2), bunch.cluster_cdf_list[0][j], np.bincount(bunch.membership)[j]), linestyle="-.")
            ax.text(bunch.mixture[0, j], np.max(y) * 1.1, "{}\n(n={})".format( np.round( bunch.mixture[0, j], 2) ,  np.bincount(bunch.membership)[j]), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})
            #ax.axvline( x = bunch.mixture[0, j],  ymin = 0, ymax = np.max(y) / (ax.get_ylim()[1] - ax.get_ylim()[0]) * weight, linestyle = '--', linewidth = 2, color = colorlist[j] )



    

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
        # Set the layout to tight to ensure the rectangle is not cut off
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
        matplotlib.pyplot.title("postnorm = {}, prenorm = {}, standard = {}\ndepth = {}, strict = {}, lenient = {}".format(list(np.round(np.sum(bunch.mixture[:, bunch.makeone_index], axis=1), 2)), np.round(sum_mixture_beforenormalization, 2), kwargs["MAKEONE_STANDARD"][0], kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize=12)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (raw) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
        matplotlib.pyplot.title("prenorm = {}, standard = {}\ndepth = {}, strict = {}, lenient = {}".format( np.round( sum_mixture_beforenormalization, 2), kwargs["MAKEONE_STANDARD"][0],  kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize=12)
        

    ax.set_xticks ( np.round (np.arange (0, 1.01, 0.2), 1) )
    ax.set_xlim( 0, 1.01 )
    ax.set_ylim( 0, max_y * 1.3 )

    xtlabels = list ( ax.get_xticks () )
    ax.set_xticklabels ( [ str(i) + "\n({})".format( int ( int(kwargs["MEAN_DEPTH"]) * float (i) * 0.5)  ) for i in xtlabels ])
    ax.set_yticks ([])
    

    ax.set_xlabel("VAF x 2\n(Alt)", fontdict={"fontsize": 14})
    ax.set_ylabel("Density", fontdict={"fontsize": 14})

    ax.legend(loc = 'best')
    matplotlib.pyplot.show()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)


################################################################################################################################################################################################


def drawfigure_1d_soft(bunch, np_vaf, output_filename, sum_mixture, **kwargs):
    import random
    from scipy.stats import kde

    membership = bunch.membership
    mixture = bunch.mixture
    membership_p_normalize = bunch.membership_p_normalize

    matplotlib.pyplot.style.use("seaborn-white")
        

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includefp == True) & (bunch.fp_index != -1):
        colorlist[bunch.fp_index] = Gr_10[8]

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]
    matplotlib.rcParams["text.color"] = "blue"

    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    fig.subplots_adjust ( bottom = 0.15, top = 0.85)
    for axis in ['left', 'right', 'top', 'bottom']:
        ax.spines[ axis ].set_color('blue')
    
    if (kwargs["STEP"] <= ( kwargs["SOFT_LENIENT"] )) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (lenient) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (raw) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)

    if kwargs["STEP"] <= (kwargs["SOFT_LENIENT"] ):  # normalization 한 것
        ax.set_title("prenorm = {}, standard = {}\ndepth = {}, strict = {}, lenient = {}".format(np.round (sum_mixture, 2) , kwargs["MAKEONE_STANDARD"][0],  kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize=12)
    else:
        ax.set_title("prenorm = {}, standard = {}\ndepth = {}, strict = {}, lenient = {}".format( np.round (sum_mixture, 2), kwargs["MAKEONE_STANDARD"][0],  kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize=12)


    # 기존에 그렸던 histogram

    for j in sorted(list(set(bunch.membership))):
        if j == bunch.fp_index:  # FP를 따로 그려줌
            sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, j], "rwidth": 0.8}, color=Gr_10[10],
                         kde=False, bins=50, label="cluster FP (mixture = {})".format(  str(round((mixture[0][j]), 2))), ax = ax)

        elif j in bunch.makeone_index:
            if j not in bunch.lowvafcluster_index [0]:
                sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, j], "linewidth": 2, "edgecolor": "red"}, kde_kws={
                            "linewidth": 5, "color": "gray"}, color=colorlist [j], kde=False, bins=50, label="cluster {} ({}, cdf = {}, n = {})".format(  j,  round((mixture[0][j]), 2), bunch.cluster_cdf_list[0][j], np.bincount(bunch.membership)[j] ), ax = ax)
            else:
                sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, j], "linewidth": 1.4, "edgecolor": "black"}, kde_kws={
                            "linewidth": 5, "color": "gray"}, color=colorlist [j], kde=False, bins=50, label="(lowvaf)cluster {} ({}, cdf = {}, n = {})".format(  j,  round((mixture[0][j]), 2), bunch.cluster_cdf_list[0][j], np.bincount(bunch.membership)[j] ), ax = ax)

        else:
            sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, j], "rwidth": 0.8}, color=colorlist[j],
                         kde=False, bins=50, label="cluster {}  (mixture = {})".format(j,  str(round((mixture[0][j]), 2))), ax = ax )
    ax.set_ylabel ("count (weighted)", fontdict={"fontsize": 14} )


    # for j in sorted(list(set(bunch.membership))):
    #     np_vaf_new_index, np_vaf_new = extract.npvaf(np_vaf, bunch.membership, j)           # extract membership1  == clone_num
    #     mean_x = np.mean (np_vaf_new [:, 0] * 2)

    #     vaf_list = []
    #     for k in range (0, len(bunch.membership)):
    #         if random.random() < membership_p_normalize[k, j]:
    #             vaf_list.append ( np_vaf[k, 0] * 2 )
        
    #     if (len ( set(vaf_list) ) <= 1 ):
    #         continue

    #     x = np.linspace(0, 1.5, 601)
    #     kde_np_vaf_new = kde.gaussian_kde( vaf_list )                
    #     weight = len(np_vaf_new) / kwargs["NUM_MUTATION"]  
    #     y = kde_np_vaf_new(x) * weight
    #     x_max, y_max = x [np.argmax(y)], y [np.argmax(y)]

    #     if j == bunch.fp_index:  # FP를 따로 그려줌
    #         ax.plot(x, y, color = Gr_10[10], linewidth=5, label="FP ({}, n = {})".format( np.round( bunch.mixture[0, j], 2), np.bincount(bunch.membership)[j]))
    #         ax.text(0, np.max(y) + 0.1, "FP({})".format (mean_x), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})
    #         #ax.axvline( x = 0,  ymin = 0, ymax =  y_max / (ax.get_ylim()[1] - ax.get_ylim()[0]) , linestyle = '--', linewidth = 2, color = Gr_10[10] )
    #     elif j in bunch.makeone_index:
    #         ax.plot(x, y, color = colorlist[j], linewidth=5, label="clone {} ({}, cdf = {}, n = {})".format(j, np.round( bunch.mixture[0, j], 2), bunch.cluster_cdf_list[0][j], np.bincount(bunch.membership)[j]) )
    #         if j not in bunch.lowvafcluster_index [0]:
    #             ax.fill_between(x, y, color = colorlist[j], alpha = 0.5, label = None )
    #         ax.text(bunch.mixture[0, j], np.max(y) + 0.1, "{}\n(n={})".format( np.round ( bunch.mixture[0,  j], 2),  np.bincount(bunch.membership)[j]), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})
    #         #ax.axvline( x = bunch.mixture[0, j],  ymin = 0, ymax = y_max / (ax.get_ylim()[1] - ax.get_ylim()[0]) , linestyle = '-', linewidth = 2, color = colorlist[j] )
    #     else:
    #         ax.plot(x, y, color=colorlist[j], linewidth=5, label="clone {} ({}, cdf = {}, n = {})".format(j, np.round( bunch.mixture[0, j], 2), bunch.cluster_cdf_list[0][j], np.bincount(bunch.membership)[j]), linestyle="-.")
    #         ax.text(bunch.mixture[0, j], np.max(y) + 0.1, "{}\n(n={})".format( np.round( bunch.mixture[0, j], 2) ,  np.bincount(bunch.membership)[j]), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})
    #         #ax.axvline( x = bunch.mixture[0, j],  ymin = 0, ymax = y_max / (ax.get_ylim()[1] - ax.get_ylim()[0]) , linestyle = '--', linewidth = 2, color = colorlist[j] )
    # ax.set_ylabel ("Density", fontdict={"fontsize": 14} )



    ax.set_xticks ( np.round (np.arange (0, 1.01, 0.2), 1) )
    ax.set_xlim( 0, 1.01 )
    ax.set_ylim( 0, ax.get_ylim()[1] * 1.3 )
    ax.set_yticks ([])

    xtlabels = list ( ax.get_xticks () )
    ax.set_xticklabels ( [ str(i) + "\n({})".format( int ( int(kwargs["MEAN_DEPTH"]) * float (i) * 0.5)  ) for i in xtlabels ])
    

    ax.set_xlabel("VAF x 2\n(Alt)", fontdict={"fontsize": 14})
    ax.legend(loc='best')

    matplotlib.pyplot.show()
    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)



################################################################################################################################################################################################



def drawfigure_2d(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}
    sum_mixture_beforenormalization = np.round ( sum_mixture_beforenormalization , 3)

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    fig.subplots_adjust ( bottom = 0.15, top = 0.85)
    ax.set_xlim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_ylim(0, np.max(np_vaf[:, :]) * 2.1)
    # matplotlib.pyplot.xlabel("VAF x 2 of the Sample 1", fontdict={"fontsize": 14})
    # matplotlib.pyplot.ylabel("VAF x 2 of the Sample 2", fontdict={"fontsize": 14})
    
    xtlabels = list ( ax.get_xticks () )
    ax.set_xticklabels ( [ str( round(i, 1) ) + "\n({})".format( int ( int(kwargs["MEAN_DEPTH"]) * float (i) * 0.5)  ) for i in xtlabels ])
    ytlabels = list ( ax.get_yticks () )
    ax.set_yticklabels ( [ str( round(i, 1) ) + "\n({})".format( int ( int(kwargs["MEAN_DEPTH"]) * float (i) * 0.5)  ) for i in ytlabels ])

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
        # Set the layout to tight to ensure the rectangle is not cut off
        #matplotlib.pyplot.tight_layout()        
        matplotlib.pyplot.title("postnorm = {}, prenorm = {}, stan = {}\ndepth = {}, strict = {}, lenient = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ), sum_mixture_beforenormalization, kwargs["MAKEONE_STANDARD"][1], kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize = 12)
        
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (raw) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)        
        matplotlib.pyplot.title("prenorm = {}, standard = {}\ndepth = {}, strict = {}, lenient = {}".format(   sum_mixture_beforenormalization, kwargs["MAKEONE_STANDARD"][1], kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize = 12)


    matplotlib.pyplot.style.use("seaborn-white")
   
    for k in range (len ( bunch.membership)):
        if (bunch.fp_index == bunch.membership[k]):
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = Gr_10[10]  )
        else:
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = colorlist[ samplename_dict [ bunch.membership[k] ]]  )


    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        # Drawe the centroid based on mixture
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        if x_mean < 0.02:
            ax.text(x_mean, y_mean, "{0}".format( [y_mean]), verticalalignment='top', horizontalalignment='right', fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        elif y_mean < 0.02:
            ax.text(x_mean, -0.02, "{0}".format( [x_mean]), verticalalignment='top', horizontalalignment='center',  fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        else:
            ax.text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )


        if (bunch.fp_index == sample_index):
            ax.scatter(x_mean, y_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.lowvafcluster_index [1]:
                if sample_index in bunch.makeone_index:
                    ax.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=400, label="(lowvaf)cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
                else:
                    ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=200, label="(lowvaf)cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                if sample_index in bunch.makeone_index:
                    ax.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=400, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
                else:
                    ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        else:
            if sample_index in bunch.lowvafcluster_index [1]:
                ax.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="(lowvaf)cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                ax.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))


    ax.legend(loc='upper right')
    matplotlib.pyplot.show()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)



################################################################################################################################################################################################

def drawfigure_2d_soft(bunch, np_vaf, output_filename, sum_mixture, **kwargs):
    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-white")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-white") 
    elif kwargs["OPTION"] in ["fp", "fp"]:
        matplotlib.pyplot.style.use("Solarize_Light2")

    NUM_MUTATION = len(bunch.membership)
    NUM_CLONE = kwargs["NUM_CLONE"]
    NUM_BLOCK = kwargs["NUM_BLOCK"]

    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    # font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    # font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    # for font in font_dirs:
    #     matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]
    matplotlib.rcParams["text.color"] = "blue"
    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    fig.subplots_adjust ( bottom = 0.15, top = 0.85)

    if bunch.includefp == True:
        fp_color_num = samplename_dict[bunch.fp_index]
        colorlist[fp_color_num] = Gr_10[10]

    dimensionreduction = ""
    if bunch.mixture.shape[0] > 2:  # More than 3 samples
        dimensionreduction = "SVD"

    if dimensionreduction == "SVD":
        print("SVD → 2D")
        tsvd = TruncatedSVD(n_components=2)
        tsvd.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = tsvd.transform(np.concatenate(  [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = tsvd.transform(np.concatenate( [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        ax.set_xlim ( np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0])  * 2.1)
        ax.set_ylim ( np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1 )
        ax.set_xlabel("SVD1")
        ax.set_ylabel("SVD2")
    elif dimensionreduction == "PCA":
        print("PCA → 2D")
        pca = PCA(n_components=2)
        pca.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = pca.transform(np.concatenate( [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = pca.transform(np.concatenate(  [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        ax.set_xlim ( np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0])  * 2.1)
        ax.set_ylim ( np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1 )
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
    else:
        ax.set_xlim ( 0,  np.max(np_vaf[:, :]) * 2.1 )
        ax.set_ylim ( 0,  np.max(np_vaf[:, :]) * 2.1 )
        ax.set_xlabel("Feature 1 : VAF x 2 of the Sample 1")
        ax.set_ylabel("Feature 2 : VAF x 2 of the Sample 2")


    if (kwargs["STEP"] <= ( kwargs["SOFT_LENIENT"] )) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (lenient) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (raw) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)    


    if kwargs["STEP"] <= (kwargs["SOFT_LENIENT"] ):  # normalization 한 것
        ax.set_title("prenorm = {}, standard = {}\ndepth = {}, strict = {}, lenient = {}".format(np.round (sum_mixture, 2) , kwargs["MAKEONE_STANDARD"][1],  kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize=12)
    else:
        ax.set_title("prenorm = {}, standard = {}\ndepth = {}, strict = {}, lenient = {}".format( np.round (sum_mixture, 2), kwargs["MAKEONE_STANDARD"][1],  kwargs ["MEAN_DEPTH"], bunch.checkall_strict, bunch.checkall_lenient ), fontsize=12)


    if bunch.includefp == "Yes":
        fp_color_num = samplename_dict[bunch.fp_index]
        colorlist[fp_color_num] = Gr_10[10]

    
    
    for j in range(0, NUM_CLONE):
        if j == bunch.fp_index:
            for k in bunch.fp_member_index:
                ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=Gr_10[10])
        else:
            for k in range(NUM_MUTATION):
                if bunch.membership_p_normalize[k, j] > 0.8:
                    ax.scatter( np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[samplename_dict[j]]])
                elif bunch.membership_p_normalize[k, j] > 0.1:
                    ax.scatter( np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=bunch.membership_p_normalize[k, j], color=[colorlist[samplename_dict[j]]])
                else: 
                    continue

    for sample_index, sample_key in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)

        if x_mean < 0.02:
            ax.text(x_mean, y_mean, "{0}".format( [y_mean]), verticalalignment='top', horizontalalignment='right', fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        elif y_mean < 0.02:
            ax.text(x_mean, -0.02, "{0}".format( [x_mean]), verticalalignment='top', horizontalalignment='center',  fontdict = {"fontsize": 14, "fontweight" : "bold"} )    
        else:
            ax.text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )

        if (bunch.fp_index == sample_index):
            ax.scatter(x_mean, y_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.lowvafcluster_index [1]:
                if sample_index in bunch.makeone_index:
                    ax.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=200, label="(lowvaf)cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
                else:
                    ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100, label="(lowvaf)cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                if sample_index in bunch.makeone_index:
                    ax.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
                else:
                    ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            if sample_index in bunch.lowvafcluster_index [1]:
                ax.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="(lowvaf)cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                ax.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

    xtlabels = list ( ax.get_xticks () )
    ax.set_xticklabels ( [ str( round(i, 1) ) + "\n({})".format( int ( int(kwargs["MEAN_DEPTH"]) * float (i) * 0.5)  ) for i in xtlabels ])
    ytlabels = list ( ax.get_yticks () )
    ax.set_yticklabels ( [ str( round(i, 1) ) + "\n({})".format( int ( int(kwargs["MEAN_DEPTH"]) * float (i) * 0.5)  ) for i in ytlabels ])

    ax.legend(loc='upper right')
    matplotlib.pyplot.show()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)


################################################################################################################################################################################################

def drawfigure_3d(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    from mpl_toolkits.mplot3d import Axes3D
    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    #fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    fig = matplotlib.pyplot.figure( figsize = (6, 6) )
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_ylim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_zlim(0, np.max(np_vaf[:, :]) * 2.1)


    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    else:
        matplotlib.pyplot.suptitle("clone{}.{}-{} (raw) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}, prenorm = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ), sum_mixture_beforenormalization), fontsize = 12)


    # if (bunch.makeone_index == []) | (bunch.makeone_index == None) | (bunch.checkall_strict == False):  
    #     matplotlib.pyplot.style.use("Solarize_Light2")
    # else:
    matplotlib.pyplot.style.use("seaborn-white")

    if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
        import matplotlib.patches as patches
        fig_width, fig_height = fig.get_size_inches()
        rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
        fig.patches.append(rect)
        # Set the layout to tight to ensure the rectangle is not cut off

    
    #matplotlib.pyplot.scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, color=[colorlist[samplename_dict[k]] for k in bunch.membership])
    for k in range (len ( bunch.membership)):
        if (bunch.fp_index == bunch.membership[k]):
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, np_vaf[k, 2] * 2, color = Gr_10[10]  )
        else:
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, np_vaf[k, 2] * 2, color = colorlist[ samplename_dict [ bunch.membership[k] ]]  )


    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        # Drawe the centroid based on mixture
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        z_mean = round(bunch.mixture[2][sample_index], 2)

        ax.text(x_mean, y_mean, z_mean, "{0}".format( [x_mean, y_mean, z_mean]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )

        if (bunch.fp_index == sample_index):
            ax.scatter(x_mean, y_mean, z_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                ax.scatter(x_mean, y_mean, z_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=400, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                ax.scatter(x_mean, y_mean, z_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            ax.scatter(x_mean, y_mean, z_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))


    ax.legend(loc='upper right')
    matplotlib.pyplot.show()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)





def drawfigure_3d_SVD(bunch, np_vaf, output_filename, sum_mixture_beforenormalization, **kwargs):
    from mpl_toolkits.mplot3d import Axes3D
    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    fig, ax = matplotlib.pyplot.subplots(ncols=1, nrows = 1, figsize=(6, 6))
    fig.subplots_adjust ( bottom = 0.15, top = 0.85)
    ax.set_xlim(0, np.max(np_vaf[:, :]) * 2.1)
    ax.set_ylim(0, np.max(np_vaf[:, :]) * 2.1)

    tsvd = TruncatedSVD(n_components=2)
    tsvd.fit( np_vaf )

    np_vaf, mixture = tsvd.transform( np.concatenate(  [np_vaf, bunch.mixture.T]) )[:-kwargs["NUM_CLONE"]], tsvd.transform(np.concatenate( [np_vaf, bunch.mixture.T]) )[-kwargs["NUM_CLONE"]:].T
    #np_vaf = tsvd.transform( np.concatenate(  [np_vaf, bunch.mixture.T]) )
    matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
    matplotlib.pyplot.style.use("seaborn-white")


    if bunch.hard_or_soft  == "hard":
        if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
            matplotlib.pyplot.suptitle("clone{}.{}-{} (postnormalization) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
        else:
            matplotlib.pyplot.suptitle("clone{}.{}-{} (raw) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
        matplotlib.pyplot.title("sum_child = {}, prenorm = {}\nmakewhat = {}, strict = {}, lenient = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ), sum_mixture_beforenormalization, kwargs ["MAKEWHAT"], bunch.checkall_strict , bunch.checkall_lenient ), fontsize = 12)

        if (kwargs["STEP"] <= (kwargs["COMPULSORY_NORMALIZATION"] - 1)) :
            import matplotlib.patches as patches
            fig_width, fig_height = fig.get_size_inches()
            rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
            fig.patches.append(rect)
            # Set the layout to tight to ensure the rectangle is not cut off

    elif bunch.hard_or_soft  == "soft":
        if (kwargs["STEP"] <= (kwargs["SOFT_LENIENT"] - 1)) :
            matplotlib.pyplot.suptitle("clone{}.{}-{} (lenient) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
        else:
            matplotlib.pyplot.suptitle("clone{}.{}-{} (raw) : {}".format(kwargs["NUM_CLONE_ITER"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
        matplotlib.pyplot.title("sum_child = {}, prenorm = {}\nmakewhat = {}, strict = {}, lenient = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ), sum_mixture_beforenormalization, kwargs ["MAKEWHAT"], bunch.checkall_strict , bunch.checkall_lenient ), fontsize = 12)


        if (kwargs["STEP"] <= (kwargs["SOFT_LENIENT"] - 1)) :
            import matplotlib.patches as patches
            fig_width, fig_height = fig.get_size_inches()
            rect = patches.Rectangle((0, 0), fig_width * fig.dpi, fig_height * fig.dpi, linewidth=8, edgecolor="#C23373", facecolor='none')
            fig.patches.append(rect)
            # Set the layout to tight to ensure the rectangle is not cut off
            matplotlib.pyplot.tight_layout()


    
    for k in range (len ( bunch.membership)):
        if (bunch.fp_index == bunch.membership[k]):
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = Gr_10[10]  )
        else:
            ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, color = colorlist[ samplename_dict [ bunch.membership[k] ]]  )


    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= mixture.shape[1]):
            continue

        # Drawe the centroid based on mixture
        x_mean = round(mixture[0][sample_index], 2)
        y_mean = round(mixture[1][sample_index], 2)

        ax.text(x_mean, y_mean, "{0}".format( [ round(bunch.mixture[0][sample_index],2), round(bunch.mixture[1][sample_index],2),  round(bunch.mixture[2][sample_index],2)  ]), verticalalignment='top', horizontalalignment='center', fontdict = {"fontsize": 16, "fontweight" : "bold"} )

        if (bunch.fp_index == sample_index):
            ax.scatter(x_mean, y_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))            
        elif (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                ax.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=400, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            ax.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

 
    ax.legend(loc='upper right')
    matplotlib.pyplot.show()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)
