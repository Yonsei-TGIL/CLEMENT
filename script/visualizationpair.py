import palettable
import matplotlib
import seaborn as sns
import numpy as np
from scipy.stats import kde
import extract
from sklearn.decomposition import TruncatedSVD, PCA


def drawfigure_1d(membership1, sample_dict_rev, membership2, mixture2, output_filename, np_vaf, **kwargs):
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    bdo = palettable.lightbartlein.diverging.BlueDarkOrange18_18.mpl_colors
    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    # font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    # font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    # for font in font_dirs:
    #     matplotlib.font_manager.fontManager.addfont(font)
    
    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]


    maxmaxmax_NUM_CLONE = np.max(membership2) + 1
    fig, ax = matplotlib.pyplot.subplots(ncols=2, figsize=(15, 6))

    ax[0].set_title("ANSWER FIGURE OF MRS (n = {0})".format(len(membership1)), fontdict = {"fontsize" : 20})
    ax[1].set_title("MYANSWER_NUM_CLONE : {0} (n = {1})".format(maxmaxmax_NUM_CLONE, len(membership2)) , fontdict = {"fontsize" : 20})

    max_y = 0

    x = np.linspace(0, 2, 200)

    for k in sorted(list(set(membership1))): 
        np_vaf_new_index, np_vaf_new = extract.npvaf(np_vaf, membership1, k) 
        kde_np_vaf_new = kde.gaussian_kde(np_vaf_new[:, 0] * 2)
        
        weight = len(np_vaf_new) / len(np_vaf) 
        y = kde_np_vaf_new(x) * weight
        if max_y < np.max(y):
            max_y = np.max(y)

        ax[0].plot(x, y, color=colorlist[k], label=sample_dict_rev[k])
        ax[0].text(np.argmax(y) / 100, np.max(y) * 1.2,
                   "{0}".format(np.argmax(y) / 100), verticalalignment='top')

        np_vaf_new_index, np_vaf_new = extract.npvaf(np_vaf, membership2, k)
        kde_np_vaf_new = kde.gaussian_kde(np_vaf_new[:, 0] * 2)
        y = kde_np_vaf_new(x)
        if max_y < np.max(y):
            max_y = np.max(y)

        ax[1].plot(x, y, color=colorlist[k], label="cluster {0}".format(k))
        ax[1].text(np.argmax(y) / 100, np.max(y) * 1.2, "{0}".format(np.argmax(y) / 100), verticalalignment='top', fontdict = {"fontsize": 16, "fontweight" : "bold"})

    ax[0].axis([0,  np.max(np_vaf[:, :]) * 2.1,  0,  max_y * 1.3])
    ax[1].axis([0,  np.max(np_vaf[:, :]) * 2.1,  0,  max_y * 1.3])

    ax[0].legend()
    ax[1].legend()

    ax[0].set_xlabel("Mixture ( = VAF x 2)")
    ax[1].set_xlabel("Mixture ( = VAF x 2)")
    ax[0].set_ylabel("Density")
    ax[1].set_ylabel("Density")

    matplotlib.pyplot.savefig(output_filename)



def drawfigure_2d(membership_left, mixture_left, membership_right, mixture_right, score_df, output_filename, fig1title, fig2title, np_vaf, includeoutlier,  makeone_index, dimensionreduction="None", **kwargs):
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    bdo = palettable.lightbartlein.diverging.BlueDarkOrange18_18.mpl_colors
    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]


    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    #print (matplotlib.font_manager.FontProperties(fname = font).get_name())

    matplotlib.pyplot.style.use("seaborn-white")
    matplotlib.rcParams["font.family"] = kwargs["FONT_FAMILY"]

    if (includeoutlier == True ) &  ("FP" in list(score_df["answer"])):
        colorlist[ np.where(score_df["answer"] == "FP")[0][0] ] = Gr_10[9]        # Draw FP in black

    if mixture_right.shape[0] > 2:  # more than 3 samples
        dimensionreduction = "SVD"

    maxmaxmax_NUM_CLONE = np.max(membership_right) + 1
    fig, ax = matplotlib.pyplot.subplots(ncols=2, figsize=(15, 6))

    if dimensionreduction == "SVD":
        print("SVD → 2D")
        tsvd = TruncatedSVD(n_components=2)
        tsvd.fit(np_vaf)
        np_vaf = tsvd.transform(np_vaf)
        ax[0].axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) *  2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        ax[1].axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) *  2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        ax[0].set_xlabel("SVD1")
        ax[1].set_xlabel("SVD1")
        ax[0].set_ylabel("SVD2")
        ax[1].set_ylabel("SVD2")
        ax[0].set_title("{}".format(fig1title) ,  fontdict = {"fontsize" : 20} )
        ax[1].set_title("{}".format(fig2title) , fontdict = {"fontsize" : 20} )
        
    elif dimensionreduction == "PCA":
        print("PCA → 2D")
        pca = PCA(n_components=2)
        pca.fit(np_vaf)
        np_vaf = pca.transform(np_vaf)
        ax[0].axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) *  2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        ax[1].axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) *  2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        ax[0].set_xlabel("PC1")
        ax[1].set_xlabel("PC1")
        ax[0].set_ylabel("PC2")
        ax[1].set_ylabel("PC2")
        ax[0].set_title("{}".format(fig1title) ,  fontdict = {"fontsize" : 20} )
        ax[0].set_title("\n\nNUM_CLONE = {}".format( mixture_left.shape[1]) , loc = 'right', fontdict = {"fontsize" : 12} )
        ax[1].set_title("{}".format(fig2title) , fontdict = {"fontsize" : 20} )
        ax[1].set_title("\n\nNUM_CLONE = {}".format( mixture_right.shape[1]) , loc = 'right', fontdict = {"fontsize" : 12} )
        
    else:
        ax[0].set_title("{}".format(fig1title) ,  fontdict = {"fontsize" : 20} )
        ax[1].set_title("{}".format(fig2title) , fontdict = {"fontsize" : 20} )
        ax[0].axis([0,  np.max(np_vaf[:, :]) * 2.1,   0,  np.max(np_vaf[:, :]) * 2.1])
        ax[1].axis([0,  np.max(np_vaf[:, :]) * 2.1,   0,  np.max(np_vaf[:, :]) * 2.1])
        # ax[0].set_xlabel("Mixture ( = VAF x 2) of Sample 1")
        # ax[1].set_xlabel("Mixture ( = VAF x 2) of Sample 1")
        # ax[0].set_ylabel("Mixture ( = VAF x 2) of Sample 2")
        # ax[1].set_ylabel("Mixture ( = VAF x 2) of Sample 2")

    

    # Scatter plot
    cluster_unmatched = {}
    color_unmatched = score_df.shape[0]   # 이것부터 시작
    for k in range(len(membership_left)):
        try:
            i =  np.where(score_df["answer"] == membership_left[k])[0][0] 
            ax[0].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[ i ]])
        except:
            i = color_unmatched
            ax[0].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[ i ]])
            color_unmatched += 1
            print ("What happened in visualizationpair.py?\n{}\n{}\n{}".format ( membership_left[k], score_df, i ) )
    for k in range(len(membership_right)):
        try:
            i = np.where(score_df["predicted"] == membership_right[k])[0][0] 
            if (includeoutlier == True) & (membership_right[k] == mixture_right.shape[1] - 1):   #FP
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color= Gr_10[10])
            else:    
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[ i ]])
        except:
            i = color_unmatched
            ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[ i ]])
            color_unmatched += 1
            print ("What happened in visualizationpair.py?\n{}\n{}\n{}".format ( membership_right[k], score_df, i ) )



    if (dimensionreduction != "SVD") & ( dimensionreduction != "PCA" ):
        #  Left figure : answer
        x_mean_list, y_mean_list = [], []
        for samplename_left_character in list(np.unique(membership_left)):
            x_mean = round(np.mean(np_vaf[[x for x in range( len(membership_left)) if membership_left[x] == samplename_left_character]][:, 0] * 2), 2)
            y_mean = round(np.mean(np_vaf[[x for x in range( len(membership_left)) if membership_left[x] == samplename_left_character]][:, 1] * 2), 2)
            x_mean_list.append(x_mean)
            y_mean_list.append(y_mean)

            ax[0].text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', fontdict = {"fontsize": 16, "fontweight" : "bold"})
            i = np.where(score_df["answer"] == samplename_left_character)[0][0] 
            ax[0].scatter(x_mean, y_mean, marker='s', color=colorlist[ i ], edgecolor='black', s=100,
                        label=str(samplename_left_character) + " : " + str(list(membership_left).count(samplename_left_character)))
        ax[0].legend()
        ax[0].text (np.max(np_vaf[:, :]), np.max(np_vaf[:, :]) * 1.95, "sum = [{},{}]".format( round( np.sum ( np.array(x_mean_list) ) , 2) , round( np.sum( np.array(y_mean_list) ), 2)  ) ,  ha = 'center', va = 'top', fontdict = {"fontsize" : 12} )
                    #set_title("{}\nsum = [{},{}]".format(fig1title, round( np.sum ( np.array(x_mean_list) ) , 2) , round( np.sum( np.array(y_mean_list) ), 2)  ) ,  fontdict = {"fontsize" : 20} )
        

        # Right figure : my dataset
        xx = mixture_right.shape[1] 
        x_mean_list, y_mean_list = [], []

        for sample_index in range(xx): 
            x_mean = mixture_right[0][sample_index]
            y_mean = mixture_right[1][sample_index]
            ax[1].text(x_mean, y_mean, "{0}".format([x_mean, y_mean]), verticalalignment='top', fontdict = {"fontsize": 16, "fontweight" : "bold"})
            try:
                i = np.where(score_df["predicted"] == sample_index)[0][0] 
            except:
                continue
            
            if (x_mean == 0) & (y_mean == 0):  # FP
                ax[1].scatter(x_mean, y_mean, marker='s', color= Gr_10[10], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)))            
                continue

            x_mean_list.append(x_mean)
            y_mean_list.append(y_mean)
            if (makeone_index != []) & (makeone_index != None):
                if sample_index in makeone_index:
                    ax[1].scatter(x_mean, y_mean, marker='*', color=colorlist[i], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)) )
                else:
                    ax[1].scatter(x_mean, y_mean, marker='s', color=colorlist[i], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)) )
            else:
                ax[1].scatter(x_mean, y_mean, marker='s', color=colorlist[i], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)) )

        ax[1].legend()
        ax[1].text (np.max(np_vaf[:, :]), np.max(np_vaf[:, :]) * 1.95, "sum = [{},{}]".format(round( np.sum ( np.array(x_mean_list) ) , 2) , round( np.sum( np.array(y_mean_list) ), 2)  ) ,  ha = 'center', va = 'top', fontdict = {"fontsize" : 12} )

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)
