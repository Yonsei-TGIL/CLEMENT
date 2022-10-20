import palettable
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import kde
import pandas as pd
import scipy
import extract
from sklearn.decomposition import TruncatedSVD, PCA
import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)


def drawfigure_1d(membership, mixture, membership_p_normalize, output_suptitle, output_filename, np_vaf, samplename_dict, includeoutlier, fp_index):
    NUM_MUTATION = membership_p_normalize.shape[0]
    NUM_CLONE = membership_p_normalize.shape[1]
    NUM_BLOCK = np_vaf.shape[1]
    
    
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    bdo = palettable.lightbartlein.diverging.BlueDarkOrange18_18.mpl_colors
    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    if includeoutlier == True:
        colorlist[ fp_index ] = Gr_10[8]        # Outlier (= FP) : Black
    
    plt.figure(figsize=(6, 6))
    plt.xlabel("Mixture ( = VAF x 2)", fontdict = {"fontsize" : 14})
    plt.ylabel("Density", fontdict = {"fontsize" : 14})

    plt.xlim([0,  1])

    maxmaxmax_NUM_CLONE = np.max(membership) + 1
    plt.suptitle(output_suptitle, fontsize = 20)
    plt.style.use("seaborn-white")
     
    for sample_index, sample_key in enumerate(samplename_dict):
        sample_value = samplename_dict[sample_key]
        j = sample_index
        
        # Pass FP
        if includeoutlier == True:
            if sample_index == fp_index: 
                break
        
        sns.distplot(pd.DataFrame(np_vaf[:,0] * 2, columns = ["vaf"])["vaf"], hist_kws = {"weights" : membership_p_normalize[:,j]}, 
                            color = colorlist[sample_value], kde = False, bins=50, label = "soft : cluster" + str(sample_key) + " (mean weighted vaf = " + str(round((mixture[0][j]) / 2, 2)) + ")")
    plt. legend()


    if output_filename != "NotSave":
        plt.savefig(output_filename)
    plt.show()



def drawfigure_2d(membership, mixture, membership_p_normalize,  output_suptitle, output_filename, np_vaf, samplename_dict, includeoutlier , fp_index, dimensionreduction="None"):
    import matplotlib

    NUM_MUTATION = len(membership)
    NUM_CLONE = mixture.shape[1]
    NUM_BLOCK = np_vaf.shape[1]
    
    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    bdo = palettable.lightbartlein.diverging.BlueDarkOrange18_18.mpl_colors
    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    
    fig, ax = matplotlib.pyplot.subplots (figsize = (6, 6))

    if mixture.shape[0] > 2:  # If 3D or more
        dimensionreduction = "SVD"

    if dimensionreduction == "SVD":
        print("SVD → 2D")
        tsvd = TruncatedSVD(n_components=2)
        tsvd.fit(np.concatenate([np_vaf, mixture]))               # Dimension's of mixture should be reduced also
        np_vaf = tsvd.transform(np.concatenate([np_vaf, mixture]))[:-NUM_BLOCK]
        mixture = tsvd.transform(np.concatenate([np_vaf, mixture]))[-NUM_BLOCK:]
        ax.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        ax.set_xlabel("SVD1")
        ax.set_ylabel("SVD2")
    elif dimensionreduction == "PCA":
        print("PCA → 2D")
        pca = PCA(n_components=2)
        pca.fit(np.concatenate([np_vaf, mixture]))               # Dimension's of mixture should be reduced also
        np_vaf = pca.transform(np.concatenate([np_vaf, mixture]))[:-NUM_BLOCK]
        mixture = pca.transform(np.concatenate([np_vaf, mixture]))[-NUM_BLOCK:]
        ax.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) *  2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
    else:
        ax.axis([0,  np.max(np_vaf[:, :]) * 2.1,  0,  np.max(np_vaf[:, :]) * 2.1])
        ax.set_xlabel("Mixture ( = VAF x 2) of Sample 1")
        ax.set_ylabel("Mixture ( = VAF x 2) of Sample 2")


    fig.suptitle(output_suptitle, fontsize = 20)
    plt.style.use("seaborn-white")



    if includeoutlier == True:
        outlier_color_num = samplename_dict[ fp_index ]       # fp_index (Outlier)
        colorlist [ outlier_color_num ] = Gr_10[8]


    # Paint each point 
    for j in range (0, NUM_CLONE):
        for k in range (NUM_MUTATION):
            if membership_p_normalize[k,j] > 0.8:
                ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[samplename_dict[j]]], s = 140)
            elif membership_p_normalize[k,j] > 0.1:
                ax.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=membership_p_normalize[k,j], color=[colorlist[samplename_dict[j]]], s = 170)
            else:        # If contribution is less than 10%, than pass it
                continue

    for sample_index, sample_key in enumerate(samplename_dict):
        sample_value = samplename_dict[sample_key]
        
        if sample_key not in set(membership):
            continue

        # Based on mixture information
        x_mean = round ( mixture[0][sample_index], 2) 
        y_mean = round ( mixture[1][sample_index], 2)
        #ax.text(x_mean, y_mean, "{0}".format([x_mean, y_mean]), verticalalignment='top', ha = "center")
        ax.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', linewidth = 2, s=400, alpha = 0.8, label="soft : cluster" + str(sample_index))

        fig.legend()

    if output_filename != "NotSave":
        fig.savefig(output_filename)
    matplotlib.pyplot.show()