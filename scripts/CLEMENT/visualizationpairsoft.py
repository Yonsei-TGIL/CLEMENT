import palettable
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import kde
import extract
from sklearn.decomposition import TruncatedSVD, PCA
import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)



def drawfigure_2d(membership1, sample_dict_rev, membership2, membership2_outside, membership2_p_normalize, membership2_p_normalize_new, mixture2, 
                                axis_index, df_inside_index, df_outside_index, output_filename, np_vaf,  includeoutlier,  dimensionreduction="None"):
    NUM_MUTATION = membership2_p_normalize.shape[0]
    NUM_CLONE = membership2_p_normalize.shape[1]
    NUM_BLOCK = np_vaf.shape[1]

    vivid_10 = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
    bdo = palettable.lightbartlein.diverging.BlueDarkOrange18_18.mpl_colors
    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    if includeoutlier == "Yes":
        colorlist[np.max(membership2)] = Gr_10[16]        # Outlier (= FP) : Black

    maxmaxmax_NUM_CLONE = np.max(membership2) + 1
    fig, ax = plt.subplots(ncols=2, figsize=(15, 6))

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
    else:
        ax[0].set_title("ANSWER FIGURE OF MRS (n = {0})".format(len(membership1)))
        ax[1].set_title("MYANSWER_NUM_CLONE : {0} (n = {1})".format(maxmaxmax_NUM_CLONE, len(membership2)))
        ax[0].axis([0,  np.max(np_vaf[:, :]) * 2.1,   0,  np.max(np_vaf[:, :]) * 2.1])
        ax[1].axis([0,  np.max(np_vaf[:, :]) * 2.1,   0,  np.max(np_vaf[:, :]) * 2.1])
        ax[0].set_xlabel("Feature 1 : VAF x 2 of Block 1")
        ax[1].set_xlabel("Feature 1 : VAF x 2 of Block 1")
        ax[0].set_ylabel("Feature 2 : VAF x 2 of Block 2")
        ax[1].set_ylabel("Feature 2 : VAF x 2 of Block 2")

    # Left scatter : just draw it
    ax[0].scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, alpha=0.6,   color=[colorlist[k] for k in membership1])


    # Right figure: draw child clone
    for j in range(membership2_p_normalize_new.shape[1] - 1):
        for k in range (NUM_MUTATION):
            if k not in df_inside_index:
                if membership2_p_normalize_new[k,j] > 0.8:
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color= colorlist[j] )
                elif membership2_p_normalize_new[k,j] > 0.1:
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=membership2_p_normalize[k,j], color= colorlist[j] )
                else:        # If contribution is less than 10%, than pass it
                    continue
    # Right figure: draw axis clone as black
    for j in axis_index:
        for k in range (NUM_MUTATION):
            if membership2_p_normalize[k,j] > 0.8:
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6, color= Gr_10[16] )
            elif membership2_p_normalize[k,j] > 0.1:
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=membership2_p_normalize[k,j], color= Gr_10[16] )
            else:        # If contribution is less than 10%, than pass it
                continue
    # Right figure: draw parent clone 
    if includeoutlier == "No":
        for k in range (NUM_MUTATION):
            if k in df_outside_index:
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=[colorlist[ membership2[k] ]] )
    else:
        for k in range (NUM_MUTATION):
            if k in df_outside_index:
                if membership2[k] != np.max(membership2):       # If not outlier
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=[colorlist[ membership2[k] ]] )
                else:
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=Gr_10[16] )
    # Draw inside outlier black
    for k in range (NUM_MUTATION):
        if k in df_inside_index:
            if membership2[k] == np.max(membership2):
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=Gr_10[16] )



    if (dimensionreduction != "SVD") & ( dimensionreduction != "PCA" ):
        #  Left figure : Draw answer (square)
        for sample_index in range(int(np.max(membership1)) + 1):
            # Based on membership & np_vaf information
            x_mean = round(np.mean(np_vaf[[x for x in range( len(membership1)) if membership1[x] == sample_index]][:, 0] * 2), 2)
            y_mean = round(np.mean(np_vaf[[x for x in range( len(membership1)) if membership1[x] == sample_index]][:, 1] * 2), 2)

            ax[0].text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top')
            ax[0].scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100,
                        label=str(sample_dict_rev[sample_index]) + " : " + str(list(membership1).count(sample_index)))
            ax[0].legend()

        #  Right figure : Draw my EM decomposition (square)
        xx = maxmaxmax_NUM_CLONE - 1 if "outlier" in output_filename else maxmaxmax_NUM_CLONE

        for sample_index in range(xx):   # Do not paint square
            # Based on mixture information
            x_mean = mixture2[0][sample_index]
            y_mean = mixture2[1][sample_index]
            ax[1].text(x_mean, y_mean, "{0}".format([x_mean, y_mean]), verticalalignment='top')
            ax[1].scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100,
                        label="cluster" + str(sample_index) + " : " + str(list(membership2).count(sample_index)))

            ax[1].legend()

    if output_filename != "NotSave":
        plt.savefig(output_filename)
    plt.show()