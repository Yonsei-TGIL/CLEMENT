import palettable
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import kde
import extract
from sklearn.decomposition import TruncatedSVD, PCA
import warnings
warnings.simplefilter (action = 'ignore', category = FutureWarning)



# Block 1, 2:  data point를 2차원 평면상에 그려보기
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
        colorlist[np.max(membership2)] = Gr_10[16]        # Outlier는 까만색으로 지정해준다

    # if includeoutlier == "Yes":
    #     outlier_color_num = samplename_dict[np.max(list(samplename_dict.keys()))]       # 맨 마지막 번호의 색깔번호 (Outlier 번호)
    #     colorlist [samplename_dict[np.max(list(samplename_dict.keys()))]] = Gr_10[4]

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
        # ax[0].set_xlabel("Feature 1 : VAF x 2 of Block 1")
        # ax[1].set_xlabel("Feature 1 : VAF x 2 of Block 1")
        # ax[0].set_ylabel("Feature 2 : VAF x 2 of Block 2")
        # ax[1].set_ylabel("Feature 2 : VAF x 2 of Block 2")

    # 왼쪽 scatter는 그냥 그리면 된다
    ax[0].scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, alpha=0.6,   color=[colorlist[k] for k in membership1])
    #ax[1].scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, alpha=0.6,   color=[colorlist[k] for k in membership2])

    # 오른쪽 scatter가 그리기 어렵다. 일단 child clone을 맨 먼저 그리자
    for j in range(membership2_p_normalize_new.shape[1] - 1):  #  Child clone만 그리기
        for k in range (NUM_MUTATION):
            if k not in df_inside_index:
                if membership2_p_normalize_new[k,j] > 0.8:
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color= colorlist[j] )
                elif membership2_p_normalize_new[k,j] > 0.1:
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=membership2_p_normalize[k,j], color= colorlist[j] )
                else:        # 10%도 기여하지 못한다면 칠하지 말고 넘어가자
                    continue
    # 오른쪽 그림의 axis clone을 까맣게 soft visualization하기
    for j in axis_index:
        for k in range (NUM_MUTATION):
            if membership2_p_normalize[k,j] > 0.8:
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6, color= Gr_10[16] )
            elif membership2_p_normalize[k,j] > 0.1:
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=membership2_p_normalize[k,j], color= Gr_10[16] )
            else:        # 10%도 기여하지 못한다면 칠하지 말고 넘어가자
                continue
    # parent 를 그린다.  이때 outlier의 존재 유무에 따라 조심해야한다
    if includeoutlier == "No":
        for k in range (NUM_MUTATION):
            if k in df_outside_index:
                ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=[colorlist[ membership2[k] ]] )
    else:
        for k in range (NUM_MUTATION):
            if k in df_outside_index:
                if membership2[k] != np.max(membership2):       # Outlier가 아니면
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=[colorlist[ membership2[k] ]] )
                else:
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=Gr_10[16] )
    # insider outlier를  까맣게 그린다
    for k in range (NUM_MUTATION):
        if k in df_inside_index:
            if membership2[k] == np.max(membership2):
                    ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=0.6,   color=Gr_10[16] )



    if (dimensionreduction != "SVD") & ( dimensionreduction != "PCA" ):
        #  왼쪽 그림 :  정답을 그리기 (square )
        for sample_index in range(int(np.max(membership1)) + 1):

            # membership & np_vaf 정보를 바탕으로
            x_mean = round(np.mean(np_vaf[[x for x in range( len(membership1)) if membership1[x] == sample_index]][:, 0] * 2), 2)
            y_mean = round(np.mean(np_vaf[[x for x in range( len(membership1)) if membership1[x] == sample_index]][:, 1] * 2), 2)

            ax[0].text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top')
            ax[0].scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100,
                        label=str(sample_dict_rev[sample_index]) + " : " + str(list(membership1).count(sample_index)))
            ax[0].legend()

        #  오른쪽 그림 : 내가 구한 EM 을 그리기 (square)
        xx = maxmaxmax_NUM_CLONE - 1 if "outlier" in output_filename else maxmaxmax_NUM_CLONE

        for sample_index in range(xx):   # square 그려주지 말자
            # mixture 정보를 바탕으로
            x_mean = mixture2[0][sample_index]
            y_mean = mixture2[1][sample_index]
            ax[1].text(x_mean, y_mean, "{0}".format([x_mean, y_mean]), verticalalignment='top')
            ax[1].scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100,
                        label="cluster" + str(sample_index) + " : " + str(list(membership2).count(sample_index)))

            ax[1].legend()

    if output_filename != "NotSave":
        plt.savefig(output_filename)
    plt.show()