import palettable
import matplotlib
import seaborn as sns
import numpy as np
from scipy.stats import kde
import extract
from sklearn.decomposition import TruncatedSVD, PCA


def drawfigure_1d(membership1, sample_dict_rev, membership2, mixture2, output_filename, np_vaf):
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

    #matplotlib.pyplot.style.use("seaborn-whitegrid")
    matplotlib.rcParams["font.family"] = 'arial'


    maxmaxmax_NUM_CLONE = np.max(membership2) + 1
    fig, ax = matplotlib.pyplot.subplots(ncols=2, figsize=(15, 6))

    ax[0].set_title("ANSWER FIGURE OF MRS (n = {0})".format(len(membership1)), fontdict = {"fontsize" : 20})
    ax[1].set_title("MYANSWER_NUM_CLONE : {0} (n = {1})".format(maxmaxmax_NUM_CLONE, len(membership2)) , fontdict = {"fontsize" : 20})

    max_y = 0

    x = np.linspace(0, 2, 200)

    for k in sorted(list(set(membership1))):        # 각 clone number를 따로 돈다
        np_vaf_new_index, np_vaf_new = extract.npvaf(np_vaf, membership1, k)           # membership1  == clone_num 인 것만 추려옴
        kde_np_vaf_new = kde.gaussian_kde(np_vaf_new[:, 0] * 2)
        
        weight = len(np_vaf_new) / len(np_vaf)  # weight를 줘서 total density 적분값이 1이 되도록 하자
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


# Block 1, 2:  data point를 2차원 평면상에 그려보기
def drawfigure_2d(membership_left, mixture_left, membership_right, mixture_right, score_df, output_filename, fig1title, fig2title, np_vaf, includeoutlier,  makeone_index, dimensionreduction="None"):
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
    matplotlib.rcParams["font.family"] = 'arial'

    if (includeoutlier == True ) &  ("FP" in list(score_df["answer"])):
        #colorlist[ np.where(score_df["predicted"] == np.max(membership_right))[0][0] ] = Gr_10[16]        # Outlier는 까만색으로 지정해준다
        colorlist[ np.where(score_df["answer"] == "FP")[0][0] ] = Gr_10[16]        # Outlier는 까만색으로 지정해준다

    if mixture_right.shape[0] > 2:  # 3차원 이상이면
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
        ax[0].set_title("\n\nNUM_CLONE = {}".format( mixture_left.shape[1]) , loc = 'right', fontdict = {"fontsize" : 12} )
        ax[1].set_title("{}".format(fig2title) , fontdict = {"fontsize" : 20} )
        ax[1].set_title("\n\nNUM_CLONE = {}".format( mixture_right.shape[1]) , loc = 'right', fontdict = {"fontsize" : 12} )
        
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
        ax[0].set_title("\n\nNUM_CLONE = {}".format( mixture_left.shape[1]) , loc = 'right', fontdict = {"fontsize" : 12} )
        ax[1].set_title("{}".format(fig2title) , fontdict = {"fontsize" : 20} )
        ax[1].set_title("\n\nNUM_CLONE = {}".format( mixture_right.shape[1]) , loc = 'right', fontdict = {"fontsize" : 12} )
        #matplotlib.pyplot.title(r'\fontsize{30pt}{3em}\selectfont{}{Mean WRFv3.5 LHF\r}{\fontsize{18pt}{3em}\selectfont{}(September 16 - October 30, 2012)}')

        ax[0].axis([0,  np.max(np_vaf[:, :]) * 2.1,   0,  np.max(np_vaf[:, :]) * 2.1])
        ax[1].axis([0,  np.max(np_vaf[:, :]) * 2.1,   0,  np.max(np_vaf[:, :]) * 2.1])
        ax[0].set_xlabel("Mixture ( = VAF x 2) of Sample 1")
        ax[1].set_xlabel("Mixture ( = VAF x 2) of Sample 1")
        ax[0].set_ylabel("Mixture ( = VAF x 2) of Sample 2")
        ax[1].set_ylabel("Mixture ( = VAF x 2) of Sample 2")

    # ax[0].scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, alpha=0.6,   color=[colorlist[k] for k in membership_left])
    # ax[1].scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, alpha=0.6,   color=[colorlist[k] for k in membership_right])


    # Scatter plot을 각각의 cluster에 대해서 그리기
    for k in range(len(membership_left)):
        try:
            i =  np.where(score_df["answer"] == membership_left[k])[0][0]    # 가장 먼저 V2가 나오는 index
            ax[0].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[ i ]])
        except:
            print ("무슨 문제가 있는거지?\n", membership_left[k], score_df)
    for k in range(len(membership_right)):
        i = np.where(score_df["predicted"] == membership_right[k])[0][0]    # 가장 먼저 3가 나오는 index
        ax[1].scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[ i ]])



    if (dimensionreduction != "SVD") & ( dimensionreduction != "PCA" ):
        #  왼쪽 그림 :  정답을 그리기 (square )
        for samplename_left_character in list(np.unique(membership_left)):
            # membership & np_vaf 정보를 바탕으로
            x_mean = round(np.mean(np_vaf[[x for x in range( len(membership_left)) if membership_left[x] == samplename_left_character]][:, 0] * 2), 3)
            y_mean = round(np.mean(np_vaf[[x for x in range( len(membership_left)) if membership_left[x] == samplename_left_character]][:, 1] * 2), 3)

            ax[0].text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', fontdict = {"fontsize": 16, "fontweight" : "bold"})
            i = np.where(score_df["answer"] == samplename_left_character)[0][0]    # 가장 먼저 V2가 나오는 index
            ax[0].scatter(x_mean, y_mean, marker='s', color=colorlist[ i ], edgecolor='black', s=100,
                        label=str(samplename_left_character) + " : " + str(list(membership_left).count(samplename_left_character)))
            ax[0].legend()


        #  오른쪽 그림 : 내가 구한 EM 을 그리기 (square)
        ##xx = mixture_right.shape[1] - 1 if includeoutlier == "Yes" else mixture_right.shape[1] 
        xx = mixture_right.shape[1] 

        for sample_index in range(xx):   # square 그려주지 말자
            # mixture 정보를 바탕으로
            x_mean = mixture_right[0][sample_index]
            y_mean = mixture_right[1][sample_index]
            ax[1].text(x_mean, y_mean, "{0}".format([x_mean, y_mean]), verticalalignment='top', fontdict = {"fontsize": 16, "fontweight" : "bold"})
            try:
                i = np.where(score_df["predicted"] == sample_index)[0][0]    # 가장 먼저 1이 나오는 index
            except:
                continue
            # ax[1].scatter(x_mean, y_mean, marker='s', color=colorlist[ i ], edgecolor='black', s=100,
            #               label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)))

            if (makeone_index != []) & (makeone_index != None):
                if sample_index in makeone_index:
                    matplotlib.pyplot.scatter(x_mean, y_mean, marker='*', color=colorlist[i], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)) )
                else:
                    matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color=colorlist[i], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)) )
            else:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color=colorlist[i], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(membership_right).count(sample_index)) )



        ax[1].legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)
