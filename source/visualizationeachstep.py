import palettable
import matplotlib
import extract
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import kde
from sklearn.decomposition import TruncatedSVD, PCA

# Block 1, 2:  data point를 2차원 평면상에 그려보기
def drawfigure_1d_hard(bunch, np_vaf, output_filename, **kwargs):
    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-whitegrid")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-darkgrid")

    if (bunch.makeone_index == []) | (bunch.makeone_index == None):       # 합쳐서 1 만드는 게 없을 경우 우울하게 만들어주자
        matplotlib.pyplot.style.use("Solarize_Light2")
    if (bunch.makeone_prenormalization == False):
        matplotlib.pyplot.style.use("Solarize_Light2")

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includefp == True) & (bunch.fp_index != -1):
        colorlist[bunch.fp_index] = Gr_10[8]

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'

    matplotlib.pyplot.figure(figsize=(6, 6))

    max_y = 0

    x = np.linspace(0, 2, 200)

    for k in sorted(list(set(bunch.membership))):        # 각 clone number를 따로 돈다
        np_vaf_new_index, np_vaf_new = extract.npvaf(np_vaf, bunch.membership, k)           # membership1  == clone_num 인 것만 추려옴

        try:             # 0만 가득 차 있는 경우 kde를 못 그린다
            kde_np_vaf_new = kde.gaussian_kde(np_vaf_new[:, 0] * 2)
            weight = len(np_vaf_new) / kwargs["NUM_MUTATION"]  # weight를 줘서 total density 적분값이 1이 되도록 하자
            y = kde_np_vaf_new(x) * weight

            if max_y < np.max(y):
                max_y = np.max(y)

            if k in bunch.makeone_index:
                matplotlib.pyplot.plot(x, y, color=colorlist[k], linewidth=5, label="clone {}".format(k))
            else:
                matplotlib.pyplot.plot(x, y, color=colorlist[k], linewidth=2, label="clone {}".format(k), linestyle="-.")
            matplotlib.pyplot.text(bunch.mixture[0, k], np.max(y) * 1.1, "{} (n={})".format(bunch.mixture[0,  k],  np.bincount(bunch.membership)[k]), ha="center", fontdict = {"fontsize": 16, "fontweight" : "bold"})

        except:
            continue

    matplotlib.pyplot.suptitle("clone{}.{}-{},  Total prob = {}".format(kwargs["NUM_CLONE"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    # if kwargs["OPTION"] in ["Hard", "hard", "fp", "fp"]:
    #     matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    # elif kwargs["OPTION"] in ["Soft", "soft"]:
    #     matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(kwargs["STEP"], round(bunch.likelihood)), fontsize=20)

    matplotlib.pyplot.title("Sum_child = {}".format(list(np.round(np.sum(bunch.mixture[:, bunch.makeone_index], axis=1), 3))), fontsize=12)

    matplotlib.pyplot.axis([0,  np.max(np_vaf[:, :]) * 2.1,  0,  max_y * 1.3])
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xlabel("VAF x 2 of the Sample 1", fontdict={"fontsize": 14})
    matplotlib.pyplot.ylabel("Density", fontdict={"fontsize": 14})

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)





def drawfigure_1d_soft(bunch, np_vaf, output_filename, **kwargs):
    membership = bunch.membership
    mixture = bunch.mixture
    membership_p_normalize = bunch.membership_p_normalize

    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-whitegrid")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-darkgrid")

    if (bunch.makeone_index == []) | (bunch.makeone_index == None):       # 합쳐서 1 만드는 게 없을 경우 우울하게 만들어주자
        matplotlib.pyplot.style.use("Solarize_Light2")
    if (bunch.makeone_prenormalization == False):
        matplotlib.pyplot.style.use("Solarize_Light2")

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includefp == True) & (bunch.fp_index != -1):
        colorlist[bunch.fp_index] = Gr_10[8]

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(
        fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'

    matplotlib.pyplot.figure(figsize=(6, 6))
    matplotlib.pyplot.xlim(0, 1)
    matplotlib.pyplot.xlabel("VAF x 2 (= Mixture)")
    matplotlib.pyplot.ylabel("count (weighted)")

    matplotlib.pyplot.suptitle("clone{}.{}-{},  Total prob = {}".format(kwargs["NUM_CLONE"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    # if kwargs["OPTION"] in ["Hard", "hard", "fp", "fp"]:
    #     matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    # elif kwargs["OPTION"] in ["Soft", "soft"]:
    #     matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(kwargs["STEP"], round(bunch.likelihood)), fontsize=20)

    matplotlib.pyplot.title("Sum_child = {}".format(list(np.round(np.sum(bunch.mixture[:, bunch.makeone_index], axis=1), 3))), fontsize=12)

    for k in sorted(list(set(bunch.membership))):
        if k in bunch.makeone_index:
            sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, k], "linewidth": 1.4, "edgecolor": "black"}, kde_kws={
                         "linewidth": 5, "color": "gray"}, color=colorlist[k], kde=False, bins=50, label="cluster {}  (weighted vaf = {}, mixture = {})".format(k,  str(round((mixture[0][k]) / 2, 2)), str(round((mixture[0][k]), 2))))
        else:
            sns.distplot(pd.DataFrame(np_vaf[:, 0] * 2, columns=["vaf"])["vaf"], hist_kws={"weights": membership_p_normalize[:, k], "rwidth": 0.8}, color=colorlist[k],
                         kde=False, bins=50, label="cluster {}  (weighted vaf = {}, mixture = {})".format(k,  str(round((mixture[0][k]) / 2, 2)), str(round((mixture[0][k]), 2))))
    matplotlib.pyplot. legend()

    matplotlib.pyplot.xlabel("VAF * 2 (= Mixture)")
    matplotlib.pyplot.ylabel("count (weighted)")

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)







def drawfigure_2d(bunch, np_vaf, output_filename, **kwargs):
    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-white")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-darkgrid")

    if (bunch.makeone_index == []) | (bunch.makeone_index == None):       # 합쳐서 1 만드는 게 없을 경우 우울하게 만들어주자
        matplotlib.pyplot.style.use("Solarize_Light2")
    if (bunch.makeone_prenormalization == False):
        matplotlib.pyplot.style.use("Solarize_Light2")

    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'

    matplotlib.pyplot.figure(figsize=(6, 6))
    matplotlib.pyplot.axis( [0,  np.max(np_vaf[:, :]) * 2.1, 0,  np.max(np_vaf[:, :]) * 2.1])
    matplotlib.pyplot.xlabel("VAF x 2 of the Sample 1", fontdict={"fontsize": 14})
    matplotlib.pyplot.ylabel("VAF x 2 of the Sample 2", fontdict={"fontsize": 14})

    matplotlib.pyplot.suptitle("clone{}.{}-{},  Total prob = {}".format(kwargs["NUM_CLONE"], kwargs["TRIAL"], kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) )), fontsize = 12)

    if (bunch.includefp == True) & (bunch.fp_index == -1):
        print("무슨 일이길래 fp 있다고 하면서 -1일지? {}\t{}\t".format(bunch.fp_index_record, bunch.makeone_index_record))

    if (bunch.includefp == True) & (bunch.fp_index != -1):
        fp_color_num = samplename_dict[bunch.fp_index]
        colorlist[fp_color_num] = Gr_10[10]

    matplotlib.pyplot.scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, color=[colorlist[samplename_dict[k]] for k in bunch.membership])

    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue
        # mixture 정보를 바탕으로
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        matplotlib.pyplot.text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', fontdict = {"fontsize": 16, "fontweight" : "bold"} )

        if (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            matplotlib.pyplot.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        matplotlib.pyplot.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)





# Block 1, 2:  data point를 2차원 평면상에 그려보기
def drawfigure_2d_soft(bunch, np_vaf, output_filename, **kwargs):
    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-white")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        #matplotlib.pyplot.style.use("seaborn-darkgrid")   # 디버깅용이면 이게 편한데
        matplotlib.pyplot.style.use("seaborn-white")    # figure 작업하려면 이게 편하다
    elif kwargs["OPTION"] in ["fp", "fp"]:
        matplotlib.pyplot.style.use("Solarize_Light2")

    NUM_MUTATION = len(bunch.membership)
    NUM_CLONE = kwargs["NUM_CLONE"]
    NUM_BLOCK = kwargs["NUM_BLOCK"]

    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'

    if bunch.includefp == True:
        fp_color_num = samplename_dict[bunch.fp_index]
        colorlist[fp_color_num] = Gr_10[10]

    matplotlib.pyplot.figure(figsize=(6, 6))

    dimensionreduction = ""
    if bunch.mixture.shape[0] > 2:  # 3차원 이상이면
        dimensionreduction = "SVD"

    if dimensionreduction == "SVD":
        print("SVD → 2D")
        tsvd = TruncatedSVD(n_components=2)
        # mixture 도 square로 찍어주려면 차원축소 변환을 잘 해줘야 한다
        tsvd.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = tsvd.transform(np.concatenate(  [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = tsvd.transform(np.concatenate( [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0]) * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        matplotlib.pyplot.xlabel("SVD1")
        matplotlib.pyplot.ylabel("SVD2")
    elif dimensionreduction == "PCA":
        print("PCA → 2D")
        pca = PCA(n_components=2)
        # mixture 도 square로 찍어주려면 차원축소 변환을 잘 해줘야 한다
        pca.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = pca.transform(np.concatenate( [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = pca.transform(np.concatenate(  [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0])  * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        matplotlib.pyplot.xlabel("PC1")
        matplotlib.pyplot.ylabel("PC2")
    else:
        matplotlib.pyplot.axis( [0,  np.max(np_vaf[:, :]) * 2.1,  0,  np.max(np_vaf[:, :]) * 2.1] )
        matplotlib.pyplot.xlabel("Feature 1 : VAF x 2 of the Sample 1")
        matplotlib.pyplot.ylabel("Feature 2 : VAF x 2 of the Sample 2")


    matplotlib.pyplot.suptitle("clone{}.{}-{},  Total prob = {}".format(kwargs["NUM_CLONE"], kwargs["TRIAL"], kwargs["STEP_TOTAL"], round(bunch.likelihood)), fontsize=20)
    matplotlib.pyplot.title("sum_child = {}".format( list ( np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ) ), fontsize = 12)

    if bunch.includefp == "Yes":
        fp_color_num = samplename_dict[bunch.fp_index]
        colorlist[fp_color_num] = Gr_10[10]

    # 각 점을 일일이 색칠하기
    
    
    for j in range(0, NUM_CLONE):
        if j == bunch.fp_index:
            for k in bunch.fp_member_index:
                matplotlib.pyplot.scatter(np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=Gr_10[10])
        else:
            for k in range(NUM_MUTATION):
                if bunch.membership_p_normalize[k, j] > 0.8:
                    matplotlib.pyplot.scatter( np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[samplename_dict[j]]])
                elif bunch.membership_p_normalize[k, j] > 0.1:
                    matplotlib.pyplot.scatter( np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=bunch.membership_p_normalize[k, j], color=[colorlist[samplename_dict[j]]])
                else:        # 10%도 기여하지 못한다면 칠하지 말고 넘어가자
                    continue

    for sample_index, sample_key in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        # mixture 정보를 바탕으로 (얘는 앞부터 순차적으로 해야 하니까 sample_index가 맞다)
        # 얘도 2차원으로 차원축소했으니 걱정없다
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        matplotlib.pyplot.text(x_mean, y_mean, "{0}".format( [x_mean, y_mean]), verticalalignment='top', fontdict = {"fontsize": 16, "fontweight" : "bold"})

        if (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            matplotlib.pyplot.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        # membership & np_vaf 정보를 바탕으로
        # x_mean = round(np.mean(np_vaf[[x for x in range( len(membership)) if membership[x] == sample_key]][:, 0] * 2), 2)
        # y_mean = round(np.mean(np_vaf[[x for x in range( len(membership)) if membership[x] == sample_key]][:, 1] * 2), 2)
        # matplotlib.pyplot.text(x_mean, y_mean, "{0}".format([x_mean, y_mean]), verticalalignment='top')
        # matplotlib.pyplot.scatter(x_mean, y_mean, marker='^', color=colorlist[sample_value],  s=150,  label="hard : cluster" + str(sample_key) + " : " + str(list(membership).count(sample_key)))
        matplotlib.pyplot.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)
