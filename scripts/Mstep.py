# Block 1, 2:  data point를 2차원 평면상에 그려보기
def drawfigure_1d_hard(bunch, np_vaf, output_filename, **kwargs):
    import palettable
    import extract
    from scipy.stats import kde
    import matplotlib
    import numpy as np
    import seaborn as sns

    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-whitegrid")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-darkgrid")
        matplotlib.pyplot.style.use("seaborn-white")

    if (bunch.makeone_index == []) | (bunch.makeone_index == None):       # 합쳐서 1 만드는 게 없을 경우 우울하게 만들어주자
        matplotlib.pyplot.style.use("Solarize_Light2")

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includeoutlier == True) & (bunch.fp_index != -1):
        colorlist[bunch.fp_index] = Gr_10[8]

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(
        fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'

    matplotlib.pyplot.figure(figsize=(6, 6))

    max_y = 0

    x = np.linspace(0, 2, 200)

    for k in sorted(list(set(bunch.membership))):        # 각 clone number를 따로 돈다
        np_vaf_new_index, np_vaf_new = extract.npvaf(
            np_vaf, bunch.membership, k)           # membership1  == clone_num 인 것만 추려옴

        try:             # 0만 가득 차 있는 경우 kde를 못 그린다
            kde_np_vaf_new = kde.gaussian_kde(np_vaf_new[:, 0] * 2)
            y = kde_np_vaf_new(x)

            if max_y < np.max(y):
                max_y = np.max(y)

            if k in bunch.makeone_index:
                matplotlib.pyplot.plot(x, y, color=colorlist[k], linewidth=5, label="clone {}".format(k))
            else:
                matplotlib.pyplot.plot( x, y, color=colorlist[k], linewidth=2, label="clone {}".format(k), linestyle="-.")
            matplotlib.pyplot.text(bunch.mixture[0, k], np.max(y) * 1.1, "{} (n={})".format(
                bunch.mixture[0,  k],  np.bincount(bunch.membership)[k]), ha="center", fontsize=15)

        except:
            continue

    if kwargs["OPTION"] in ["Hard", "hard", "Outlier", "outlier"]:
        matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(
            kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(
            kwargs["STEP"], round(bunch.likelihood)), fontsize=20)

    matplotlib.pyplot.title("Sum_child = {}".format(list(np.round(
        np.sum(bunch.mixture[:, bunch.makeone_index], axis=1), 3))), fontsize=12)

    matplotlib.pyplot.axis([0,  np.max(np_vaf[:, :]) * 2.1,  0,  max_y * 1.3])
    matplotlib.pyplot.legend()
    matplotlib.pyplot.xlabel("VAF x 2 of the Sample 1",
                             fontdict={"fontsize": 14})
    matplotlib.pyplot.ylabel("Density", fontdict={"fontsize": 14})

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)


def drawfigure_1d_soft(bunch, np_vaf, output_filename, **kwargs):
    import palettable
    import extract
    from scipy.stats import kde
    import matplotlib
    import numpy as np
    import pandas as pd
    import seaborn as sns

    membership = bunch.membership
    mixture = bunch.mixture
    membership_p_normalize = bunch.membership_p_normalize

    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-whitegrid")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-darkgrid")
        matplotlib.pyplot.style.use("seaborn-white")

    if (bunch.makeone_index == []) | (bunch.makeone_index == None):       # 합쳐서 1 만드는 게 없을 경우 우울하게 만들어주자
        matplotlib.pyplot.style.use("Solarize_Light2")

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]
    if (bunch.includeoutlier == True) & (bunch.fp_index != -1):
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

    if kwargs["OPTION"] in ["Hard", "hard", "Outlier", "outlier"]:
        matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(
            kwargs["STEP"], round(bunch.likelihood)), fontsize=20)
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(
            kwargs["STEP"], round(bunch.likelihood)), fontsize=20)

    matplotlib.pyplot.title("Sum_child = {}".format(list(np.round(
        np.sum(bunch.mixture[:, bunch.makeone_index], axis=1), 3))), fontsize=12)

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


# Block 1, 2:  data point를 2차원 평면상에 그려보기
def drawfigure_2d(bunch, np_vaf, output_filename, **kwargs):
    import palettable
    import matplotlib
    import numpy as np
    import seaborn as sns

    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-white")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-darkgrid")

    if (bunch.makeone_index == []) | (bunch.makeone_index == None):       # 합쳐서 1 만드는 게 없을 경우 우울하게 만들어주자
        matplotlib.pyplot.style.use("Solarize_Light2")

    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(
        fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'

    matplotlib.pyplot.figure(figsize=(6, 6))
    matplotlib.pyplot.axis(
        [0,  np.max(np_vaf[:, :]) * 2.1, 0,  np.max(np_vaf[:, :]) * 2.1])
    matplotlib.pyplot.xlabel("VAF x 2 of the Sample 1",
                             fontdict={"fontsize": 14})
    matplotlib.pyplot.ylabel("VAF x 2 of the Sample 2",
                             fontdict={"fontsize": 14})

    # if kwargs["OPTION"] in ["Hard", "hard", "Outlier", "outlier"]:
    #     matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(kwargs["STEP"], round(bunch.likelihood) ), fontsize = 20)
    # elif kwargs["OPTION"] in ["Soft", "soft"]:
    #     matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(kwargs["STEP"], round(bunch.likelihood) ), fontsize = 20)

    # matplotlib.pyplot.title("Sum_child = {}".format( list (np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) )), fontsize = 12)

    if (bunch.includeoutlier == True) & (bunch.fp_index == -1):
        print("무슨 일이길래 outlier 있다고 하면서 -1일지? {}\t{}\t".format(
            bunch.fp_index_record, bunch.makeone_index_record))

    if (bunch.includeoutlier == True) & (bunch.fp_index != -1):
        # 맨 마지막 번호의 색깔번호 (Outlier 번호)
        outlier_color_num = samplename_dict[bunch.fp_index]
        colorlist[outlier_color_num] = Gr_10[10]

    matplotlib.pyplot.scatter(np_vaf[:, 0] * 2, np_vaf[:, 1] * 2, color=[
                              colorlist[samplename_dict[k]] for k in bunch.membership])

    for sample_index, sample in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue
        # mixture 정보를 바탕으로
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        matplotlib.pyplot.text(x_mean, y_mean, "{0}".format(
            [x_mean, y_mean]), verticalalignment='top', fontsize=15)

        if (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(
                    sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100, label="cluster" + str(
                    sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            matplotlib.pyplot.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(
                sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        matplotlib.pyplot.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)


# Block 1, 2:  data point를 2차원 평면상에 그려보기
def drawfigure_2d_soft(bunch, np_vaf, output_filename, **kwargs):
    import palettable
    import matplotlib
    import numpy as np
    import seaborn as sns
    from sklearn.decomposition import TruncatedSVD, PCA

    if kwargs["OPTION"] in ["Hard", "hard"]:
        matplotlib.pyplot.style.use("seaborn-white")
    elif kwargs["OPTION"] in ["Soft", "soft"]:
        matplotlib.pyplot.style.use("seaborn-white")
    elif kwargs["OPTION"] in ["Outlier", "outlier"]:
        matplotlib.pyplot.style.use("Solarize_Light2")

    NUM_MUTATION = len(bunch.membership)
    NUM_CLONE = kwargs["NUM_CLONE"]
    NUM_BLOCK = kwargs["NUM_BLOCK"]

    samplename_dict = {k: k for k in range(0, bunch.mixture.shape[1])}

    tabl = palettable.tableau.Tableau_20.mpl_colors
    Gr_10 = palettable.scientific.sequential.GrayC_20.mpl_colors
    colorlist = [i for i in tabl]

    font_dir = "/home/goldpm1/miniconda3/envs/cnvpytor/lib/python3.7/site-packages/matplotlib/mpl-data/fonts/ttf/"
    font_dirs = matplotlib.font_manager.findSystemFonts(
        fontpaths=font_dir, fontext='ttf')
    for font in font_dirs:
        matplotlib.font_manager.fontManager.addfont(font)
    matplotlib.rcParams["font.family"] = 'arial'

    if bunch.includeoutlier == True:
        # 맨 마지막 번호의 색깔번호 (Outlier 번호)
        outlier_color_num = samplename_dict[bunch.fp_index]
        colorlist[outlier_color_num] = Gr_10[10]

    matplotlib.pyplot.figure(figsize=(6, 6))

    dimensionreduction = ""
    if bunch.mixture.shape[0] > 2:  # 3차원 이상이면
        dimensionreduction = "SVD"

    if dimensionreduction == "SVD":
        print("SVD → 2D")
        tsvd = TruncatedSVD(n_components=2)
        # mixture 도 square로 찍어주려면 차원축소 변환을 잘 해줘야 한다
        tsvd.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = tsvd.transform(np.concatenate(
            [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = tsvd.transform(np.concatenate(
            [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0])
                               * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        matplotlib.pyplot.xlabel("SVD1")
        matplotlib.pyplot.ylabel("SVD2")
    elif dimensionreduction == "PCA":
        print("PCA → 2D")
        pca = PCA(n_components=2)
        # mixture 도 square로 찍어주려면 차원축소 변환을 잘 해줘야 한다
        pca.fit(np.concatenate([np_vaf, bunch.mixture]))
        np_vaf = pca.transform(np.concatenate(
            [np_vaf, bunch.mixture]))[:-NUM_BLOCK]
        mixture = pca.transform(np.concatenate(
            [np_vaf, bunch.mixture]))[-NUM_BLOCK:]
        matplotlib.pyplot.axis([np.min(np_vaf[:, 0]) * 2.1,  np.max(np_vaf[:, 0])
                               * 2.1,  np.min(np_vaf[:, 1]) * 2.1,  np.max(np_vaf[:, 1]) * 2.1])
        matplotlib.pyplot.xlabel("PC1")
        matplotlib.pyplot.ylabel("PC2")
    else:
        matplotlib.pyplot.axis(
            [0,  np.max(np_vaf[:, :]) * 2.1,  0,  np.max(np_vaf[:, :]) * 2.1])
        matplotlib.pyplot.xlabel("Feature 1 : VAF x 2 of the Sample 1")
        matplotlib.pyplot.ylabel("Feature 2 : VAF x 2 of the Sample 2")

    # if kwargs["OPTION"] in ["Soft", "soft"]:
    #     matplotlib.pyplot.suptitle("Step {},  Total prob = {}".format(kwargs["STEP"], round(bunch.likelihood) ), fontsize = 20)

    # matplotlib.pyplot.title("Sum_child = {}".format( list ( np.round( np.sum (bunch.mixture[:, bunch.makeone_index], axis = 1), 3) ) ), fontsize = 12)

    if bunch.includeoutlier == "Yes":
        # 맨 마지막 번호의 색깔번호 (Outlier 번호)
        outlier_color_num = samplename_dict[bunch.fp_index]
        colorlist[outlier_color_num] = Gr_10[10]

    # 각 점을 일일이 색칠하기
    #print (NUM_CLONE, bunch.membership_p_normalize [ np.where(bunch.membership == (NUM_CLONE - 1))[0]])
    for j in range(0, NUM_CLONE):
        if j == bunch.fp_index:
            for k in bunch.outlier_index:
                matplotlib.pyplot.scatter(
                    np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=Gr_10[10], s = 150)
        else:
            for k in range(NUM_MUTATION):
                if bunch.membership_p_normalize[k, j] > 0.8:
                    matplotlib.pyplot.scatter(
                        np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=1, color=[colorlist[samplename_dict[j]]], s = 150)
                elif bunch.membership_p_normalize[k, j] > 0.1:
                    matplotlib.pyplot.scatter(
                        np_vaf[k, 0] * 2, np_vaf[k, 1] * 2, alpha=bunch.membership_p_normalize[k, j], color=[colorlist[samplename_dict[j]]], s = 150)
                else:        # 10%도 기여하지 못한다면 칠하지 말고 넘어가자
                    continue

    for sample_index, sample_key in enumerate(samplename_dict):
        if (sample_index >= bunch.mixture.shape[1]):
            continue

        # mixture 정보를 바탕으로 (얘는 앞부터 순차적으로 해야 하니까 sample_index가 맞다)
        # 얘도 2차원으로 차원축소했으니 걱정없다
        x_mean = round(bunch.mixture[0][sample_index], 2)
        y_mean = round(bunch.mixture[1][sample_index], 2)
        matplotlib.pyplot.text(x_mean, y_mean, "{0}".format(
            [x_mean, y_mean]), verticalalignment='top', fontsize=15)

        if (bunch.makeone_index != []) & (bunch.makeone_index != None):
            if sample_index in bunch.makeone_index:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='*', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(
                    sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
            else:
                matplotlib.pyplot.scatter(x_mean, y_mean, marker='s', color=colorlist[sample_index], edgecolor='black', s=100, label="cluster" + str(
                    sample_index) + " : " + str(list(bunch.membership).count(sample_index)))
        else:
            matplotlib.pyplot.scatter(x_mean, y_mean, marker='+', color=colorlist[sample_index], edgecolor='black', s=200, label="cluster" + str(
                sample_index) + " : " + str(list(bunch.membership).count(sample_index)))

        # membership & np_vaf 정보를 바탕으로
        # x_mean = round(np.mean(np_vaf[[x for x in range( len(membership)) if membership[x] == sample_key]][:, 0] * 2), 2)
        # y_mean = round(np.mean(np_vaf[[x for x in range( len(membership)) if membership[x] == sample_key]][:, 1] * 2), 2)
        # matplotlib.pyplot.text(x_mean, y_mean, "{0}".format([x_mean, y_mean]), verticalalignment='top')
        # matplotlib.pyplot.scatter(x_mean, y_mean, marker='^', color=colorlist[sample_value],  s=150,  label="hard : cluster" + str(sample_key) + " : " + str(list(membership).count(sample_key)))
        matplotlib.pyplot.legend()

    if output_filename != "NotSave":
        matplotlib.pyplot.savefig(output_filename)


def main(df, np_vaf, step, option, **kwargs):  # 새로운 MIXTURE 정하는 과정
    import math
    import isparent
    import visualizationsinglesoft
    import numpy as np

    NUM_BLOCK, NUM_CLONE, NUM_MUTATIN = kwargs["NUM_BLOCK"], kwargs["NUM_CLONE"], kwargs["NUM_MUTATION"]
    NUM_MUTATION = kwargs["RANDOM_PICK"]

    kwargs["OPTION"] = option

    if option in ["Hard", "hard"]:
        ############################### HARD CLUSTERING ##############################
        for j in range(NUM_CLONE):
            # membership == j 인 index를 구하기
            ind_list = np.where(step.membership == j)[0]
            for i in range(NUM_BLOCK):
                sum_depth, sum_alt = 0, 0
                for ind in ind_list:       # depth, alt를 다 더하기
                    sum_depth = sum_depth + df[ind][i]["depth"]
                    sum_alt = sum_alt + df[ind][i]["alt"]

                #print (i, j, sum_depth, sum_alt)
                # j번째 clone만을 생각한 이상적인 분율을 일단 assign
                step.mixture[i][j] = round(
                    (sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0

        # Block당 Mixture의 합이 1이 되도록 재조정
        if kwargs["adjustment"] in ["True", "true", True]:
            for i in range(NUM_BLOCK):
                sum = np.sum(step.mixture[i])
                # 만약 sum이 0이면 분모에 0이 들어갈 수 있으니까...
                step.mixture[i] = np.round(
                    step.mixture[i] / sum, 4) if sum != 0 else 0

        elif kwargs["adjustment"] in ["Half", "half"]:
            previous_fp_index = step.fp_index
            step.makeone_index, p_list, step.fp_index = isparent.makeone(step, **kwargs)

            if step.fp_index != -1:  # Outlier clone이 있다면
                step.includeoutlier = True
                step.outlier_index = list(np.where(step.membership == step.fp_index)[0])
                if previous_fp_index != step.fp_index:
                    step.likelihood = -9999999
                    #print ("FP clone을 갑자기 바꿨기 때문에 likelihood 계산이 불가하다")

            else:   # Outlier clone이 없다면
                step.includeoutlier = False
                step.outlier_index = []

            #print ("step_index = {}\tstep.fp_index={}\tstep.includeoutlier={}".format(kwargs["STEP_TOTAL"], step.fp_index , step.includeoutlier))

            if NUM_CLONE == 1:
                step.mixture = np.array([[1.0]] * kwargs["NUM_BLOCK"])

            if step.makeone_index == []:
                step.likelihood = -9999999
                #print ("여러 조건 (parent, child) 때문에 합쳐서 1을 못 만든다\t {}".format( list(step.mixture)) )

            if kwargs["STEP"] <= 4:
                if step.makeone_index != []:
                    for i in range(NUM_BLOCK):
                        sum = 0
                        for j in range(NUM_CLONE):
                            if j in step.makeone_index:                # 1과 가장 가까운 clone들만 붙잡아줌
                                sum = sum + step.mixture[i][j]
                        # 만약 sum이 0이면 분모에 0이 들어갈 수 있으니까...
                        step.mixture[i] = np.round(step.mixture[i] / sum, 4) if sum != 0 else 0

            if (kwargs["VERBOSE"] in ["True", "T"]) | (int(str(kwargs["VERBOSE"])) >= 1):
                if (kwargs["NUM_BLOCK"] == 1):
                    drawfigure_1d_hard(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], **kwargs)
                if (kwargs["NUM_BLOCK"] >= 2):
                    drawfigure_2d(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(hard)." + kwargs["IMAGE_FORMAT"], **kwargs)
        #############################################################################

    ################################ SOFT CLUSTERING ##############################

    if option in ["Soft", "soft"]:
        #print ("\t\ta. Mixture (before soft clustering) : {}". format(list(step.mixture)))

        makeone_index_i = []
        for k in range(NUM_MUTATION):
            if step.membership[k] in step.makeone_index:
                makeone_index_i.append(k)

        for j in range(NUM_CLONE):
            # child가 아닌 애들 (parent, fp)은 mixutre도 그냥 hard처럼 mixture를 계산해준다
            if j not in step.makeone_index:
                for i in range(NUM_BLOCK):
                    sum_depth, sum_alt = 0, 0
                    # depth, alt를 다 더하기
                    for ind in np.where(step.membership == j)[0]:
                        sum_depth = sum_depth + df[ind][i]["depth"]
                        sum_alt = sum_alt + df[ind][i]["alt"]
                    # j번째 clone만을 생각한 이상적인 분율을 일단 assign
                    step.mixture[i][j] = round(
                        (sum_alt * 2) / sum_depth, 2) if sum_depth != 0 else 0

            elif j in step.makeone_index:  # child clone들만 soft clustering 한다
                for i in range(NUM_BLOCK):   # 각 block에 대해 weighted mean을 계산
                    vaf, weight = np.zeros(NUM_MUTATION, dtype="float"), np.zeros(
                        NUM_MUTATION, dtype="float")
                    for k in range(NUM_MUTATION):
                        vaf[k] = int(df[k][i]["alt"]) / int(df[k][i]["depth"])

                        if step.membership[k] in step.makeone_index:
                            weight[k] = math.pow(10, step.membership_p[k][j])

                    step.mixture[i][j] = round(np.average(
                        vaf[makeone_index_i], weights=weight[makeone_index_i]), 4) * 2

        #print ("\t\tb. Mixture (after soft clustering) : {}". format(list(step.mixture)))

        # Block당 Mixture의 합이 1이 되도록 재조정
        if kwargs["adjustment"] in ["True", "true", True]:
            for i in range(NUM_BLOCK):
                sum = np.sum(step.mixture[i])
                # 만약 sum이 0이면 분모에 0이 들어갈 수 있으니까...
                step.mixture[i] = np.round(
                    step.mixture[i] / sum, 4) if sum != 0 else 0

        elif kwargs["adjustment"] in ["Half", "half"]:
            if step.fp_index != -1:  # Outlier clone이 있다면
                step.includeoutlier = True
                step.outlier_index = list(np.where(step.membership == step.fp_index)[0])
            else:   # Outlier clone이 없다면
                step.includeoutlier = False
                step.outlier_index = []

            if NUM_CLONE == 1:
                step.mixture = np.array([[1.0]] * kwargs["NUM_BLOCK"])

            if step.makeone_index == []:
                step.likelihood = -9999999
                # print ("여러 조건 (parent, child) 때문에 합쳐서 1을 못 만든다\t", list(mixture))

        if kwargs["STEP"] <= 5:
            if step.makeone_index != []:
                for i in range(NUM_BLOCK):
                    sum = 0
                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:                # 1과 가장 가까운 clone들만 붙잡아줌
                            sum = sum + step.mixture[i][j]
                    # 만약 sum이 0이면 분모에 0이 들어갈 수 있으니까...

                    for j in range(NUM_CLONE):
                        if j in step.makeone_index:                # 1과 가장 가까운 clone들만 붙잡아줌
                            step.mixture[i][j] = np.round(
                                step.mixture[i][j] / sum, 4) if sum != 0 else 0

        #print ("\t\tc. Mixture (after normalization) : {}". format(list(step.mixture)))

        if (kwargs["VERBOSE"] in ["True", "T"]) | (int(str(kwargs["VERBOSE"])) >= 1):
            step.membership_p_normalize = np.zeros(
                (NUM_MUTATION, step.membership_p.shape[1]), dtype="float64")
            for k in range(NUM_MUTATION):
                if k in step.outlier_index:
                    step.membership_p_normalize[k] = np.zeros(
                        step.membership_p_normalize.shape[1], dtype="float64")  # 1 (FP_index) 0 0 0 0 0 으로 만들어준다
                    # 1 (FP_index) 0 0 0 0 0 으로 만들어준다
                    step.membership_p_normalize[k][step.fp_index] = 1
                else:
                    step.membership_p_normalize[k] = np.round(np.power(10, step.membership_p[k])/np.power(10, step.membership_p[k]).sum(axis=0, keepdims=1), 2)   # 로그를 취해으나 다시 지수를 취해준다
                    if step.fp_index != -1:     # fp가 있을 때에만...  fp가 없는데도 -1 column이 0으로 되버리는 참사
                        # 0  (FP_index) 0.2 0.7 0.1 0 0 으로 만들어준다
                        step.membership_p_normalize[k][step.fp_index] = 0

            if (kwargs["NUM_BLOCK"] == 1):
                drawfigure_1d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], **kwargs)
            if (kwargs["NUM_BLOCK"] == 2):
                drawfigure_2d_soft(step, np_vaf, kwargs["CLEMENT_DIR"] + "/trial/clone" + str(kwargs["NUM_CLONE_NOMINAL"]) + "." + str(kwargs["TRIAL"]) + "-" + str(kwargs["STEP_TOTAL"]) + "(soft)." + kwargs["IMAGE_FORMAT"], **kwargs)

    #############################################################################

    return step
