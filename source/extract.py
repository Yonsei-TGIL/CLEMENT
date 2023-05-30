import numpy as np
import pandas as pd
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
import scipy.stats as stats



#  Outlier 들만 추출해주기 (df, np_vaf)
def main (df, np_vaf, membership, clone_no):
    df_new = []
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :       # membership 번호가 맨 끝번이면 outlier로 처리
        df_new.append(df[i])
        df_new_index.append(i)
    return df_new, df_new_index, np_vaf[df_new_index]


# Outlier 들만 추출해주기 (np_vaf)  visualization에서 필요함
def npvaf (np_vaf, membership, clone_no):
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :       # membership 번호가 맨 끝번이면 outlier로 처리
        df_new_index.append(i)
    return df_new_index, np_vaf[df_new_index]
