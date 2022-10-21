import numpy as np
import pandas as pd
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
import scipy.stats as stats



#  Extract only outlier  (df, np_vaf)
def main (df, np_vaf, membership, clone_no):
    df_new = []
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :  
        df_new.append(df[i])
        df_new_index.append(i)
    return df_new, df_new_index, np_vaf[df_new_index]


# Extract only outlier  (np_vaf) :  needed in visualization
def npvaf (np_vaf, membership, clone_no):
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :  
        df_new_index.append(i)
    return df_new_index, np_vaf[df_new_index]
