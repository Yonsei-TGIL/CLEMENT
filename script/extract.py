def main (df, np_vaf, membership, clone_no):
    df_new = []
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :   
        df_new.append(df[i])
        df_new_index.append(i)
    return df_new, df_new_index, np_vaf[df_new_index]


# visualization
def npvaf (np_vaf, membership, clone_no):
    df_new_index= []
    for i in [i for i in range(len(membership)) if membership[i] == clone_no] :       
        df_new_index.append(i)
    return df_new_index, np_vaf[df_new_index]
