import numpy as np
import pandas as pd 
import random

def makedf ( **kwargs ):
    global input_containpos, mutation_id, membership, membership_answer, inputdf, df,  np_vaf, np_BQ, samplename_dict_CharacterToNum, samplename_dict_NumToCharacter

    input_containpos = pd.read_csv( kwargs["INPUT_TSV"],  header = None, sep = "\t") 
    
    if input_containpos.shape[1] == 3: #  If 4th column (BQ) is absent
        input_containpos.columns = ["pos", "sample", "info"]
        input_containpos.astype ({"pos":"str", "sample":"str", "info":"str"})
    elif input_containpos.shape[1] == 4: #  If 4th column (BQ) is present
        input_containpos.columns = ["pos", "sample", "info", "BQ"]
        input_containpos.astype ({"pos":"str", "sample":"str", "info":"str", "BQ":"str"})
    
    
    input_containpos ["cha1"] = "child"         # child, parent
    input_containpos ["cha2"] = "space"       # axis (FN), space 
    samplename_dict_CharacterToNum = {}
    samplename_dict_NumToCharacter = {}

    if kwargs["NUM_MUTATION"] == -1:
        kwargs["NUM_MUTATION"] = input_containpos.shape[0]
        kwargs["RANDOM_PICK"] = input_containpos.shape[0]

    np_vaf = np.zeros(( kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]), dtype = 'float')
    np_BQ = np.zeros(( kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]), dtype = 'float')
    np_BQ.fill(20)  # default BQ is 20
    inputdf = pd.DataFrame (np.zeros(( kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]), dtype = 'object'), columns = ['block' + str(i + 1) for i in range(kwargs["NUM_BLOCK"])])
    mutation_id = []
    membership = []
    depth_list = []
    


    # input_containpos (n * 3)  :  ID (chr_pos), membmership_answer, genomic_profile (Depth1,Alt1,Depth2,Alt2,...,Depth_n,Alt_n)

    depth_col = [[]] * int(len(input_containpos.iloc[0][2].split(","))/2)
    depth_row = []
    for row in range( kwargs["NUM_MUTATION"] ):
        depth_row_mini = []
        mutation_id.append( str(input_containpos.iloc[row][0]) )            # "pos"
        membership.append( str(input_containpos.iloc[row][1]) )           # "sample"
        if "," in str(input_containpos.iloc[row][1]) :
            input_containpos.loc[row,"cha1"] = "parent"

        if str(input_containpos.iloc[row][1]) not in samplename_dict_CharacterToNum.keys():
            samplename_dict_CharacterToNum[str(input_containpos.iloc[row][1])] = int (len(samplename_dict_CharacterToNum))      # {'other': 0, 'V5': 1, 'V3': 2, 'V1': 3}   

        rmv_bracket=input_containpos.iloc[row][2].split(",")     # 3rd column:  112,27,104,9  
        for i in range(0, len(rmv_bracket), 2 ):
            depth = int(rmv_bracket[i])
            alt = int(rmv_bracket[i+1])
            ref = depth - alt

            col = int(i / 2)

            if depth == 0:
                np_vaf[row][col] = 0
                inputdf.iloc[row][col] = "0:0:0"
            else:   
                np_vaf[row][col] = round (alt / depth , 2)
                inputdf.iloc[row][col] = str(depth) + ":" + str(ref) + ":" + str(alt)
                depth_row_mini.append(depth)
                depth_col[col].append(depth)

            if "BQ" in input_containpos.columns: #  If 4th column (BQ) is present
                if kwargs["NUM_BLOCK"] == 1:
                    BQ_input = [input_containpos.iloc[row][3]]     
                else:
                    BQ_input = input_containpos.iloc[row][3].split(",")   # 4th column : 20,20
                for i in range (0, len(BQ_input) ):
                    np_BQ [row][i] = int ( BQ_input[i] )
                    if int ( BQ_input[i] ) == 0:     # Set BQ=20 in axis (FN) variant
                        np_BQ [row][i] = 20 

        depth_row.append (depth_row_mini)


    # Substitute in case of depth = 0  ("0:0:0") to average depth of the other samples.
    for row in range( kwargs["NUM_MUTATION"] ):
        for  i in range(0, len(rmv_bracket), 2 ):
            col = int(i / 2)
            if inputdf.iloc[row][col] == "0:0:0":
                inputdf.iloc[row][col] = str(round(np.mean(depth_col[col]))) + ":" + str(round(np.mean(depth_col[col]))) + ":0"
                input_containpos.loc[row,"cha2"] = "axis"
        depth_list.append(np.mean(depth_row[row]))

    df = [[None] * kwargs["NUM_BLOCK"] for i in range(inputdf.shape[0])]
    for row in range (inputdf.shape[0]):
        for col in range ( kwargs["NUM_BLOCK"] ):
            df[row][col] = {"depth":int(inputdf.iloc[row][col].split(":")[0]), "ref":int(inputdf.iloc[row][col].split(":")[1]), "alt":int(inputdf.iloc[row][col].split(":")[2])}
            if df[row][col]["depth"] == 0:
                print (df[row][col], row, col)

    return kwargs

   


def random_pick_fun(**kwargs):
    global input_containpos, mutation_id, membership, membership_answer, inputdf, df,  np_vaf, np_BQ, samplename_dict_CharacterToNum, samplename_dict_NumToCharacter

    random.seed(kwargs["RANDOM_SEED"])
    random_index = sorted( list ( range (0, kwargs["NUM_MUTATION"])  ))  
    
    input_containpos =  input_containpos.iloc[random_index]
    inputdf  = inputdf.iloc[random_index]
    df = [df[i] for i in random_index]
    np_vaf = np_vaf [random_index]
    membership_answer = [membership[i] for i in random_index]
    mutation_id = [mutation_id[i] for i in random_index]


    t = pd.DataFrame(np_vaf, columns = ["block{0}".format(i) for i in range( kwargs["NUM_BLOCK"])], index = mutation_id)
    t["membership_answer"] = pd.Series(membership_answer, index = mutation_id)
    t.to_csv ("{0}/npvaf.txt".format( kwargs["CLEMENT_DIR"] ), index = True, header=True, sep = "\t")


    samplename_dict_CharacterToNum, cnt = {}, 0
    for k in membership_answer:
        if k not in samplename_dict_CharacterToNum.keys():
            samplename_dict_CharacterToNum[k] = cnt
            cnt = cnt + 1
    
    kwargs["samplename_dict_CharacterToNum"] = samplename_dict_CharacterToNum

    return kwargs




def main (**kwargs):
    global input_containpos, mutation_id, membership, membership_answer, inputdf, df,  np_vaf, np_BQ, samplename_dict_CharacterToNum, samplename_dict_NumToCharacter
    
    kwargs = makedf ( **kwargs )
    
    kwargs = random_pick_fun(**kwargs)

    return (inputdf, df, np_vaf, np_BQ, membership_answer, mutation_id, kwargs )
