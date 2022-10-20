import numpy as np
import pandas as pd 


def makedf (**kwargs):
    global NUM_BLOCK_INPUT, NUM_BLOCK, RANDOM_PICK, AXIS_RATIO, PARENT_RATIO, FP_RATIO, INPUT_TSV, mixture_answer, membership_answer, samplename_dict, DEPTH_CUTOFF
    global  df, inputdf, input_containpos, np_vaf, membership, mutation_id, depth_list, NUM_MUTATION, samplename_dict
    global membership_answer, mutation_id, mixture_answer, parent_type, parent_type_selected

    input_containpos = pd.read_csv(INPUT_TSV,  header = None, names =["pos", "sample", "info"], sep = "\t") 
    input_containpos ["cha1"] = "exclusive"  # Exclusive ( = child = independent clone),  Parent ( = Superclone = mixture of two subclone: if "," is found in inputdata
    input_containpos ["cha2"] = "space"       # Space (all samples have depth >= DEPTH_CUTOFF  or Axis if at least depth of 1 sample is zero)
    samplename_dict = {}
    NUM_MUTATION = input_containpos.shape[0]

    np_vaf = np.zeros((NUM_MUTATION, NUM_BLOCK_INPUT), dtype = 'float')
    inputdf = pd.DataFrame (np.zeros((NUM_MUTATION, NUM_BLOCK_INPUT), dtype = 'object'), columns = ['block' + str(i + 1) for i in range(NUM_BLOCK_INPUT)])
    mutation_id = []
    membership = []
    depth_list = []
    
    #print ( input_containpos )

    # INPUT_FORMAT (n * 3) :   ID (chr_pos), membmership(정답 set 일 경우),  NUM_BLOCK_INPUT(3)만큼의 depth, alt 정보

    depth_col = [[]] * int(len(input_containpos.iloc[0][2].split(","))/2)
    depth_row = []
    for row in range(NUM_MUTATION):
        depth_row_mini = []
        mutation_id.append( str(input_containpos.iloc[row][0]) )            # "pos"
        membership.append( str(input_containpos.iloc[row][1]) )           # "sample"
        if "," in str(input_containpos.iloc[row][1]) :
            input_containpos.loc[row,"cha1"] = "parent"

        if str(input_containpos.iloc[row][1]) not in samplename_dict.keys():
            samplename_dict[str(input_containpos.iloc[row][1])] = int (len(samplename_dict))      # {'other': 0, 'V5': 1, 'V3': 2, 'V1': 3}           # To numericalize each sample name (e.g. V5 → 1)

        rmv_bracket=input_containpos.iloc[row][2].split(",")
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
        depth_row.append (depth_row_mini)

    # "0.0.0"을 그대로 놔둘 수 없다.  평균 depth로 갈음해서 바꿔 넣는다  (alt는 0으로 유지)

    for row in range(NUM_MUTATION):
        for  i in range(0, len(rmv_bracket), 2 ):
            col = int(i / 2)
            if inputdf.iloc[row][col] == "0:0:0":
                inputdf.iloc[row][col] = str(round(np.mean(depth_col[col]))) + ":" + str(round(np.mean(depth_col[col]))) + ":0"
                input_containpos.loc[row,"cha2"] = "axis"
        depth_list.append(np.mean(depth_row[row]))

    df = [[None] * NUM_BLOCK for i in range(inputdf.shape[0])]
    for row in range (inputdf.shape[0]):
        for col in range (NUM_BLOCK):
            df[row][col] = {"depth":int(inputdf.iloc[row][col].split(":")[0]), "ref":int(inputdf.iloc[row][col].split(":")[1]), "alt":int(inputdf.iloc[row][col].split(":")[2])}
            if df[row][col]["depth"] == 0:
                print (df[row][col], row, col)




    # Print count table from the Input data of which depth exceeds DEPTH_CUTOFF
    df_counts = pd.DataFrame (np.unique ( [membership[i] for i in [i for i in range (0, NUM_MUTATION) if (  (depth_list[i] > DEPTH_CUTOFF) ) ] ], return_counts = True ))

    for s_index, s in enumerate(df_counts.iloc[0]):
        samplename_dict[s] = s_index

    for i in range(df_counts.shape[0]):
        for j in range(df_counts.shape[1]) :
            print (df_counts.iloc[i][j], end ="\t")
        print ("")


    
    
    


def random_pick_fun(**kwargs):
    import random
    global NUM_BLOCK_INPUT, NUM_BLOCK, RANDOM_PICK, mixture_answer, membership_answer, DEPTH_CUTOFF
    global  df, inputdf, input_containpos, np_vaf, membership, mutation_id, depth_list, NUM_MUTATION, NUM_CLONE, NUM_BLOCK, samplename_dict
    global membership_answer, mutation_id, mixture_answer, random_index

    random.seed(kwargs["RANDOM_SEED"])
    random_index = sorted( random.sample (range (0, NUM_MUTATION), kwargs["RANDOM_PICK"])  )             # 다 합치면 RADOM_PICK 개수가 되겠지

    # Pick only RANDOM_SEED out of total input
    input_containpos =  input_containpos.iloc[random_index]
    inputdf  = inputdf.iloc[random_index]
    df = [df[i] for i in random_index]
    np_vaf = np_vaf[random_index]
    membership_answer = [membership[i] for i in random_index]
    mutation_id = [mutation_id[i] for i in random_index]


def print_record (**kwargs):
    global NUM_BLOCK_INPUT, NUM_BLOCK, RANDOM_PICK, mixture_answer, membership_answer, DEPTH_CUTOFF
    global  df, inputdf, input_containpos, np_vaf, membership, mutation_id, depth_list, NUM_MUTATION, NUM_CLONE, NUM_BLOCK, samplename_dict
    global membership_answer, mutation_id, mixture_answer

    # print np_vaf + membership to csv  ( PICK only RANDOM_PICK )
    t = pd.DataFrame(np_vaf, columns = ["block{0}".format(i) for i in range(NUM_BLOCK_INPUT)], index = mutation_id)
    t["membership_answer"] = pd.Series(membership_answer, index = mutation_id)
    t.to_csv ("{}.npvaf".format( kwargs["NPVAF_DIR"] ), index = True, header=True, sep = "\t")


    samplename_dict, cnt = {}, 0
    for k in membership_answer:
        if k not in samplename_dict.keys():
            samplename_dict[k] = cnt
            cnt = cnt + 1
      
    # Set Mixture_answer based on Membership & NP_VAF information
    NUM_CLONE = len(set(membership_answer))
    mixture_answer = np.zeros ((NUM_BLOCK, NUM_CLONE), dtype = 'float')     # Resetting mixture_answer to zero
    for i in range(NUM_BLOCK):
        for j in range(NUM_CLONE):
            mixture_answer[i][j] = round(np.mean(np_vaf[[x  for x in range(len(membership_answer)) if membership_answer[x] == list(samplename_dict.keys())[j]   ]] [: , i] * 2), 3)
    



def main (**kwargs):
    global NUM_BLOCK_INPUT, NUM_BLOCK, RANDOM_PICK, INPUT_TSV, mixture_answer, membership_answer, DEPTH_CUTOFF
    global  df, inputdf, input_containpos, np_vaf, membership, mutation_id, depth_list, NUM_MUTATION, NUM_CLONE, NUM_BLOCK, samplename_dict
    global membership_answer, mutation_id, mixture_answer

    NUM_BLOCK_INPUT = kwargs["NUM_BLOCK_INPUT"] 
    NUM_BLOCK = kwargs["NUM_BLOCK"] 
    RANDOM_PICK = kwargs["RANDOM_PICK"]
    DEPTH_CUTOFF = kwargs["DEPTH_CUTOFF"]
    INPUT_TSV = kwargs["INPUT_TSV"]

    
    makedf(**kwargs)
    
    if kwargs["USE_ALL"]  == False:
        random_pick_fun(**kwargs)
        
    print_record(**kwargs)

            
    return (inputdf, df, np_vaf, membership_answer, mixture_answer,  mutation_id, samplename_dict)
