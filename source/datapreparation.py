import numpy as np
import pandas as pd 
import random

def makedf ( **kwargs ):
    global input_containpos, mutation_id, membership, membership_answer, inputdf, df,  np_vaf, np_BQ, samplename_dict_CharacterToNum, samplename_dict_NumToCharacter

    input_containpos = pd.read_csv( kwargs["INPUT_TSV"],  header = None, sep = "\t") 
    
    if input_containpos.shape[1] == 3: #  3번쨰 column (BQ)가 없다면
        input_containpos.columns = ["pos", "sample", "info"]
        input_containpos.astype ({"pos":"str", "sample":"str", "info":"str"})
    elif input_containpos.shape[1] == 4: #  3번쨰 column (BQ)가 있다면
        input_containpos.columns = ["pos", "sample", "info", "BQ"]
        input_containpos.astype ({"pos":"str", "sample":"str", "info":"str", "BQ":"str"})
    
    
    input_containpos ["cha1"] = "child"  # monoclone이면 child, 둘 이상의 clone이 합쳐진거면 parent
    input_containpos ["cha2"] = "space"       # 축 상에 있으면 axis, 공간 상에 있으면 space
    samplename_dict_CharacterToNum = {}
    samplename_dict_NumToCharacter = {}

    if kwargs["NUM_MUTATION"] == -1:
        kwargs["NUM_MUTATION"] = input_containpos.shape[0]
        kwargs["RANDOM_PICK"] = input_containpos.shape[0]

    np_vaf = np.zeros(( kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]), dtype = 'float')
    np_BQ = np.zeros(( kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]), dtype = 'float')       # BQ를 담은 것. 초기값으로 20으로 setting해준다
    np_BQ.fill(20)
    inputdf = pd.DataFrame (np.zeros(( kwargs["NUM_MUTATION"], kwargs["NUM_BLOCK"]), dtype = 'object'), columns = ['block' + str(i + 1) for i in range(kwargs["NUM_BLOCK"])])
    mutation_id = []
    membership = []
    depth_list = []
    


    # input 형식은 n * 3 으로 구성 :   ID (chr_pos), membmership(정답 set 일 경우),  NUM_BLOCK_INPUT(3)만큼의 depth, alt 정보

    depth_col = [[]] * int(len(input_containpos.iloc[0][2].split(","))/2)
    depth_row = []
    for row in range( kwargs["NUM_MUTATION"] ):
        depth_row_mini = []
        mutation_id.append( str(input_containpos.iloc[row][0]) )            # "pos"
        membership.append( str(input_containpos.iloc[row][1]) )           # "sample"
        if "," in str(input_containpos.iloc[row][1]) :
            input_containpos.loc[row,"cha1"] = "parent"

        if str(input_containpos.iloc[row][1]) not in samplename_dict_CharacterToNum.keys():
            samplename_dict_CharacterToNum[str(input_containpos.iloc[row][1])] = int (len(samplename_dict_CharacterToNum))      # {'other': 0, 'V5': 1, 'V3': 2, 'V1': 3}           # 각 sample name을 숫자화시킴

        #rmv_bracket = re.sub("[\[\] ]", '', str(input_containpos.iloc[row][2])).split(",")            # [194, 25, 193, 66, 0, 0] 라고 되어 있는데 bracket과 한 칸 공백을 지움
        rmv_bracket=input_containpos.iloc[row][2].split(",")     # 2번째 column의 예시  112,27,104,9  
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

            if "BQ" in input_containpos.columns: #  3번쨰 column (BQ)가 있다면
                if kwargs["NUM_BLOCK"] == 1:
                    BQ_input = [input_containpos.iloc[row][3]]     # 1D는 ","가 없으니까 splitdㅣ 에러난다
                else:
                    BQ_input = input_containpos.iloc[row][3].split(",")     # 3번째 column (BQ)의 예시  20,20
                for i in range (0, len(BQ_input) ):
                    np_BQ [row][i] = int ( BQ_input[i] )
                    if int ( BQ_input[i] ) == 0:       # AXIS mutation의 경우 0,0이고 BQ 0일 테니까 default 20으로 넣어준다
                        np_BQ [row][i] = 20 

        depth_row.append (depth_row_mini)

    # "0.0.0"을 그대로 놔둘 수 없다.  평균 depth로 갈음해서 바꿔 넣는다  (alt는 0으로 유지)

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

    # RANDOM하게 n개만 뽑기
    random.seed(kwargs["RANDOM_SEED"])
    random_index = sorted( list ( range (0, kwargs["NUM_MUTATION"])  ))  
    
    input_containpos =  input_containpos.iloc[random_index]
    inputdf  = inputdf.iloc[random_index]
    df = [df[i] for i in random_index]
    np_vaf = np_vaf [random_index]
    membership_answer = [membership[i] for i in random_index]
    mutation_id = [mutation_id[i] for i in random_index]

    #np_vaf + membership를 df 형식으로 하고 RANDOM_PICK개만 출력 
    t = pd.DataFrame(np_vaf, columns = ["block{0}".format(i) for i in range( kwargs["NUM_BLOCK"])], index = mutation_id)
    t["membership_answer"] = pd.Series(membership_answer, index = mutation_id)
    t.to_csv ("{0}/npvaf.txt".format( kwargs["CLEMENT_DIR"] ), index = True, header=True, sep = "\t")


    samplename_dict_CharacterToNum, cnt = {}, 0
    for k in membership_answer:
        if k not in samplename_dict_CharacterToNum.keys():
            samplename_dict_CharacterToNum[k] = cnt
            cnt = cnt + 1
    

    return kwargs




def main (**kwargs):
    global input_containpos, mutation_id, membership, membership_answer, inputdf, df,  np_vaf, np_BQ, samplename_dict_CharacterToNum, samplename_dict_NumToCharacter
    
    kwargs = makedf ( **kwargs )
    
    kwargs = random_pick_fun(**kwargs)

    return (inputdf, df, np_vaf, np_BQ, membership_answer, mutation_id, samplename_dict_CharacterToNum, kwargs )
