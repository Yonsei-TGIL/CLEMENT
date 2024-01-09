def main (INPUT_TSV):
    with open (INPUT_TSV, "r", encoding = "utf8") as input_file:
        line = input_file.readlines()[0].rstrip("\n")
        if "VCF" in line:
            INPUT_FILETYPE = "VCF"
            while True:
                line = input_file.readline().rstrip("\n")
                if "#CHROM" in line:
                    NUM_BLOCK = len(line.split("\t")) - 9
                    break
        else:
            INPUT_FILETYPE = "TSV"
            if len(line.split("\t")) == 4:  # If 4th column (BQ) is present
                NUM_BLOCK = int(len(line.split("\t")[-2].split(",")) / 2)         
            elif len(line.split("\t")) == 3:  # If 4th column (BQ) is absent
                NUM_BLOCK = int(len(line.split("\t")[-1].split(",")) / 2)
        return INPUT_FILETYPE, NUM_BLOCK

def sexdetermination (INPUT_TSV):
    import pandas as pd
    inputdf = pd.read_csv(INPUT_TSV, sep = "\t")
    contains_y = bool ( inputdf.iloc[:, 0].str.contains('Y').any() )
    
    if contains_y == True:
        return "M"
    else:
        return "F"

    
