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
            NUM_BLOCK = int(len(line.split("\t")[-1].split(",")) / 2)
        return INPUT_FILETYPE, NUM_BLOCK
    