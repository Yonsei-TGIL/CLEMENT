library ("fclust")

args <- commandArgs(trailingOnly = TRUE)
TSV1 = args[1]
TSV2 = args[2]
OUTPUT_ARI = args[3]
transform  = args[4]

# print (TSV1)
# print (TSV2)
# print (OUTPUT_ARI)




ANSWER_MEMBERSHIP = as.factor ( scan( TSV1, sep="\t", quiet=TRUE ) ) 

FUZZY_MATRIX =  read.table( TSV2, header = TRUE, row.names = 1, sep = "\t")
colnames(FUZZY_MATRIX) <- gsub("^X", "", colnames(FUZZY_MATRIX))

# 적용할 함수 정의 (10을 밑으로 하는 지수 연산)
exponential_function <- function(x) {
  return(10^x)
}

# sapply를 사용하여 함수를 각 열에 적용하여 새로운 데이터프레임 생성
if (transform == "False"){
    FUZZY_MATRIX_10 = as.data.frame(sapply(FUZZY_MATRIX, exponential_function))
    FUZZY_MATRIX  = FUZZY_MATRIX_10
}

# print (ANSWER_MEMBERSHIP [0:20])
# print (FUZZY_MATRIX [0:20, ])

ari.f = ARI.F(VC = ANSWER_MEMBERSHIP, U = FUZZY_MATRIX)

writeLines( as.character(ari.f), OUTPUT_ARI)