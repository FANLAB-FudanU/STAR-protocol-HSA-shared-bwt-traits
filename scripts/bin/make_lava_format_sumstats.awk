#!/usr/bin/env awk -f
# Input file should not have header
# Parameters: min_info, min_maf

BEGIN{
    print "SNP\tCHR\tPOS\tA1\tA2\tBETA\tSE\tP\tN";
}

$12>=min_info {
    EAF = $7; 
    if (EAF + 0 < 0.5) {
        maf = EAF + 0
    } else {
        maf = (1 - EAF) + 0
    }

    if (maf >= min_maf) {
        SNP=$1; CHR=$2; POS=$3; A1=$4; A2=$5;
        BETA=$8; SE=$9; P=$10; N=$11;
        printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", SNP, CHR, POS, A1, A2, BETA, SE, P, N);
    }
}



