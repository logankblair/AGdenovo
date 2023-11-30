library(rtracklayer)
library(GenomicFeatures)

a=import.gff("/Users/loganblair/Desktop/seq/gtf_files/gtf_21-9-29/gtf_files/1tpm_29asmbl/refcom_29_1tpm/unfiltered_dn_1tpm.gtf")
b=import.gff("/Users/loganblair/Desktop/seq/gtf_files/dmel-all-r6.34.gtf")

export(c(a,b),"/Users/loganblair/Desktop/seq/gtf_files/allgenes_MU_MSU_combined.gtf")

