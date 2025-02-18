#this script takes Kallisto count matrices (in a single directory) and gtfs as input
#output is scaled count matrix

setwd("~/Desktop/github1/")

library(tximport)
library(sleuth)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)

#################################################
#####getting table of transcript to gene idx#####
#################################################
#6.41 dmel appended to TACO output
allgtf=read.table(file="6.41_all.gtf", sep='\t')

allgtf_tr=allgtf[grep('RNA',allgtf$V3 ),]
allgtf_tr$gene_id=str_split_fixed(allgtf_tr$V9, ";", 3)[,1]

allgtf_tr$transcript_id=str_split_fixed(allgtf_tr$V9, ";", 4)[,c(-1,-2,-4,-5)]
allgtf_tr$gene_id=gsub("gene_id ","",allgtf_tr$gene_id)
allgtf_tr$transcript_id=gsub("transcript_id ","",allgtf_tr$transcript_id)
transcript_idx=data.frame(allgtf_tr$transcript_id,allgtf_tr$gene_id)

#add de novo names
allgtf_dn=allgtf[grep('assembly',allgtf$V2 ),]
allgtf_dn=allgtf_dn[allgtf_dn$V3=='transcript',]
allgtf_dn$gene_id=sapply(strsplit(gsub('gene_id ', '', allgtf_dn$V9), ';'),`[`, 6)
allgtf_dn$transcript_id=sapply(strsplit(gsub('transcript_id ', '', allgtf_dn$V9), ';'),`[`, 9)
transcript_idx1=data.frame(allgtf_dn$transcript_id,allgtf_dn$gene_id)

#combine
colnames(transcript_idx)=c('TXNAME', 'GENEID')
colnames(transcript_idx1)=c('TXNAME', 'GENEID')
transcript_idx=rbind(transcript_idx, transcript_idx1)
transcript_idx=unique(transcript_idx)

#######################################################

abund153=read.table(file="kallisto_tables/DGRP/153_abundance.tsv", header=T)
alltr=as.data.frame(abund153$target_id)
colnames(alltr)="transcript_name"

all_tr_dn=as.data.frame(alltr[grep("^TU", alltr$transcript_name),])
all_tr_RNA=as.data.frame(alltr[grep("^FBtr", alltr$transcript_name),])
all_tr_CDS=as.data.frame(alltr[grep("-P", alltr$transcript_name),])
all_tr_CDS[-(grep(("^TU|^FBtr|-P"), alltr$transcript_name)),]

colnames(all_tr_CDS)='TXNAME'
all_tr_CDS$GENEID=NA
all_tr_CDS$TXNAME=gsub("-P", "#", all_tr_CDS$TXNAME)
all_tr_CDS$GENEID=str_split_fixed(all_tr_CDS$TXNAME, "#",2)[,1]
all_tr_CDS$TXNAME=gsub("#", "-P", all_tr_CDS$TXNAME)

dgrp_list=read.table(file="lines.txt", sep='\t')
dgrp_list$V1=gsub("R", "", dgrp_list$V1)


a=list.files(path = "kallisto_tables/DGRP/", pattern = "tsv$")
transcript_idx$TXNAME=gsub(' ', '', transcript_idx$TXNAME)
transcript_idx=rbind(transcript_idx, all_tr_CDS)

###########
#IDs come from fasta header of all transcripts.fa, from flybase
IDs=read.delim(file="tables/flybase_transcripts.facol6.txt", sep="\t", header = F)
names(IDs)="id"
IDs$id=gsub("REFSEQ:","", IDs$id)
IDs[,1]=str_split_fixed(IDs[,1], ":", 2)[,2]
IDs[,1]=str_split_fixed(IDs[,1], ",", 2)[,1]
IDs_1=read.delim(file="tables/flybase_transcripts.facol6.txt", sep="\t", header = F)
names(IDs_1)="id"
IDs_1$id=gsub("FlyBase_Annotation_IDs:","#", IDs_1$id)
IDs_1[,1]=str_split_fixed(IDs_1[,1], "#", 2)[,2]
IDs[,2]=str_split_fixed(IDs_1[,1], ",", 2)[,1]
Header=read.delim(file="tables/flybase_transcripts_just_headers.tab", header = F)
Header$V1=gsub("parent=","#", Header$V1)
Header$V1=str_split_fixed(Header$V1, "#", 2)[,2]
IDs[,3]=str_split_fixed(Header$V1, ";", 2)[,1]
Header=read.delim(file="tables/flybase_transcripts_just_headers.tab", header = F)
Header$V1=gsub("name=","#", Header$V1)
Header$V1=str_split_fixed(Header$V1, "#", 2)[,2]
IDs[,4]=str_split_fixed(Header$V1, ";", 2)[,1]





txi <- tximport(paste0('kallisto_tables/DGRP/',a), type = "kallisto",
                tx2gene = transcript_idx, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")
lengthscaledcounts=as.data.frame(txi$abundance)
#change IDS

#################
transcript_idx1=as.data.frame(c(IDs$id,IDs$V2,IDs$V4))
names(transcript_idx1)="TXNAME"
transcript_idx1$GENEID=c(IDs$V3,IDs$V3,IDs$V3)
rownames(lengthscaledcounts)=gsub(" ", "", rownames(lengthscaledcounts))

##for three dashes
extra_char_row1=as.data.frame(str_split_fixed(IDs$V4,"-",  4)[,4])
extra_char_loc1=grep("R", extra_char_row1[,1])

extra_char_df1=IDs[extra_char_loc1,]
extra_char_df_columns1=as.data.frame(str_split_fixed(extra_char_df1$V4, '-',4)[,1:3])
extra_char_df1$V4=paste(extra_char_df_columns1[,1],extra_char_df_columns1[,2], 
                        extra_char_df_columns1[,3],sep='-')

#IDs$V4=str_split_fixed(IDs$V4, '-', 3)[,1]
#for two dashes
extra_char_row=as.data.frame(str_split_fixed(IDs$V4,"-",  3)[,3])
extra_char_loc=grep("R", extra_char_row[,1])

extra_char_df=IDs[extra_char_loc,]
extra_char_df_columns=as.data.frame(str_split_fixed(extra_char_df$V4, '-',3)[,1:2])
extra_char_df$V4=paste(extra_char_df_columns[,1],extra_char_df_columns[,2],sep='-')
extra_char_df$V4
extra_char_df1$V4
IDs$V4=str_split_fixed(IDs$V4, '-', 2)[,1]
IDs[extra_char_loc,4]=extra_char_df$V4
IDs[extra_char_loc1,4]=extra_char_df1$V4

no_id=which(IDs$V4%in%row.names(lengthscaledcounts)==F)

####some just weren't matching so I manualy looked them up
IDs$V4[no_id]
no_id1=as.data.frame(row.names(RNA_counts)[which(row.names(RNA_counts)%in%IDs$V4==F)])


rownames(lengthscaledcounts)=gsub("beta-PheRS","FBgn0039175",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("PIG-P","FBgn0039405",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("alpha-PheRS","FBgn0030007",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("mei-P26","FBgn0026206",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("GlcAT-P","FBgn0036144",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("pre-mod(mdg4)-P","FBgn0266175",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("Est-P","FBgn0000594",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("mei-P22","FBgn0016036",rownames(lengthscaledcounts))
rownames(lengthscaledcounts)=gsub("ND-PDSW","FBgn0021967",rownames(lengthscaledcounts))



IDs=IDs[,3:4]
IDs=IDs[!duplicated(IDs[,1]),]
####replacing gene names with fbgn id
for(i in 1:nrow(lengthscaledcounts)){
  i_row=rownames(lengthscaledcounts)[i]
  i_row=IDs[IDs$V4%in%i_row,1]
  if(length(i_row)>0){
    rownames(lengthscaledcounts)[i]=i_row
  }
}

colnames(lengthscaledcounts)=c('R153','R217','R229','R287','R304','R320','R338','R352','R357','R359','R360','R370','R380','R399','R517','R530','R563','R630','R703','R761','R805','R812','R822','R85','R850','R88','R900','R911','R93','ed10','iso1')
rownames(lengthscaledcounts) <- gsub(" ", "", rownames(lengthscaledcounts))
write.table(lengthscaledcounts, "kallisto_expression.tsv", quote=F, sep='\t')

