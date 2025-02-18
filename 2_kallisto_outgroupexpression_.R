#this script takes Kallisto count matrices (in a single directory) and gtfs as input
#output is scaled count matrix

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager") 
#BiocManager::install("tximportData")
#BiocManager::install("tximport")
#BiocManager::install("rhdf5")
options(scipen=999)

setwd("~/Desktop/github1/")

library(tximport)
library(sleuth)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)

##
final_screen_list=read.delim(file='MOD_list.txt',
                             sep='\t',stringsAsFactors = F, header=T)
unan_out_list=unan_out_list=read.table('MSD_list.txt', sep='\t')



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


a=list.files(path = "kallisto_tables/outgroups", pattern = "tsv$")
a
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
Header=read.delim(file="flybase_transcripts_just_headers.tab", header = F)
Header$V1=gsub("name=","#", Header$V1)
Header$V1=str_split_fixed(Header$V1, "#", 2)[,2]
IDs[,4]=str_split_fixed(Header$V1, ";", 2)[,1]

transcript_idx1$FBGN
transcript_idx1$TXNAME
transcript_idx1$GENEID

txi <- tximport(paste0('kallisto_tables/outgroups/',a), type = "kallisto",
                tx2gene = transcript_idx, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")


#############################################
lengthscaledcounts=as.data.frame(txi$abundance)

#############################################

lengthscaledtpm=as.data.frame(txi$abundance)

#colnames(lengthscaledcounts)=c('R153','R217','R229','R287','R304','R320','R338','R352','R357','R359','R360','R370','R380','R399','R517','R530','R563','R630','R703','R761','R805','R812','R822','R85','R850','R88','R900','R911','R93')
c(a)
colnames(lengthscaledcounts)=c(a)

#lengthscaledcounts$R517_mean=apply(lengthscaledcounts[,1:3], 1, mean)
#lengthscaledcounts$w501_mean=apply(lengthscaledcounts[,4:6], 1, mean)
rownames(lengthscaledcounts)=gsub(" ", "", rownames(lengthscaledcounts))

unann_exp=lengthscaledcounts[which(rownames(lengthscaledcounts)%in%unan_out_list[,1]),]
unann_exp$class='UA'
dn_exp=lengthscaledcounts[which(rownames(lengthscaledcounts)%in%final_screen_list[,1]),]
dn_exp$class='DN'
exp_df=rbind(unann_exp, dn_exp)
exp_df1=exp_df
colnames(exp_df)=gsub('.tsv', "", colnames(exp_df))
colnames(exp_df)=c("lara10","116","z_tai18","w501", "class")

exp_df=pivot_longer(exp_df,
  cols = 'lara10':'w501', 
  names_to = "type",
  values_to = "value"
)

exp_df$species='sim'
exp_df$species[exp_df$type=='z_tai18']='yak'

ggplot(exp_df, aes(x=type, y=value+0.01, fill=species))+
  geom_boxplot()+
  scale_y_log10()+
  facet_wrap(~class)+
  scale_fill_manual(values = c('grey50', 'white'))+
  theme_classic()+
  theme(text=element_text(size=25),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(labels=c('116', 'lara10','w501', 'tai18'))+
  ylab("TPM + 0.01")+
  xlab("line")

ggplot(exp_df[exp_df$class=='DN',], aes(x=type, y=value+0.01, fill=species))+
  geom_violin()+
  #scale_y_log10()+
  scale_fill_manual(values = c('grey50', 'white'))+
  theme_classic()+
  theme(text=element_text(size=25),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(labels=c('S1', 'S2','S3', 'Y1'))+
  ylab("TPM")+
  xlab("line")+
  geom_hline(yintercept = 0.2,linetype='dotted')

ggplot(exp_df[exp_df$class=='UA',], aes(x=type, y=value+0.01, fill=species))+
  geom_boxplot()+  
  scale_y_log10()+
  scale_fill_manual(values = c('grey50', 'white'))+
  theme_classic()+
  theme(text=element_text(size=25),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(labels=c('SZ116', 'lara10','w501', 'tai18'))+
  ylab("MSD TPM + 0.01")+
  xlab("line")

exp_df1$max=apply(exp_df1[,1:4], 1, max)

write.table(exp_df1, "tables/unannotated_outgroup_expression.tsv", sep='\t')

#testing Peng et al transcripts
candidate_AG_DN=c('FBgn0053244','FBgn0053245','FBgn0259973',
                  'FBgn0261006','FBgn0263090','FBgn0263095','FBgn0265424',
                  'FBgn0265739','FBgn0265798','FBgn0265842','FBgn0265900',
                  'FBgn0267179','FBgn0267215','FBgn0267268','FBgn0267926')

candidate_AG_DN_exp=lengthscaledcounts[rownames(lengthscaledcounts)%in%candidate_AG_DN,]

candidate_AG_DN_exp$max=apply(candidate_AG_DN_exp, 1, max)
row.names(candidate_AG_DN_exp)[which(candidate_AG_DN_exp$max<0.2)]



