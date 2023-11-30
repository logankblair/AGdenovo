#selecting blast hits between d mel and outgroup species to make final filtered list of dn genes
library(stringr)
library(dplyr)
library(ggplot2)
library(shadowtext)

fnames <- list.files("/Users/loganblair/Desktop/seq/filtering/stringtie_blast", full.names = T)
fnames[1]
flist=lapply(fnames,function(x) read.table(x,header=FALSE,sep='\t',quote=""))

#read each blast table into list
for(i in 1:length(fnames)){
  flist[[i]]$screened_name=str_split_fixed(string = fnames[i], pattern = "/", 10)[,8]
  
  #remove columns that I don't need from the BLAST table
  flist[[i]][,5:12]=NULL 
}


#read gtf of all transcripts, including unannotated transcripts
unfiltered_gtf=read.table(
  file="/Users/loganblair/Desktop/seq/gtf_files/gtf_21-9-29/filtered/genes_unfiltered_dn_1tpm.gtf", sep='\t')

#add several transcript name columns
unfiltered_gtf$GID=str_split_fixed(unfiltered_gtf$V9,';',9)[,6]
unfiltered_gtf$GID=gsub(" gene_id ","",unfiltered_gtf$GID)
unfiltered_gtf$TID=str_split_fixed(unfiltered_gtf$V9,';',13)[,9]
unfiltered_gtf$TID=gsub(" transcript_id ","",unfiltered_gtf$TID)
blast_df=bind_rows(flist, .id = "column_label")

blast_df=blast_df[blast_df$V3>80,] #require of 80% alignment
blast_df=blast_df[blast_df$V4>50,] #require alignment to be over 50bp
blastt_df1=blast_df
screen_df=data.frame(1:length(unique(unfiltered_gtf$TID)))
screen_df$TID=unfiltered_gtf$TID

for(i in 1:length(fnames)){
  screen_df[,(i+2)]=NA
  names(screen_df)[i+2]=print(str_split_fixed(string = fnames[i], pattern = "/", 10)[,8])
}

rownames(screen_df)=screen_df$TID
screen_df=screen_df[,-1]
screen_df=screen_df[,-1]
fnames1=str_split_fixed(string = fnames, pattern = "/", 10)[,8]

for(i in 1:length(fnames)){
  blast_df_n=blast_df[blast_df["screened_name"]==fnames1[i],]
  screen_df[,i]=unlist(lapply(1:length(rownames(unfiltered_gtf)), function(x) length(which(rownames(screen_df)[x]==blast_df_n$V1))))
}

grep(c('chromosome|intergenic|gene_extended'), colnames(screen_df))
screen_df_genomic=screen_df[,grep(c('chromosome|intergenic|gene_extended'), colnames(screen_df))]

#write.table(screen_df, file = "/Users/loganblair/Desktop/seq/filtering/screen_df1", sep='\t')


screen_df$`ALL_RALL.to.dmel-chromosome-r6.34`=NULL
screen_df$`ALL_RALL.to.dmel-intergenic-r6.34`=NULL
screen_df$`ALL_RALL.to.dsim-all-intergenic-r2.02.fasta`=NULL
screen_df$`ALL_RALL.to.dsim-all-chromosome-r2.02`=NULL
screen_df$`ALL_RALL.to.dsim-all-gene_extended2000-r2.02.fasta`=NULL
screen_df$`ALL_RALL.to.dyak-chromosome-r1.05`=NULL

screen_df=screen_df[,-which(apply(screen_df, 2, sum)==0)]

screen_df$filtercount=apply(screen_df, 1, sum)

colnames(screen_df)=gsub(pattern = "ALL_RALL.to.", replacement = "", colnames(screen_df))
colnames(screen_df)=gsub(pattern = "Stringtie_tpm1.to.", replacement = "1", colnames(screen_df))

colnames(screen_df)
length(which(screen_df$filtercount==0))


####
kallisto_counts=read.table(
  file="/Users/loganblair/Desktop/seq/count_tables/kallisto_lengthscaledtpm",
  sep='\t', quote = "")

rownames(kallisto_counts)=gsub(" ","",rownames(kallisto_counts))
###

screen_df$GID=unfiltered_gtf[match(unfiltered_gtf$TID,rownames(screen_df)),"GID"]

screen_df$GID%in%rownames(kallisto_counts)[kallisto_counts$max_tpm>2]


####################################FILTERING###############################################################
#filtering gameplan
#1) start with transcripts > 1 tpm
#2) get rid of annotated mel things
#3) get rid of transposons
#4) outgroups older than yak
#5) yak 
#6) sim 
kallisto_counts$max_tpm=apply(kallisto_counts, 1, max)
TPM_screen_df=screen_df[screen_df$GID%in%rownames(kallisto_counts)[kallisto_counts$max_tpm>1],]

#reorder columns so important stuff visible
TPM_screen_df <- TPM_screen_df %>%
  dplyr::select(filtercount, everything())
TPM_screen_df <- TPM_screen_df %>%
  dplyr::select(GID, everything())

TPM_screen_df_1=TPM_screen_df %>% group_by(GID) %>% summarise(across(everything(), list(sum))) #collapse transcripts to genes
nrow(TPM_screen_df_1) #1-620

###tpm>1
TPM_screen_df_1=TPM_screen_df_1 %>% mutate_if(is.numeric, ~1 * (. > 0)) #collapse multiple matches into binary 0/1
TPM_screen_df_1=as.data.frame(TPM_screen_df_1)
rownames(TPM_screen_df_1)=TPM_screen_df_1$GID
TPM_screen_df_1$GID=NULL
colnames(TPM_screen_df_1)=gsub("_1", "", colnames(TPM_screen_df_1))


####annotated and formatting
dmel_annotated=c("dmel-3prime-r6.34","dmel-5prime-r6.34","dmel-all-transcript-r6.34.fasta","dmel-CDS-r6.34",
                 "dmel-miRNA-r6.34","dmel-miscRNA-r6.34","dmel-ncRNA-r6.34","dmel-pseudogene-r6.34")
outgroup_annotated_col_list=colnames(TPM_screen_df_1)[56:73]

dmel_annotated_col_list=colnames(TPM_screen_df_1)[colnames(TPM_screen_df_1)%in%dmel_annotated]
colnames(TPM_screen_df_1)[colnames(TPM_screen_df_1)%in%dmel_annotated]

annotated_col_list=c(dmel_annotated_col_list, outgroup_annotated_col_list)

annotated_col_list=annotated_col_list[-which(annotated_col_list=="dsim-intron-r2.02")]
annotated_col_list=annotated_col_list[-which(annotated_col_list=="dyak-intron-r1.05")]

filtered_df_1=TPM_screen_df_1[apply(TPM_screen_df_1[,annotated_col_list], 1, sum)>0,]
TPM_screen_df_2=TPM_screen_df_1[apply(TPM_screen_df_1[,annotated_col_list], 1, sum)==0,] #344

###transposons

TPM_screen_df_3=TPM_screen_df_2[TPM_screen_df_2$`dmel-transposon-r6.34`==0,] #341
filtered_df_2=TPM_screen_df_2[TPM_screen_df_2$`dmel-transposon-r6.34`>0,]

###distant outgroups
outgroup=c("ananassae", "mojavensis", "pseudoobscura", "persimilis", "willistoni","virilis")

outgroup_col_list=str_split_fixed(
  string =(str_split_fixed(string = colnames(TPM_screen_df_1), pattern = "_", 2)[,2]), 
  pattern="\\.",3)[,1]%in%outgroup
colnames(TPM_screen_df_4)[outgroup_col_list]

TPM_screen_df_4=TPM_screen_df_3[apply(TPM_screen_df_3[,outgroup_col_list], 1, sum)<1,] #305
filtered_df_3=TPM_screen_df_3[apply(TPM_screen_df_3[,outgroup_col_list], 1, sum)>0,]

##yak
yak_col_names=colnames(TPM_screen_df)[grep("yak", colnames(TPM_screen_df))]
yak_col_names=append(yak_col_names, colnames(TPM_screen_df)[grep("tai", colnames(TPM_screen_df))])
yak_col_names=yak_col_names[which(yak_col_names%in%"dyak-intron-r1.05"==F)]

TPM_screen_df_5=TPM_screen_df_4[apply(TPM_screen_df_4[,colnames(TPM_screen_df_1)%in%yak_col_names], 1, sum)==0,]#201
filtered_df_4=TPM_screen_df_4[apply(TPM_screen_df_4[,colnames(TPM_screen_df_1)%in%yak_col_names], 1, sum)>0,]

###sim

sim_col_names=colnames(TPM_screen_df)[grep("sim-", colnames(TPM_screen_df))]
sim_col_names=sim_col_names[which(sim_col_names%in%"dsim-intron-r2.02"==F)]
sim_col_names=append(sim_col_names, colnames(TPM_screen_df)[grep("w501", colnames(TPM_screen_df))])
sim_col_names=append(sim_col_names, colnames(TPM_screen_df)[grep("116", colnames(TPM_screen_df))])
sim_col_names=append(sim_col_names, colnames(TPM_screen_df)[grep("LARA", colnames(TPM_screen_df))])
sim_col_names

TPM_screen_df_6=TPM_screen_df_5[apply(TPM_screen_df_5[,colnames(TPM_screen_df_1)%in%sim_col_names], 1, sum)==0,]#129
filtered_df_5=TPM_screen_df_5[apply(TPM_screen_df_5[,colnames(TPM_screen_df_1)%in%sim_col_names], 1, sum)>0,]

final_list=rownames(TPM_screen_df_6)
pre_tpm_synteny_list=as.data.frame(final_list)
write.table(final_list, file='/Users/loganblair/Desktop/seq/gene_tab_files/pre_tpm_synteny_list', 
            col.names = F, row.names = F)


####how many of "distant outgroups" would have been filtered anyway?
colnames(filtered_df_3)
filtered_df_3test=filtered_df_3[,c(2,3,65:77)]
apply(filtered_df_3test, 1, sum)
#not many, but a few

#######making graph of stage-specific filtering counts
filterstep=data.frame(1)
filterstep$tpm_1=nrow(TPM_screen_df_2)
filterstep$transposon=nrow(TPM_screen_df_2)-nrow(filtered_df_2)
filterstep$outgroups=nrow(TPM_screen_df_2)-nrow(filtered_df_2)-nrow(filtered_df_3)
filterstep$yakuba=nrow(TPM_screen_df_2)-nrow(filtered_df_2)-nrow(filtered_df_3)-nrow(filtered_df_4)
filterstep$simulans=nrow(TPM_screen_df_2)-nrow(filtered_df_2)-nrow(filtered_df_3)-nrow(filtered_df_4)-nrow(filtered_df_5)

filterstep=as.data.frame(t(filterstep))
filterstep$name=rownames(filterstep)
filterstep=filterstep[-1,]

ggplot(filterstep, aes(x=name[order(name)], y=V1, group = 1))+
  geom_bar(stat="identity", color='black', fill='grey50')+
  theme_classic()+
  scale_y_continuous(limits = c(0,600), expand = c(0, 0))+
  scale_x_discrete(labels=c("-under 1 tpm","-transposons","-distant outgroups","-D. yakuba","-D. simulans (final)"),)+
  xlab("filtering step")+
  ylab("# transcripts")+
  geom_text(aes(label=V1), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(axis.text = element_text(size=11))
#, axis.text.x = element_text(angle = 10, vjust = 0.95, hjust=1))

##################################genomes##################################
#repeating first part of script in genomes folder instead
#kind of the inverse of the first step - genes that do align are not filtered 
#goal is to figure out how many of dn genes have alignable orthologous sequence in outgroups

fnames_g <- list.files("/Users/loganblair/Desktop/seq/filtering/genome_blast/", full.names = T)
flist_g=lapply(fnames_g,function(x) read.table(x,header=FALSE,sep='\t',quote=""))

#get names and get rid of extraneous columns
for(i in 1:length(fnames_g)){
  flist_g[[i]]$screened_name=str_split_fixed(string = fnames_g[i], pattern = "_", 4)[,4]
  flist_g[[i]][,5:12]=NULL
}
blast_df_g=bind_rows(flist_g, .id = "column_label")

####################################################################################
#super important - setting parameters for identifying orthologous sequence
blast_df_g=blast_df_g[blast_df_g$V3>80,]
blast_df_g=blast_df_g[blast_df_g$V4>200,]

screen_df_g=data.frame(1:length(unique(unfiltered_gtf$TID)))
screen_df_g$TID=unfiltered_gtf$TID

for(i in 1:length(fnames_g)){
  screen_df_g[,(i+2)]=NA
  names(screen_df_g)[i+2]=print(str_split_fixed(string = fnames_g[i], pattern = "_", 4)[,4])
}

#

rownames(screen_df_g)=screen_df_g$TID
screen_df_g=screen_df_g[,-1]
screen_df_g=screen_df_g[,-1]

fnames1_g=str_split_fixed(string = fnames_g, pattern = "_", 4)[,4]

for(i in 1:length(fnames_g)){
  blast_df_g_n=blast_df_g[blast_df_g["screened_name"]==fnames1_g[i],]
  screen_df_g[,i]=unlist(lapply(1:length(rownames(unfiltered_gtf)), function(x) length(which(rownames(screen_df_g)[x]==blast_df_g_n$V1))))
}

screen_df$GID=unfiltered_gtf[match(unfiltered_gtf$TID,rownames(screen_df)),"GID"]

##
screen_df_g$GID=unfiltered_gtf[match(unfiltered_gtf$TID,rownames(screen_df_g)),"GID"]
screen_df_g=screen_df_g %>% group_by(GID) %>% summarise(across(everything(), list(sum))) #collapse transcripts to genes
screen_df_g$nomatch=apply(screen_df_g[,2:ncol(screen_df_g)], 1, sum)
screen_df_g=as.data.frame(screen_df_g)
cantffind=unfiltered_gtf[unfiltered_gtf$GID%in%screen_df_g[which(screen_df_g[screen_df_g$nomatch==0,'GID']%in%final_list==T),'GID'],]

#total counts
length(which(screen_df_g[screen_df_g$yakuba.genome>0,'GID']%in%final_list))
length(which(screen_df_g[screen_df_g$simulans.genome_1>0,'GID']%in%final_list))

####All transcripts are present in at least one genome
final_list%in%cantffind$GID
final_list[which(final_list%in%cantffind$GID==T)] #1 gene not found in other genomes
final_list=final_list[-which(final_list%in%cantffind$GID==T)]
####
unfiltered_gtf_final=unfiltered_gtf[unfiltered_gtf$GID%in%final_list,]
mit_gene=which(unfiltered_gtf_final$V1=='mitochondrion_genome')
mit_gene #no mitochondria gene

final_list=final_list[!(final_list%in%unfiltered_gtf_final[mit_gene,'GID'])]
dn_postmatch_pretpm=as.data.frame(final_list)

#####add Yige's synteny analysis for filtering
UA_synteny=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/UA_synteny_results.txt', header = T)
DN_synteny=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/DN_synteny_results.txt',
                      stringsAsFactors = F, header = T)

dn_postmatch_pretpm1=dn_postmatch_pretpm[dn_postmatch_pretpm$final_list%in%DN_synteny[which(DN_synteny$status=='confirmed'),'Query_gene_ID'],]
#length(dn_postmatch_pretpm$final_list)
#length(dn_postmatch_pretpm1)
#remove 8 genes that failed synteny

#module 2: tpm filtering
#some DN candidates passed match-based filtering and synteny, but have enough short reads mapping to give substantial TPM counts in outgroups
#this script is stolen from the other kalisto TPM script, but instead of making another file I just re run it here
#####################################
#this script takes Kallisto count matrixs (in a single directory) and gtfs as input
#output is scaled count matrix

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager") 
#BiocManager::install("tximportData")
#BiocManager::install("tximport")
#BiocManager::install("rhdf5")
options(scipen=999)

#install.packages("devtools")

#devtools::install_github("pachterlab/sleuth")
#library(tximportData)
library(tximport)
library(sleuth)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)

##
final_screen_list=pre_tpm_synteny_list[,1]
final_screen_list

#####getting table of transcript to gene idx#####

allgtf=read.table(file="/Users/loganblair/Desktop/seq/gtf_files/gtf_21-9-29/6.41_allgtf.gtf", sep='\t')

allgtf_tr=allgtf[grep('RNA',allgtf$V3 ),]
grep('FBgn0287768', allgtf_tr$V9)

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
grep("FBgn0287768", transcript_idx$GENEID)
grep("FBgn0036516", transcript_idx$GENEID)

transcript_idx=unique(transcript_idx)

#######################################################
setwd('/Users/loganblair/Desktop/seq/count_tables/kal_quant/kalisto_quant2/')

abund153=read.table(file="153_all_abundance.tsv", header=T)
alltr=as.data.frame(abund153$target_id)
#write.table(abund153$target_id, file='/Users/loganblair/Desktop/seq/genenames.csv', quote=F, row.names = F)
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

dgrp_list=read.table(file="/Users/loganblair/Desktop/lines.txt", sep='\t')
dgrp_list$V1=gsub("R", "", dgrp_list$V1)

########move to directory with abundance.tsv kalisto output from different libraries
setwd('/Users/loganblair/Desktop/seq/bleh/tissues_dataset/sim_RNAseq/old_samples/')
a=list.files(pattern = "\\.tsv$")
transcript_idx$TXNAME=gsub(' ', '', transcript_idx$TXNAME)
transcript_idx=rbind(transcript_idx, all_tr_CDS)

#IDs come from fasta header of all transcripts.fa
IDs=read.delim(file="/Users/loganblair/Desktop/seq/genomes/col6", sep="\t", header = F)
names(IDs)="id"
IDs$id=gsub("REFSEQ:","", IDs$id)
IDs[,1]=str_split_fixed(IDs[,1], ":", 2)[,2]
IDs[,1]=str_split_fixed(IDs[,1], ",", 2)[,1]
IDs_1=read.delim(file="/Users/loganblair/Desktop/seq/genomes/col6", sep="\t", header = F)
names(IDs_1)="id"
IDs_1$id=gsub("FlyBase_Annotation_IDs:","#", IDs_1$id)
IDs_1[,1]=str_split_fixed(IDs_1[,1], "#", 2)[,2]
IDs[,2]=str_split_fixed(IDs_1[,1], ",", 2)[,1]
Header=read.delim(file="/Users/loganblair/Desktop/seq/genomes/just_headers.tab", header = F)
Header$V1=gsub("parent=","#", Header$V1)
Header$V1=str_split_fixed(Header$V1, "#", 2)[,2]
IDs[,3]=str_split_fixed(Header$V1, ";", 2)[,1]
Header=read.delim(file="/Users/loganblair/Desktop/seq/genomes/just_headers.tab", header = F)
Header$V1=gsub("name=","#", Header$V1)
Header$V1=str_split_fixed(Header$V1, "#", 2)[,2]
IDs[,4]=str_split_fixed(Header$V1, ";", 2)[,1]

txi <- tximport(a, type = "kallisto",
                tx2gene = transcript_idx, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")

lengthscaledcounts=as.data.frame(txi$abundance)
lengthscaledtpm=as.data.frame(txi$abundance)

c(a)
colnames(lengthscaledcounts)=c(a)
rownames(lengthscaledcounts)=gsub(" ", "", rownames(lengthscaledcounts))

##############################################
#going back to the blast filtering dataframe to re check previous filtering stuff
non_mel_col_names=c(yak_col_names,sim_col_names)
outgroup_unannotated_match=(unlist(lapply(1:nrow(TPM_screen_df_4), function(x) sum(TPM_screen_df_4[x,non_mel_col_names]))))
outgroup_unannotated_match
which(outgroup_unannotated_match==0)
which(outgroup_unannotated_match>0)

MSU_match=TPM_screen_df_4[which(outgroup_unannotated_match>0),]
MU_match=TPM_screen_df_4[which(outgroup_unannotated_match==0),]

#test - were the right things filtered?
#(unlist(lapply(1:nrow(MSU_match), function(x) sum(MSU_match[x,outgroup_names]))))
#(unlist(lapply(1:nrow(MU_match), function(x) sum(MU_match[x,outgroup_names]))))

which((unlist(lapply(1:nrow(MSU_match), function(x) sum(MSU_match[x,c(2,3,74,75,76,77)]))))==0)
#4 filtered because of extra yakuba tissues datasets, outside of the male reproductive system

MSU_list1=rownames(TPM_screen_df_4)[which(outgroup_unannotated_match>0)]
#got final melanogaster subgroup, non mel list

##################
#apply TPM screening here!!!

read.delim(file="/Users/loganblair/Desktop/seq/Rplots/paper1_v2/unannotated_outgroup_expression.tsv", sep='\t')
lengthscaledcounts$max=apply(lengthscaledcounts, 1, max)
lengthscaledcounts_MU=lengthscaledcounts[rownames(lengthscaledcounts)%in%final_list,]
nrow(lengthscaledcounts_MU)
#final_list1=rownames(lengthscaledcounts_MU)[which(lengthscaledcounts_MU$max<0.2)]


#####unannotated list required starting here
unann_exp=lengthscaledcounts[which(rownames(lengthscaledcounts)%in%MSU_list1),]
unann_exp$class='UA'
dn_exp=lengthscaledcounts[which(rownames(lengthscaledcounts)%in%pre_tpm_synteny_list$final_list),]
dn_exp$class='DN'
exp_df=rbind(unann_exp, dn_exp)
colnames(exp_df)=gsub('.tsv', "", colnames(exp_df))
exp_df$max=NULL
colnames(exp_df)=c("lara10","116","z_tai18","w501", "class")
exp_df1=exp_df

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
  scale_x_discrete(labels=c('S1', 'S2','S3', 'Y1'))+
  ylab("TPM+0.01")+
  xlab("line")

exp_df1$max=apply(exp_df1[,1:4], 1, max)

length(which(exp_df1[exp_df1$class=='DN', 'max']>1))
length(which(exp_df1[exp_df1$class=='DN', 'max']>0.5))
length(which(exp_df1[exp_df1$class=='DN', 'max']>0.2))
length(which(exp_df1[exp_df1$class=='DN', 'max']<0.1))
length(which(exp_df1[exp_df1$class=='DN', 'max']==0))

test_exp_df1=exp_df1[rownames(exp_df1)%in%final_list1,]
length(which(test_exp_df1[test_exp_df1$class=='DN', 'max']==0))
length(which(test_exp_df1[test_exp_df1$class=='DN', 'max']<.1))


unann_exp$n_zero=apply(unann_exp[,c(1,2,4)], 1, function(x) length(which(x<0.01)))
unann_exp$n_one=apply(unann_exp[,c(1,2,4)], 1, function(x) length(which(x>2)))
unann_exp2=unann_exp[unann_exp$n_zero>0.1,]

write.table(exp_df1, "/Users/loganblair/Desktop/seq/Rplots/paper1_v2/unannotated_outgroup_expression.tsv", sep='\t')

row.names(exp_df1[exp_df1$class=='DN',])


filter1=exp_df1[exp_df1$max>0.2,]

dn_tpm_filter_list=rownames(filter1[filter1$class=='DN',])
dn_tpm_filter_list

#at 0.2 TPM, 15 genes are filtered
length(dn_tpm_filter_list)
######################
#################################
######################

final_list1=dn_postmatch_pretpm1[!dn_postmatch_pretpm1%in%dn_tpm_filter_list]
dn_tpm_filter_list
final_list1
length(final_list1)

write.table(final_list1, file = "/Users/loganblair/Desktop/seq/final_list", sep='\t')
write.table(as.data.frame(MSU_list1), file = "/Users/loganblair/Desktop/seq/final_UA_list", sep='\t')

######################
#################################
######################

#making filtered gene plot

filter_plot=as.data.frame(c("A >1 TPM", "B - TE", "C - distant outgroups", "D - yakuba", "F - simulans"))
filter_plot$transcriptcount=NA

names(filter_plot)=c("filtering step", 'filtered gene count')

length(unique(rownames(TPM_screen_df_1)))
filter_plot[1,2]=length(unique(rownames(TPM_screen_df_2)))
filter_plot[2,2]=length(unique(rownames(TPM_screen_df_3)))
filter_plot[3,2]=length(unique(rownames(TPM_screen_df_4)))
filter_plot[4,2]=length(unique(rownames(TPM_screen_df_5)))
filter_plot[5,2]=length(unique(rownames(TPM_screen_df_6)))

ggplot(filter_plot, aes(x=`filtering step`, y=`filtered gene count`))+
  geom_bar(stat='identity', width=.8, alpha=.6, color='black')+
  scale_x_discrete(labels = c(">1 TPM", "- TE", "- distant outgroups", "- yakuba", "- simulans (final)"))+
  scale_y_continuous(expand = c(0, 0), limits=c(0,450))+
  theme_classic()

outgroup_unannotated=TPM_screen_df_3[row.names(TPM_screen_df_3)%in%row.names(TPM_screen_df_6)==F,]
outgroup_unannotated=row.names(outgroup_unannotated)
write.table(as.data.frame(outgroup_unannotated), file = '/Users/loganblair/Desktop/seq/gene_tab_files/outgroup_unannotated', row.names = F)

screen_df_genomic[,rownames(screen_df_genomic)%in%final_list1]

################
#graphs for filtering
colnames(TPM_screen_df_4)
outgroup_unannotated
final_list1
TPM_screen_df_4=TPM_screen_df_4[rownames(TPM_screen_df_4)%in%c(outgroup_unannotated,final_list1),]
species_graph=data.frame(1:nrow(TPM_screen_df_4))

species_graph=cbind( species_graph, TPM_screen_df_4[,c('trimmed_trinity_sim116','trimmed_trinity_w501', 'LARA10', 'Begunw501')])

#merging yak replicates. Filtering transcripts present in either replicate
species_graph$yak=TPM_screen_df_4$beguntai18+TPM_screen_df_4$trimmed_trinity_tai18
species_graph$yak[which(species_graph$yak>1)]=1

#make column of transcripts overlapping annotated genes

colnames(TPM_screen_df_4)[colnames(TPM_screen_df_4)=='dyak-intron-r1.05']="dya_intron" 
#since I don't want to filter based on introns, changing intron names
colnames(TPM_screen_df_4)[colnames(TPM_screen_df_4)=='dsim-intron-r2.02']="dsi_intron"

#sanity checks, annotated genes should already be filtered
species_graph$ann_sim=apply(TPM_screen_df_4[,grep('dsim', colnames(TPM_screen_df_4))], 1, sum) 
species_graph$ann_sim[which(species_graph$ann_sim>0)]=1 
species_graph$Yangyak=apply(TPM_screen_df_4[,grep('yakuba', colnames(TPM_screen_df_4))], 1, sum) 
species_graph$ann_yak=apply(TPM_screen_df_4[,grep('dyak', colnames(TPM_screen_df_4))], 1, sum) #ditto

species_graph1=species_graph

#checking species graph for "passing" MSU genes
species_graph1$outgrp_sum=apply(species_graph1[,c('trimmed_trinity_sim116', 'trimmed_trinity_w501', 'LARA10', 'yak')], 1, sum)
species_graph1$sim_sum=apply(species_graph1[,c('trimmed_trinity_sim116', 'trimmed_trinity_w501', 'LARA10')], 1, sum)

ggplot(species_graph1, aes(x=outgrp_sum))+
  geom_bar(stat='count')
nrow(final_screen_list)
length(which(species_graph1$outgrp_sum==0))
######
#rarefaction of loss of de novo following addition of simulans/yak transcriptomes

####
pie_df=as.data.frame(c("1mel_denovo", "2mel+poly_sim", "3mel+fixed_sim",  "4mel+yak", "5mel+poly_sim+yak", "6mel+fixed_sim+yak"))
pie_df$counts=NA

pie_df$counts[1]=nrow(TPM_screen_df_6)
pie_df$counts[2]=length(which(species_graph1$sim_sum<3&species_graph1$sim_sum>0&species_graph1$yak==0))
pie_df$counts[3]=length(which(species_graph1$sim_sum==3&species_graph1$yak==0))
pie_df$counts[4]=length(which(species_graph1$sim_sum==0&species_graph1$yak>0))
pie_df$counts[5]=length(which(species_graph1$sim_sum<3&species_graph1$sim_sum>0&species_graph1$yak>0))
pie_df$counts[6]=length(which(species_graph1$sim_sum==3&species_graph1$yak>0))

#pie chart

colnames(pie_df)=c('group','value')
pie_df$percent=NA
for(i in 1:nrow(pie_df)){
  pie_df$percent[i]=pie_df$value[i]/sum(pie_df$value)
}

pie_df$value[i]/sum(pie_df$value)

# Compute the position of labels
pie_df = pie_df %>% 
  dplyr::arrange(dplyr::desc(group)) %>%
  dplyr::mutate(prop = value / sum(value) *100) %>%
  dplyr::mutate(ypos = cumsum(prop)- 0.5*prop )
# Basic piechart
ggplot(pie_df, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, size = 1, color="white") +
  coord_polar("y", start=0) +
  labs(x = NULL, y = NULL, fill = NULL, 
       title = "") +  
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"),
        text=element_text(size=20),
        legend.spacing.y = unit(.2, 'cm')) +
  labs(fill='species')+
  scale_y_reverse()+
  geom_shadowtext(aes(y = ypos, label = paste(round(percent*100, digits=1 ), '%',sep="" )), color = "white", size=6)+
  scale_fill_manual(values=c("#EB6100","#002C60","cornflowerblue",
                             "grey10","grey70", "grey90") ,
                    labels=c("M", 'M+SP', 'M+SF', 'M+Y', 'M+SP+Y', 'M+SF+Y'))+
  guides(fill = guide_legend(reverse = FALSE,byrow = TRUE ))

#make list of transcripts from each class
species_graph1$type=NA

for(i in 1:nrow(species_graph1)) {
  if(species_graph1$sim_sum[i]==0&species_graph1$yak[i]==0){
    species_graph1$type[i]='M'
  }
  if(species_graph1$sim_sum[i]==0&species_graph1$yak[i]>0){
    species_graph1$type[i]='M,Y'
  }
  if(species_graph1$sim_sum[i]>0&species_graph1$sim_sum[i]<3&species_graph1$yak[i]==0){
    species_graph1$type[i]='M,SP'
  }
  if(species_graph1$sim_sum[i]>0&species_graph1$sim_sum[i]<3&species_graph1$yak[i]>0){
    species_graph1$type[i]='M,SP,Y'
  }
  if(species_graph1$sim_sum[i]==3&species_graph1$yak[i]==0){
    species_graph1$type[i]='M,SF'
  }
  if(species_graph1$sim_sum[i]==3&species_graph1$yak[i]>0){
    species_graph1$type[i]='M,SF,Y'
  }
} 

gene_list_MY=row.names(species_graph1)[species_graph1$type=='M,Y']
gene_list_MSP=row.names(species_graph1)[species_graph1$type=='M,SP']
gene_list_MSPY=row.names(species_graph1)[species_graph1$type=='M,SP,Y']
gene_list_MSF=row.names(species_graph1)[species_graph1$type=='M,SF']
gene_list_MSFY=row.names(species_graph1)[species_graph1$type=='M,SF,Y']

#123 total MSU genes from adding up all the classes
length(gene_list_MSP)+length(gene_list_MSPY)+length(gene_list_MSF)+length(gene_list_MSFY)+length(gene_list_MY)

#other species specific gene lists
gene_list_simfixed=row.names(species_graph1)[which(species_graph1$type=='M,SF,Y'|species_graph1$type=='M,SF')]
gene_list_threespecies=row.names(species_graph1)[which(species_graph1$type=='M,SF,Y'|species_graph1$type=='M,SP,Y')]
gene_list_sim1=row.names(species_graph1)[species_graph1$sim_sum==1]
gene_list_sim2=row.names(species_graph1)[species_graph1$sim_sum==2]
gene_list_sim3=row.names(species_graph1)[species_graph1$sim_sum==3]

#final MSU list
unnan_list_nogenome=row.names(species_graph1)[!species_graph1$type=='M']
unnan_list_nogenome1=as.data.frame(unnan_list_nogenome)
write.table(unnan_list_nogenome1, file='/Users/loganblair/Desktop/seq/gene_tab_files/unnan_list_nogenome', row.names = F, col.names = F)

