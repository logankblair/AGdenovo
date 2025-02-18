#selecting blast hits between d mel and outgroup species to make final filtered list of dn genes
library(stringr)
library(dplyr)
library(ggplot2)
library(shadowtext)

setwd("~/Desktop/github1/")
fnames <- list.files("filtering/stringtie_blast", full.names = T)
flist=lapply(fnames,function(x) read.table(x,header=FALSE,sep='\t',quote=""))

#read each blast table into list

nslash <- str_count(string= fnames[1], pattern = "/")

for(i in 1:length(fnames)){
  flist[[i]]$screened_name=str_split_fixed(string = fnames[i], pattern = "/", 10)[,(nslash+1)]
  
  #remove columns that I don't need from the BLAST table
  flist[[i]][,5:12]=NULL 
}


#read gtf of all transcripts, including unannotated transcripts
unfiltered_gtf=read.table(
  file="taco_unfiltered_1tpm.gtf", sep='\t')

#add several transcript name columns
unfiltered_gtf$GID=str_split_fixed(unfiltered_gtf$V9,';',9)[,6]
unfiltered_gtf$GID=gsub(" gene_id ","",unfiltered_gtf$GID)
unfiltered_gtf$TID=str_split_fixed(unfiltered_gtf$V9,';',13)[,9]
unfiltered_gtf$TID=gsub(" transcript_id ","",unfiltered_gtf$TID)
blast_df=bind_rows(flist, .id = "column_label")

blast_df=blast_df[blast_df$V3>80,] #require of 80% alignment
blast_df=blast_df[blast_df$V4>50,] #require alignment to be over 50bp
screen_df=data.frame(1:length(unique(unfiltered_gtf$TID)))
screen_df$TID=unfiltered_gtf$TID

#fill out table of hits (genes - rows) to outgroups (columns - outgroup dataset)
for(i in 1:length(fnames)){
  screen_df[,(i+2)]=NA
  names(screen_df)[i+2]=(str_split_fixed(string = fnames[i], pattern = "/", 10)[,(nslash+1)])
}

rownames(screen_df)=screen_df$TID
screen_df$X1.length.unique.unfiltered_gtf.TID..=NULL
screen_df$TID=NULL
fnames1=str_split_fixed(string = fnames, pattern = "/", 10)[,(nslash+1)]

for(i in 1:length(fnames)){
  blast_df_n=blast_df[blast_df["screened_name"]==fnames1[i],]
  screen_df[,i]=unlist(lapply(1:length(rownames(unfiltered_gtf)), function(x) length(which(rownames(screen_df)[x]==blast_df_n$V1))))
}

# I now have a dataframe with dimensions of A rows and B columns
#A = number of >1 TPM genes
#B = number of outgroup transcriptomes to filter against

# If there is a BLAST match to that particular outgroup transcriptome, that cell will be >0

# each filtering step looks for rows with cells >0, under grouped columns
# for instance, all columns corresponding to D. simulans transcriptomes will be a single step

#--------------------------------------------------------------------------------

#here, I separate out the matches to genomic regions and not the transcriptomes
#these are treated in the opposite manner - 
#genes need to have matches here, otherwise they will fail the orthology criteria
screen_df_genomic=screen_df[,grep(c('chromosome|intergenic|gene_extended'), colnames(screen_df))]

omit_cols <- c('ALL_RALL.to.dmel-chromosome-r6.34',
               'ALL_RALL.to.dmel-intergenic-r6.34',
               'ALL_RALL.to.dsim-all-intergenic-r2.02.fasta',
               'ALL_RALL.to.dsim-all-chromosome-r2.02',
               'ALL_RALL.to.dsim-all-gene_extended2000-r2.02.fasta',
               'ALL_RALL.to.dyak-chromosome-r1.05')

screen_df <- screen_df[,!colnames(screen_df)%in%omit_cols]

#no point in using columns with no matches
screen_df <- screen_df[,-which(apply(screen_df, 2, sum)==0)]

#sum of filtered cols per gene
screen_df$filtercount=apply(screen_df, 1, sum)

#cleaning up colnames
colnames(screen_df)=gsub(pattern = "ALL_RALL.to.", replacement = "", colnames(screen_df))
colnames(screen_df)=gsub(pattern = "Stringtie_tpm1.to.", replacement = "1", colnames(screen_df))

#getting TPM counts from kallisto table here
kallisto_counts=read.table(file="kallisto_expression.tsv", sep='\t', quote = "")
rownames(kallisto_counts)=gsub(" ","",rownames(kallisto_counts))
screen_df$GID=unfiltered_gtf[match(unfiltered_gtf$TID,rownames(screen_df)),"GID"]

####################################FILTERING###############################################################

#filtering gameplan
#1) start with transcripts > 1 tpm
#2) get rid of annotated mel things
#3) get rid of transposons
#4) outgroups older than yak
#5) yak 
#6) sim 

kallisto_counts$max_tpm=apply(kallisto_counts, 1, max)

#subset to genes >1TPM
TPM_screen_df=screen_df[screen_df$GID%in%rownames(kallisto_counts)[kallisto_counts$max_tpm>1],]

FILTER_1_not_enough_exp <- data.frame(filteredgenes=screen_df[!screen_df$GID%in%TPM_screen_df$GID,"GID"])

#reorder columns so important stuff visible
TPM_screen_df <- TPM_screen_df %>%
  dplyr::select(filtercount, everything())
TPM_screen_df <- TPM_screen_df %>%
  dplyr::select(GID, everything())

Passed_filter1=TPM_screen_df %>% group_by(GID) %>% summarise(across(everything(), list(sum))) #collapse transcripts to genes
nrow(Passed_filter1) #1-620
#Passed_filter1$`dyak-intron-r1.05`=NULL
#Passed_filter1$`dsim-intron-r1.05`=NULL

###tpm>1
Passed_filter1=Passed_filter1 %>% mutate_if(is.numeric, ~1 * (. > 0)) #collapse multiple matches into binary 0/1
Passed_filter1=as.data.frame(Passed_filter1)
rownames(Passed_filter1)=Passed_filter1$GID
Passed_filter1$GID=NULL
colnames(Passed_filter1)=gsub("_1", "", colnames(Passed_filter1))


####annotated and formatting
dmel_annotated=c("dmel-3prime-r6.34","dmel-5prime-r6.34","dmel-all-transcript-r6.34.fasta","dmel-CDS-r6.34",
                 "dmel-miRNA-r6.34","dmel-miscRNA-r6.34","dmel-ncRNA-r6.34","dmel-pseudogene-r6.34")

annotated_col_list=colnames(Passed_filter1)[colnames(Passed_filter1)%in%dmel_annotated]
annotated_col_list=annotated_col_list[-which(annotated_col_list=="dsim-intron-r2.02")]
annotated_col_list=annotated_col_list[-which(annotated_col_list=="dyak-intron-r1.05")]

Passed_filter2=Passed_filter1[apply(Passed_filter1[,annotated_col_list], 1, sum)==0,] #344


###transposons

Passed_filter3=Passed_filter2[Passed_filter2$`dmel-transposon-r6.34`==0,] #341
FILTER_2_transposons=Passed_filter2[Passed_filter2$`dmel-transposon-r6.34`>0,]

###distant outgroups
outgroup=c("ananassae", "mojavensis", "pseudoobscura", "persimilis", "willistoni","virilis")

outgroup_col_list=str_split_fixed(
  string =(str_split_fixed(string = colnames(Passed_filter1), pattern = "_", 2)[,2]), 
  pattern="\\.",3)[,1]%in%outgroup

Passed_filter4=Passed_filter3[apply(Passed_filter3[,outgroup_col_list], 1, sum)<1,] #305

FILTER_3_distant_outgroups=Passed_filter3[apply(Passed_filter3[,outgroup_col_list], 1, sum)>0,]

##yak
yak_col_names=colnames(TPM_screen_df)[grep("yak", colnames(TPM_screen_df))]
yak_col_names=append(yak_col_names, colnames(TPM_screen_df)[grep("tai", colnames(TPM_screen_df))])
yak_col_names=yak_col_names[which(yak_col_names%in%"dyak-intron-r1.05"==F)]
yak_col_names_notranscipt=yak_col_names[15:21]

test_Passed_filter4=Passed_filter4[,colnames(Passed_filter4)%in%yak_col_names_notranscipt]

remove_for_UA_yak<-names(which(apply(test_Passed_filter4, 1, mean)>0))

print(paste(length(which(apply(test_Passed_filter4, 1, mean)>0)), "yak transcripts filtered from annotations DNA"))

Passed_filter5=Passed_filter4[apply(Passed_filter4[,colnames(Passed_filter1)%in%yak_col_names], 1, sum)==0,]#201
FILTER_4_yakuba=Passed_filter4[apply(Passed_filter4[,colnames(Passed_filter1)%in%yak_col_names], 1, sum)>0,]

##sim
sim_col_names=colnames(TPM_screen_df)[grep("sim-", colnames(TPM_screen_df))]
sim_col_names=sim_col_names[which(sim_col_names%in%"dsim-intron-r2.02"==F)]
sim_col_names=append(sim_col_names, colnames(TPM_screen_df)[grep("w501", colnames(TPM_screen_df))])
sim_col_names=append(sim_col_names, colnames(TPM_screen_df)[grep("116", colnames(TPM_screen_df))])
sim_col_names=append(sim_col_names, colnames(TPM_screen_df)[grep("LARA", colnames(TPM_screen_df))])

sim_col_names_genomic=sim_col_names[1:9]
TPM_screen_df_6=Passed_filter5[apply(Passed_filter5[,colnames(Passed_filter1)%in%sim_col_names], 1, sum)==0,]#129
FILTER_5_sim=Passed_filter5[apply(Passed_filter5[,colnames(Passed_filter1)%in%sim_col_names], 1, sum)>0,]

Passed_filter5[apply(Passed_filter5[,colnames(Passed_filter1)%in%sim_col_names], 1, sum)>0,]

test_Passed_filter5=Passed_filter5[,colnames(Passed_filter5)%in%sim_col_names_genomic]
length(which(apply(test_Passed_filter5, 1, mean)>0))#23 yak transcripts filtered from annotations DNA
sim_ann=row.names(test_Passed_filter5)[(which(apply(test_Passed_filter5, 1, mean)>0))]
remove_for_UA_sim<-sim_ann

initial_list=rownames(TPM_screen_df_6)
pre_tpm_synteny_list=as.data.frame(initial_list)
#write.table(initial_list, file='/Users/loganblair/Desktop/seq/gene_tab_files/pre_tpm_synteny_list', 
#            col.names = F, row.names = F)

#######making graph of stage-specific filtering counts
filterstep=data.frame(1)
filterstep$tpm_1=nrow(Passed_filter2)
filterstep$transposon=nrow(Passed_filter2)-nrow(FILTER_2_transposons)
filterstep$outgroups=nrow(Passed_filter2)-nrow(FILTER_2_transposons)-nrow(FILTER_3_distant_outgroups)
filterstep$yakuba=nrow(Passed_filter2)-nrow(FILTER_2_transposons)-nrow(FILTER_3_distant_outgroups)-nrow(FILTER_4_yakuba)
filterstep$simulans=nrow(Passed_filter2)-nrow(FILTER_2_transposons)-nrow(FILTER_3_distant_outgroups)-nrow(FILTER_4_yakuba)-nrow(FILTER_5_sim)

filterstep=as.data.frame(t(filterstep))
filterstep$name=rownames(filterstep)
filterstep=filterstep[-1,]

ggplot(filterstep, aes(x=name[order(name)], y=V1, group = 1))+
  geom_bar(stat="identity", color='black', fill='grey50')+
  theme_classic()+
  scale_y_continuous(limits = c(0,650), expand = c(0, 0))+
  scale_x_discrete(labels=c("-under 1 tpm","-transposons","-distant outgroups","-D. yakuba","-D. simulans (final)"),)+
  xlab("filtering step")+
  ylab("# transcripts")+
  geom_text(aes(label=V1), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(axis.text = element_text(size=11))

##################################genomes##################################
#repeating first part of script in genomes folder instead
#kind of the inverse of the first step - genes that do align are not filtered 
#goal is to figure out how many of dn genes have alignable orthologous sequence in outgroups

fnames_g <- list.files("filtering/genome_blast/", full.names = T)
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
cantffind=unfiltered_gtf[unfiltered_gtf$GID%in%screen_df_g[which(screen_df_g[screen_df_g$nomatch==0,'GID']%in%initial_list==T),'GID'],]

#total counts
length(which(screen_df_g[screen_df_g$yakuba.genome>0,'GID']%in%initial_list))
length(which(screen_df_g[screen_df_g$simulans.genome_1>0,'GID']%in%initial_list))

####All transcripts are present in at least one genome
initial_list%in%cantffind$GID
initial_list[which(initial_list%in%cantffind$GID==T)] #1 gene not found in other genomes
initial_list=initial_list[-which(initial_list%in%cantffind$GID==T)]
####
unfiltered_gtf_final=unfiltered_gtf[unfiltered_gtf$GID%in%initial_list,]
mit_gene=which(unfiltered_gtf_final$V1=='mitochondrion_genome')
mit_gene #no mitochondria gene

initial_list=initial_list[!(initial_list%in%unfiltered_gtf_final[mit_gene,'GID'])]
dn_postmatch_pretpm=as.data.frame(initial_list)

#####add Yige's synteny analysis for filtering
setwd("~/Desktop/github1/")

UA_synteny=read.delim(file='synteny/MSD_synteny_results.txt', header = T)
DN_synteny=read.delim(file='synteny/MOD_synteny_results.txt',
                      stringsAsFactors = F, header = T)

dn_postmatch_pretpm1=dn_postmatch_pretpm[dn_postmatch_pretpm$initial_list%in%DN_synteny[which(DN_synteny$status=='confirmed'),'Query_gene_ID'],]
#length(dn_postmatch_pretpm$initial_list)
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

allgtf=read.table(file="dmel-all-r6.41.gtf", sep='\t')

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

#combine
colnames(transcript_idx)=c('TXNAME', 'GENEID')
transcript_idx=rbind(transcript_idx, transcript_idx1)
transcript_idx=unique(transcript_idx)


##############################################
#going back to the blast filtering dataframe to re check previous filtering stuff
non_mel_col_names=c(yak_col_names,sim_col_names)
outgroup_unannotated_match=(unlist(lapply(1:nrow(Passed_filter4), function(x) sum(Passed_filter4[x,non_mel_col_names]))))

MSU_match=Passed_filter4[which(outgroup_unannotated_match>0),]
MU_match=Passed_filter4[which(outgroup_unannotated_match==0),]

MSU_list1=rownames(Passed_filter4)[which(outgroup_unannotated_match>0)]


##################
#apply TPM screening here!!!
lengthscaledcounts <- read.delim(file = 'tables/unannotated_outgroup_expression.tsv')

lengthscaledcounts$max=apply(lengthscaledcounts[,1:4], 1, max)

lengthscaledcounts_MU=lengthscaledcounts[rownames(lengthscaledcounts)%in%initial_list,]
final_list1=rownames(lengthscaledcounts_MU)[which(lengthscaledcounts_MU$max<0.2)]


#####unannotated list required starting here
unann_exp=lengthscaledcounts[which(rownames(lengthscaledcounts)%in%MSU_list1),]
unann_exp$class='UA'
dn_exp=lengthscaledcounts[which(rownames(lengthscaledcounts)%in%pre_tpm_synteny_list$initial_list),]
dn_exp$class='DN'
exp_df=rbind(unann_exp, dn_exp)
colnames(exp_df)=gsub('.tsv', "", colnames(exp_df))
filter1=exp_df[exp_df$max>0.2,]

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

#Figure
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

unann_exp$n_zero=apply(unann_exp[,c(1,2,4)], 1, function(x) length(which(x<0.01)))
unann_exp$n_one=apply(unann_exp[,c(1,2,4)], 1, function(x) length(which(x>2)))
unann_exp2=unann_exp[unann_exp$n_zero>0.1,]

dn_tpm_filter_list=rownames(filter1[filter1$class=='DN',])
dn_tpm_filter_list

#at 0.2 TPM, 15 genes are filtered
length(dn_tpm_filter_list)
######################
#################################
######################
dn_postmatch_pretpm1[!dn_postmatch_pretpm1%in%dn_tpm_filter_list]
ambig_exp <- dn_postmatch_pretpm1[dn_postmatch_pretpm1%in%dn_tpm_filter_list]

#MSD list
finalMSD<-rbind(FILTER_4_yakuba, FILTER_5_sim,Passed_filter5)

finalMSD <- finalMSD[!rownames(finalMSD)%in%remove_for_UA_sim,]
finalMSD <- finalMSD[!rownames(finalMSD)%in%remove_for_UA_yak,]

finalMSD<-data.frame(MSD_genes = rownames(finalMSD)[(rownames(finalMSD)%in%rownames(lengthscaledcounts)[lengthscaledcounts$max>=0.5])])

#
write.table(finalMSD, "tables/MSD_list.txt", quote=F, col.names = F, row.names = F)
write.table(final_list1, file = "tables/MOD_list.txt", quote=F, col.names = F, row.names = F)
