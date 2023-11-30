library(rtracklayer)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(tidyr)

bidi_len=1000

RNA_counts=read.delim(file='/Users/loganblair/Desktop/seq/count_tables/lengthscaledcounts',
                      sep='\t', stringsAsFactors = F, header=T)
RNA_counts=RNA_counts[,c(1:29)]
FBGNid=read.delim(file='/Users/loganblair/Desktop/seq/count_tables/kallisto_FBGNid', sep='\t')
AG_specific_list=read.table(file='/Users/loganblair/Desktop/seq/gene_tab_files/AG_biased_gene_list', sep='\t') #AG

unnan_list_nogenome=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/unnan_list_nogenome', sep='\t',header=F)
final_screen_list=read.table(file='/Users/loganblair/Desktop/seq/final_list',
                             sep='\t',stringsAsFactors = F)

unannotated_list=c(unnan_list_nogenome[,1], final_screen_list[,1])
#####splitting annotated genes by which are expressed in AGs
min_list=c(apply(RNA_counts, 1, which.min))
RNA_counts$min_exp=NA
for(i in 1:nrow(RNA_counts)){
  RNA_counts$min_exp[i]=RNA_counts[i,min_list[i]]
}
RNA_counts$mean=apply(RNA_counts[,1:29], 1, mean)
RNA_counts$notep=apply(RNA_counts[,c(1:29)], 1, function(x) mean(x)>1)
RNA_counts$FBGN_list=1:nrow(RNA_counts)%in%grep("^F", rownames(RNA_counts))
NE_list=rownames(RNA_counts)[which(RNA_counts$FBGN_list==T&RNA_counts$mean<1)]
BE_list=rownames(RNA_counts)[which(RNA_counts$FBGN_list==T&RNA_counts$mean>1)]
BE_list=BE_list[!BE_list%in%AG_specific_list$V1]


RNA_counts$mean=apply(RNA_counts[,c(1:29)], 1, function(x) as.numeric(mean(x)))

non_exp_genelist=row.names(RNA_counts[which(RNA_counts$mean<1),])
non_exp_genelist=non_exp_genelist[grep("FB", non_exp_genelist)]

exp_genelist_no_AG=row.names(RNA_counts[which(RNA_counts$mean>1),])
exp_genelist_no_AG=exp_genelist_no_AG[!exp_genelist_no_AG%in%AG_specific_list$V1]
exp_genelist_no_AG=exp_genelist_no_AG[grep("FBgn", exp_genelist_no_AG)]

#location of all genes
gr_allgenes = import.gff('/Users/loganblair/Desktop/seq/gtf_files/dm6_dn2022gtf')

#reduce all genes to one range
gr_allgenes=unlist(reduce(split(gr_allgenes, gr_allgenes$gene_id)))

#AG genes
genes_AG=gr_allgenes[names(gr_allgenes)%in%AG_specific_list[,1],]
genes_AG=as.data.frame(genes_AG)
genes_AG$GID=AG_specific_list[,1]

#dn genes
genes_DN=gr_allgenes[names(gr_allgenes)%in%final_screen_list$GID,]

#########################
#AG specific - Julie's AG rho>0.9 set
genes_AGspecific=gr_allgenes[names(gr_allgenes)%in%AG_specific_list[,1],]

#
bidi_df=as.data.frame(1:length(gr_allgenes))
bidi_df$chr=as.character(seqnames(gr_allgenes))
bidi_df$geneid=as.character(names(gr_allgenes))
bidi_df$TSS=NA
bidi_df$strand=as.character(strand(gr_allgenes))
bidi_df=bidi_df[,-1]
bidi_df$start=as.numeric(start(gr_allgenes))
bidi_df$end=as.numeric(end(gr_allgenes))

#choose start +/- genes
for(i in 1:nrow(bidi_df)){
  if(bidi_df[i,'strand']=='+'){
    bidi_df[i,'TSS']=bidi_df[i,'start']
  }
  if(bidi_df[i,'strand']=='-'){
    bidi_df[i,'TSS']=bidi_df[i,'end']
  }
}

bidi_df$bidirectional=NA
bidi_df$bidi_plus=NA
bidi_df$bidi_minus=NA

#filling out genetype column

bidi_df$genetype=NA

#omitting mostly pseudogenes that are not in the RNA dataframe
bidi_df=bidi_df[-which(bidi_df$geneid%in% rownames(RNA_counts)==F),]

#
bidi_df$genetype[bidi_df$geneid%in% BE_list] ='BE'
bidi_df$genetype[bidi_df$geneid%in% NE_list]='NE'
bidi_df$genetype[bidi_df$geneid%in%AG_specific_list[,1]]='AG'
bidi_df$genetype[bidi_df$geneid%in% unnan_list_nogenome$V1]='UA'
bidi_df$genetype[bidi_df$geneid%in% final_screen_list[,1]]='DN'


for(i in 1:nrow(bidi_df)){
  temp_df=bidi_df[bidi_df$chr==bidi_df$chr[i],] #take genes only on same chromosome
  
  temp_df=temp_df[-which(temp_df$strand==bidi_df$strand[i]),] #look at opposite stranded only
  
  if(length(temp_df[which(temp_df$TSS>(bidi_df$TSS[i]-bidi_len)&
                   temp_df$TSS<bidi_df$TSS[i]),'geneid'])>0){
    bidi_df$bidi_plus[i]=paste(temp_df[which(temp_df$TSS>(bidi_df$TSS[i]-bidi_len)&
                                      temp_df$TSS<bidi_df$TSS[i]),'geneid'],collapse = ";")
  }
  
  if(length(temp_df[which(temp_df$TSS<(bidi_df$TSS[i]+bidi_len)&
                          temp_df$TSS>bidi_df$TSS[i]),'geneid'])>0){
    bidi_df$bidi_minus[i]=paste(temp_df[which(temp_df$TSS<(bidi_df$TSS[i]+bidi_len)&
                                               temp_df$TSS>bidi_df$TSS[i]),'geneid'],collapse = ";")
  }
}

#select start based on strandedness
for(i in 1:nrow(bidi_df)){
  if(bidi_df$strand[i]=='+'){
    bidi_df$bidirectional[i]=bidi_df$bidi_plus[i]
  }
  
  if(bidi_df$strand[i]=='-'){
    bidi_df$bidirectional[i]=bidi_df$bidi_minus[i]
  }
}

bidi_df=separate_rows(bidi_df, sep=';', bidirectional)



#at this point have one row for each gene, plus duplicate rows for genes that are paired two+ times
bidi_df$BE="0"
bidi_df$UA="0"
bidi_df$AG="0"
bidi_df$DN="0"
bidi_df$NE="0"

for(i in 1:nrow(bidi_df)){
  if(is.na(bidi_df$bidirectional[i])==F){
    tmp_genelist=str_split_fixed(string = bidi_df$bidirectional[i], pattern = ';', 4)
    tmp_genelist
    #since some genes have multiple others being transcribed in opposite direction within 1000bp cutoff
    #I am looking for any match for each class of genes
    
    bidi_df$BE[i]=length(which(tmp_genelist %in% BE_list ==T))
    bidi_df$NE[i]=length(which(tmp_genelist %in% NE_list ==T))
    bidi_df$AG[i]=length(which(tmp_genelist%in%AG_specific_list[,1]==T))

    bidi_df$UA[i]=length(which(tmp_genelist %in% unnan_list_nogenome[,1] ==T))
    bidi_df$DN[i]=length(which(tmp_genelist %in% final_screen_list$x ==T))
    }
}


bidi_df$AG=as.numeric(bidi_df$AG)
bidi_df$UA=as.numeric(bidi_df$UA)
bidi_df$BE=as.numeric(bidi_df$BE)
bidi_df$DN=as.numeric(bidi_df$DN)
bidi_df$NE=as.numeric(bidi_df$NE)

bidi_df$bidi_genetype=NA
counter_cols=which(colnames(bidi_df)=='BE'):which(colnames(bidi_df)=='NE')

#converting matched numeric to gene categories
for(i in 1:nrow(bidi_df)){
  if(sum(bidi_df[i,counter_cols]==0)){bidi_df$bidi_genetype[i]="None"}
  if(bidi_df[i,'BE']>0){bidi_df$bidi_genetype[i]='BE'}
  if(bidi_df[i,'NE']>0){bidi_df$bidi_genetype[i]='NE'}
  if(bidi_df[i,'AG']>0){bidi_df$bidi_genetype[i]='AG'}
  if(bidi_df[i,'UA']>0){bidi_df$bidi_genetype[i]='UA'}
  if(bidi_df[i,'DN']>0){bidi_df$bidi_genetype[i]='DN'}
}

bidi_counts_A=(subset(bidi_df,select = c('genetype', 'bidi_genetype')))
bidi_counts_A$counter=1:nrow(bidi_counts_A)

######double check no NAs in gene type of first gene
bidi_counts=as.data.frame(table(subset(bidi_df,select = c('genetype', 'bidi_genetype'))))

names(bidi_counts) <- c("genetype", "bidi_genetype", "Freq")

bidi_counts$prop=NA
bidi_counts[bidi_counts$genetype=='UA','prop']=bidi_counts[bidi_counts$genetype=='UA','Freq']/
  sum(bidi_counts[bidi_counts$genetype=='UA','Freq'])

bidi_counts[bidi_counts$genetype=='UA','Freq']/
  sum(bidi_counts[bidi_counts$genetype=='UA','Freq'])


bidi_counts[bidi_counts$genetype=='NE','prop']=bidi_counts[bidi_counts$genetype=='NE','Freq']/
  sum(bidi_counts[bidi_counts$genetype=='NE','Freq'])
bidi_counts[bidi_counts$genetype=='AG','prop']=bidi_counts[bidi_counts$genetype=='AG','Freq']/
  sum(bidi_counts[bidi_counts$genetype=='AG','Freq'])
bidi_counts[bidi_counts$genetype=='BE','prop']=bidi_counts[bidi_counts$genetype=='BE','Freq']/
  sum(bidi_counts[bidi_counts$genetype=='BE','Freq'])
bidi_counts[bidi_counts$genetype=='DN','prop']=bidi_counts[bidi_counts$genetype=='DN','Freq']/
  sum(bidi_counts[bidi_counts$genetype=='DN','Freq'])

bidi_counts$Freq=as.numeric(bidi_counts$Freq)
bidi_counts_freq_heatmap=bidi_counts

####################################################################################
####################################################################################
#####################Plot1##########################################################
library(ggtext)

graph1000=ggplot(bidi_counts, aes(x=genetype, y=prop, fill=bidi_genetype, label=Freq))+
  geom_bar(stat='identity', color='white', alpha=0.87, size=.1)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0), limits = c(0,1.04), breaks=pretty_breaks(n=10))+
  #scale_fill_manual(values=c("goldenrod", "#C75127","#5DB1DD","#7A65A5","grey35", "grey75"), 
  #                            labels = c("DN", "UA","AG","BE","NE", "none/unidirectional"))+
  #scale_x_discrete(labels= c("DN","UA","AG","BE", "NE"))+
  ylab("proportion")+
  xlab("1st gene in pair")+
  geom_text(size = 3, position = position_stack(vjust = 0.45), color='black')+
  labs(fill='2nd gene in pair')
  #theme(axis.text.x = element_text(colour = c("grey35","#7A65A5","#5DB1DD","#C75127", "goldenrod")))

graph1000

#get names in right order via stupid way
bidi_counts$genetype=gsub("DN","1DN",bidi_counts$genetype)
bidi_counts$bidi_genetype=gsub("DN","1DN",bidi_counts$bidi_genetype)
bidi_counts$genetype=gsub("UA","2UA",bidi_counts$genetype)
bidi_counts$bidi_genetype=gsub("UA","2UA",bidi_counts$bidi_genetype)
bidi_counts$genetype=gsub("AG","3AG",bidi_counts$genetype)
bidi_counts$bidi_genetype=gsub("AG","3AG",bidi_counts$bidi_genetype)
bidi_counts$genetype=gsub("BE","4BE",bidi_counts$genetype)
bidi_counts$bidi_genetype=gsub("BE","4BE",bidi_counts$bidi_genetype)
bidi_counts$genetype=gsub("NE","5NE",bidi_counts$genetype)
bidi_counts$bidi_genetype=gsub("NE","5NE",bidi_counts$bidi_genetype)

bidi_counts[bidi_counts$genetype=='2UA'&bidi_counts$bidi_genetype=='1DN',3:4]=NA
bidi_counts[bidi_counts$genetype=='3AG'&bidi_counts$bidi_genetype=='1DN',3:4]=NA
bidi_counts[bidi_counts$genetype=='3AG'&bidi_counts$bidi_genetype=='2UA',3:4]=NA
bidi_counts[bidi_counts$genetype=='4BE'&bidi_counts$bidi_genetype=='1DN',3:4]=NA
bidi_counts[bidi_counts$genetype=='4BE'&bidi_counts$bidi_genetype=='2UA',3:4]=NA
bidi_counts[bidi_counts$genetype=='4BE'&bidi_counts$bidi_genetype=='3AG',3:4]=NA
bidi_counts[bidi_counts$genetype=='5NE'&bidi_counts$bidi_genetype=='1DN',3:4]=NA
bidi_counts[bidi_counts$genetype=='5NE'&bidi_counts$bidi_genetype=='2UA',3:4]=NA
bidi_counts[bidi_counts$genetype=='5NE'&bidi_counts$bidi_genetype=='3AG',3:4]=NA
bidi_counts[bidi_counts$genetype=='5NE'&bidi_counts$bidi_genetype=='4BE',3:4]=NA

ggplot(bidi_counts, aes(y=genetype, fill=prop, x=bidi_genetype, label=Freq))+
  geom_tile(stat='identity', )+
  theme_classic()+
  scale_y_discrete(
    labels= c("MU","MSU", "AG", "BE", "NE"),expand=c(0,0))+
  scale_x_discrete(labels= c("MU", "MSU","AG","BE","NE","None"),
                   expand=c(0,0), position = "bottom")+
  ylab(" ")+
  xlab(" ")+
  #geom_text(size = 3, position = position_stack(vjust = 0.45), color='black')+
  labs(fill='prop/class')+
  theme(text=element_text(size=20),
        legend.text = element_text(size=20),
        axis.text.x = element_text(angle = 0, vjust =0, hjust=0.5, size=20),
        axis.text.y = element_text(size=20))+
  scale_fill_gradient2(low="blue4",
                       mid = "white",
                       high = "red", midpoint = 0.15, )+
  geom_text(size=6,aes(label = Freq 
                #,size=6.5
                )  )
########
total_gene_num=length(c(BE_list,NE_list,AG_specific_list[,1],unnan_list_nogenome[,1], final_screen_list[,1]))
total_gene_num
bidi_counts$bidi_freq=NA

bidi_counts[bidi_counts$bidi_genetype=='1DN','bidi_freq']=length(final_screen_list[,1])/total_gene_num
bidi_counts[bidi_counts$bidi_genetype=='2UA','bidi_freq']=length(unnan_list_nogenome[,1])/total_gene_num
bidi_counts[bidi_counts$bidi_genetype=='3AG','bidi_freq']=length(AG_specific_list[,1])/total_gene_num
bidi_counts[bidi_counts$bidi_genetype=='4BE','bidi_freq']=length(BE_list)/total_gene_num
bidi_counts[bidi_counts$bidi_genetype=='5NE','bidi_freq']=length(NE_list)/total_gene_num
bidi_counts$relative_freq=bidi_counts$prop/bidi_counts$bidi_freq

bidi_counts1=bidi_counts[!bidi_counts$bidi_genetype=='None',]
ggplot(bidi_counts1, aes(y=genetype, fill=relative_freq, x=bidi_genetype, label=relative_freq))+
  geom_tile(stat='identity', )+
  theme_classic()+
  scale_y_discrete(
    labels= c("MU","MSU", "AG", "BE", "NE"),
    expand=c(0,0))+
  scale_x_discrete(
    labels= c("MU", "MSU","AG","BE","NE","None"),
    expand=c(0,0), position = "bottom")+
  ylab(" ")+
  xlab(" ")+
  #geom_text(size = 3, position = position_stack(vjust = 0.45), color='black')+
  labs(fill='norm. prop/class')+
  theme(text=element_text(size=20),
        legend.text = element_text(size=20),
        axis.text.x = element_text(angle = 0, vjust =0, hjust=0.5, size=20),
        axis.text.y = element_text(size=20))+
  scale_fill_gradient2(low="blue4",
                       mid = "white",
                       high = "red", midpoint = 2, )+
  geom_text(size=6,aes(label = Freq
                       #,size=6.5
  )  )

ggplot(bidi_counts1, aes(y=genetype, fill=relative_freq, x=bidi_genetype, label=relative_freq))+
  geom_tile(stat='identity', )+
  theme_classic()+
  scale_y_discrete(
    labels= c("MU","MSU", "AG", "BE", "NE"),
    expand=c(0,0))+
  scale_x_discrete(
    labels= c("MU", "MSU","AG","BE","NE","None"),
    expand=c(0,0), position = "bottom")+
  ylab(" ")+
  xlab(" ")+
  #geom_text(size = 3, position = position_stack(vjust = 0.45), color='black')+
  labs(fill='normalized freq')+
  theme(text=element_text(size=20),
        legend.text = element_text(size=20),
        axis.text.x = element_text(angle = 0, vjust =0, hjust=0.5, size=20),
        axis.text.y = element_text(size=20))+
  scale_fill_gradient2(low="blue4",
                       mid = "white",
                       high = "red", midpoint = 2, )+
  geom_text(size=6,aes(label = round(relative_freq, digits = 2)
                       #,size=6.5
  )  )


#############################
#Calculating correlated expression of bidirectional pairs
library(tidyr)

bidi_df1=bidi_df[bidi_df$geneid%in%rownames(RNA_counts),]
bidi_df1=bidi_df1[bidi_df1$bidirectional%in%rownames(RNA_counts),]

RNA_counts$min_exp=NULL

rank_df=apply(RNA_counts[,c(1:29)], 1, rank)
apply(rank_df, 1, sum)

#using just tpm values inflates spearman rank, because some genotypes have high expression of few genes
#trying the edgeR cpm corrected values

RNA_cpm=read.delim(file='/Users/loganblair/Desktop/seq/count_tables/6.41_kallisto_edgeR_cpm', sep='\t')

rank_df1=apply(RNA_cpm, 1, rank)
apply(rank_df1, 1, sum) #looks better

bidi_df1$bidi_cor=NA
for(i in 1:nrow(bidi_df1)){
  gene1=rank(RNA_cpm[which(rownames(RNA_cpm)==bidi_df1$geneid[i]),])
  gene2=rank(RNA_cpm[which(rownames(RNA_cpm)==bidi_df1$bidirectional[i]),])
  bidi_df1$bidi_cor[i]=cor(gene1, gene2, method=c("spearman"))
}


##############################
bi_didf3=bidi_df1[bidi_df1$genetype%in%c("BE","AG","DN","UA"),]
bi_didf3=bi_didf3[bi_didf3$bidi_genetype%in%c("BE","DN","AG", "UA"),]

wil_bidi_df=bi_didf3[,colnames(bi_didf3)%in%c("bidi_cor", "genetype", "bidi_genetype")]

#reordering
wil_bidi_df$genetype=gsub("BE","1BE", wil_bidi_df$genetype)
wil_bidi_df$bidi_genetype=gsub("BE","1BE", wil_bidi_df$bidi_genetype)

#clumping MU and MSU together
wil_bidi_df$genetype=gsub("DN","UA", wil_bidi_df$genetype)
wil_bidi_df$bidi_genetype=gsub("DN","UA", wil_bidi_df$bidi_genetype)

remove_list1=c("1BE UA", "AG UA")
wil_bidi_df$paste=paste(wil_bidi_df$genetype, wil_bidi_df$bidi_genetype, sep=' ')
wil_bidi_df=wil_bidi_df[!wil_bidi_df$paste%in%remove_list1,]

ggplot(wil_bidi_df, aes(y=bidi_cor, x=genetype, fill=bidi_genetype))+
  #geom_boxplot(alpha=.75, width=.8)+
  geom_boxplot(position = position_dodge2(preserve = "single", padding = 0))+
  scale_fill_manual(values=c("#7A65A5","#5DB1DD", "darkorange3"),
                    labels=c("BE","AG","MU+MSU"))+
  theme_bw()+
  ylab("spearman's correlation coefficient")+
  xlab("1st gene class")+
  scale_x_discrete(labels=c("BE", "AG", "MU+MSU"))+
  scale_y_continuous(limits=c(-0.74,1.24))+
  labs(fill="2nd gene class")+
  theme(text=element_text(size=21))+
  geom_segment(aes(x = 3, xend = 3, y = 1.07, yend=1.1))+
  geom_segment(aes(x = 3.25, xend = 3.25, y = 1.07, yend=1.1))+
  geom_segment(aes(x = 3, xend = 3.25, y = 1.1, yend=1.1))+
  annotate("text", x = 3.125, y=1.15, label="***", size=8)



##############################
#negative binomial framing (obvs vs exp)
bidi_counts_B=as.data.frame(table(subset(bidi_df,select = c('genetype', 'bidi_genetype'))))
bidi_counts_B$genetype=gsub("UA", "MSU", bidi_counts_B$genetype)
bidi_counts_B$genetype=gsub("DN", "MU", bidi_counts_B$genetype)
bidi_counts_B$bidi_genetype=gsub("UA", "MSU", bidi_counts_B$bidi_genetype)
bidi_counts_B$bidi_genetype=gsub("DN", "MU", bidi_counts_B$bidi_genetype)

bidi_counts_D=bidi_counts_B[!bidi_counts_B$bidi_genetype=="None",]
bidi_counts_D$paste=paste(bidi_counts_D$genetype, bidi_counts_D$bidi_genetype)
bidi_counts_D1=bidi_counts_D

remove_list=c("MSU MU","AG MU","BE MU","NE MU","AG MSU","BE MSU",
              "NE MSU","BE AG","NE AG","NE BE")
double_list=c("MU MU","MSU MSU","AG AG","BE BE","NE NE")
bidi_counts_D=bidi_counts_D[!bidi_counts_D$paste%in%remove_list,]

#same class pairs get counted double, so that needs to be halved when assessing total number of
#bidirectional "trials"
total_bidi=sum(bidi_counts_D$Freq)-(sum(bidi_counts_D[bidi_counts_D$paste%in%double_list,'Freq'])/2)
#(8+3+4+6+11)+(3+14+13+7+3)+(4+13+118+26+19)+(6+7+26+1640+1151)+(11+3+19+1151+1774)/2

MU_bidi_frac=sum(bidi_counts_D[(grep("MU", bidi_counts_D$paste)),'Freq'])/(total_bidi*2)
MSU_bidi_frac=sum(bidi_counts_D[(grep("MSU", bidi_counts_D$paste)),'Freq'])/(total_bidi*2)
AG_bidi_frac=sum(bidi_counts_D[(grep("AG", bidi_counts_D$paste)),'Freq'])/(total_bidi*2)
BE_bidi_frac=sum(bidi_counts_D[(grep("BE", bidi_counts_D$paste)),'Freq'])/(total_bidi*2)
NE_bidi_frac=sum(bidi_counts_D[(grep("NE", bidi_counts_D$paste)),'Freq'])/(total_bidi*2)
MU_bidi_frac+MSU_bidi_frac+AG_bidi_frac+BE_bidi_frac+NE_bidi_frac

#making table of what proportion each gene class consists of all bidirectional genes
frac_table=as.data.frame(c("MU","MSU","AG","BE","NE"))
colnames(frac_table)[1]="class"
frac_table$fraction=c(MU_bidi_frac,MSU_bidi_frac,AG_bidi_frac,BE_bidi_frac,NE_bidi_frac)

for(i in 1:nrow(bidi_counts_D1)){
  bidi_counts_D1$expected[i]=
    total_bidi*
    frac_table[which(frac_table$class==bidi_counts_D1$genetype[i]), 'fraction']*
    frac_table[which(frac_table$class==bidi_counts_D1$bidi_genetype[i]), 'fraction']*2
}

#observed and expected should match
sum(bidi_counts_D1$expected)
sum(bidi_counts_D1$Freq)

#when pivoting to trials, not total genes, same-class pairs need to be halved
bidi_counts_D1=bidi_counts_D1[!bidi_counts_D1$paste%in%remove_list,]
bidi_counts_D1[bidi_counts_D1$paste%in%double_list,'expected']=
  bidi_counts_D1[bidi_counts_D1$paste%in%double_list,'expected']/2
bidi_counts_D1[bidi_counts_D1$paste%in%double_list,'Freq']=
  bidi_counts_D1[bidi_counts_D1$paste%in%double_list,'Freq']/2

#match
sum(bidi_counts_D1$expected)
sum(bidi_counts_D1$Freq)

#enrichment and binomial testing
bidi_counts_D1$enrichment=bidi_counts_D1$Freq/bidi_counts_D1$expected
bidi_counts_D1$binomial=NA

for(i in 1:nrow(bidi_counts_D1)){
  bidi_counts_D1$binomial[i]=binom_test(x = bidi_counts_D1$Freq[i], n = total_bidi,
                                        p = bidi_counts_D1$expected[i]/(total_bidi))$p
}

bidi_counts_D1$adj=p.adjust(bidi_counts_D1$binomial, method='holm')
bidi_counts_D1$sig=" "
bidi_counts_D1[bidi_counts_D1$adj<0.05,'sig']="*"
bidi_counts_D1[bidi_counts_D1$adj<0.01,'sig']="**"
bidi_counts_D1[bidi_counts_D1$adj<0.001,'sig']="***"

#add remove list back for graph. Creates "triangluar" look of contingency table
row_length=nrow(bidi_counts_D1)
class1_list=str_split_fixed(remove_list, " ", 2)[,1]
class2_list=str_split_fixed(remove_list, " ", 2)[,2]

for(i in 1:length(remove_list)){
  bidi_counts_D1[row_length+i,]=c(class1_list[i],class2_list[i],rep(NA, times=(ncol(bidi_counts_D1)-2)))
}

bidi_counts_D1$genetype=gsub("MU","1MU",bidi_counts_D1$genetype)
bidi_counts_D1$genetype=gsub("MSU","2MSU",bidi_counts_D1$genetype)
bidi_counts_D1$genetype=gsub("AG","3AG",bidi_counts_D1$genetype)
bidi_counts_D1$genetype=gsub("BE","4BE",bidi_counts_D1$genetype)
bidi_counts_D1$bidi_genetype=gsub("MU","1MU",bidi_counts_D1$bidi_genetype)
bidi_counts_D1$bidi_genetype=gsub("MSU","2MSU",bidi_counts_D1$bidi_genetype)
bidi_counts_D1$bidi_genetype=gsub("AG","3AG",bidi_counts_D1$bidi_genetype)
bidi_counts_D1$bidi_genetype=gsub("BE","4BE",bidi_counts_D1$bidi_genetype)

bidi_counts_D1$adj=as.numeric(bidi_counts_D1$adj)

bidi_counts_D1$enrichment=as.numeric(bidi_counts_D1$enrichment)
ggplot(bidi_counts_D1, aes(y=genetype, fill=(enrichment), x=bidi_genetype, label=enrichment))+
  geom_tile(stat='identity', color='black')+
  theme_classic()+
  scale_y_discrete(
    labels= c("MU","MSU", "AG", "BE", "NE"),
    expand=c(0,0))+
  scale_x_discrete(
    labels= c("MU", "MSU","AG","BE","NE","None"),
    expand=c(0,0), position = "bottom")+
  ylab(" ")+
  xlab(" ")+
  #geom_text(size = 3, position = position_stack(vjust = 0.45), color='black')+
  labs(fill='obs/exp')+
  theme(text=element_text(size=20),
        legend.text = element_text(size=20),
        axis.text.x = element_text(angle = 0, vjust =0, hjust=0.5, size=20),
        axis.text.y = element_text(size=20),
        legend.key.height = unit(1,"cm"))+
  scale_fill_gradient2(low="blue4",
                       mid ="white",
                       high = "red", midpoint = 0, breaks=c(0.3,1,3,10,30), trans="log")+
 geom_text(size=10,aes(label = sig))
