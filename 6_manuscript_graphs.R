library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggpmisc)
library(scales)

options(scipen=999999)
allgtf=read.delim(file='allgtf', sep='\t', header=F)

#load RNA count table from kalisto
RNA_counts=read.delim(file='kallisto_expression.tsv',
                      sep='\t', stringsAsFactors = F, header=T)

#####tau vs exp######
#load table from Leader et al 2018, showing expression of genes across different tissues in D mel
tau_df=read.delim(file="FlyAtlas_data_table", sep='\t', header=T)
row.names(tau_df)=tau_df[,1]
tau_df=tau_df[,-1]

#actual tau calculation from Yanni et al 2005. Takes a while.
for(i in 1:nrow(tau_df)){
  tau_df$tau[i]=sum(unlist(lapply(tau_df[i,], function(x) (1-(x/max(tau_df[i,]))))))/(ncol(tau_df)-1)}
for(i in 1:nrow(tau_df)){tau_df$max_tissue[i]=names(which.max(tau_df[i,c(1:(ncol(tau_df)-1))]))}

#remove a few weird genes that are not in the count table
tau_df1=tau_df[rownames(tau_df)%in%rownames(RNA_counts),]

RNA_counts_tau=RNA_counts[rownames(RNA_counts)%in%rownames(tau_df1),]
RNA_counts_tau$mean=apply(RNA_counts_tau, 1, mean)

tau_df1$RNA_mean=RNA_counts_tau[match(rownames(tau_df1), rownames(RNA_counts_tau)),'mean']
tau_df1$AG_enriched=tau_df1$max_tissue=='Male_Accessory_Glands'

#define parameters for tissue specific expression
tau_df1$AG_enriched[which(tau_df1$AG_enriched=="TRUE"&tau_df1$tau>0.9
                          &tau_df1$RNA_mean>1
                          )]="AG (specific)"
tau_df1$AG_enriched=gsub("TRUE", "AG (nonspecific/low expression)", tau_df1$AG_enriched)
tau_df1$AG_enriched=gsub("FALSE", "other tissue", tau_df1$AG_enriched)
tau_df1=tau_df1[order(tau_df1$AG_enriched, decreasing = T),]

#enrichment graph (figure S4)
options(scipen=1)
graph_tau_enrichment <- ggplot(tau_df1, aes(y=RNA_mean+1, x=tau, color=AG_enriched))+
  geom_point(size=1.5,alpha=1)+
  scale_y_log10(breaks=c(1,10,100,1000,10000,100000,1000000), limits=c(0.5,4500000))+
  theme_classic()+
  scale_x_continuous(limits=c(0.42,1.04), breaks=seq(0.4, 1, 0.1), expand = c(0,0))+
  scale_color_manual(values=c('darkblue',"#56B4E9",'grey78'), 
                     labels = c("AG", "AG (Exp+ Sp+)", "other tissue"))+
  geom_hline(yintercept = 2, linetype='dashed')+
  geom_vline(xintercept = 0.9, linetype='dashed')+
  ylab('AG exp +1 (TPM)')+
  xlab("tissue specificity (tau)")+
  annotate("text", x=0.5, y=4100000, label= "Exp+ Sp-", size=8)+
  annotate("text", x=.965, y=0.5, label= "Exp- Sp+", size=8)+
  annotate("text", x=.965, y=4100000, label= "Exp+ Sp+", size=8)+
  annotate("text", x=0.5, y=0.5, label= "Exp- Sp-", size=8)+
  labs(color="highest exp")+
  theme(legend.position = 'right',
        text=element_text(size=25))+
  guides(colour = guide_legend(override.aes = list(size=2.27)))

AG_specific_list=row.names(tau_df1[tau_df1$AG_enriched=='AG (specific)',])

#############################################################################################






#de novo and AG specific genes
RNA_counts=read.delim(file='kallisto_expression.tsv',
                      sep='\t', stringsAsFactors = F, header=T)
RNA_counts=RNA_counts[,c(1:29)]#no iso1 or ed10

final_screen_list=read.delim(file='tables/MOD_list.txt',
                                        sep='\t',stringsAsFactors = F, header=T)

unan_out_list=read.table(file='tables/MSD_list.txt', sep='\t')

#
rownames(RNA_counts)=gsub(' ','', rownames(RNA_counts))

RNA_counts$mean=apply(RNA_counts[,1:29], 1, mean)
RNA_counts$var=apply(RNA_counts[,1:29], 1, var)
RNA_counts$sd=apply(RNA_counts[,1:29], 1, sd)
RNA_counts$genetype='1'
RNA_counts$genetype[rownames(RNA_counts)%in%AG_specific_list[,1]]='2'
RNA_counts$genetype[rownames(RNA_counts)%in%final_screen_list[,1]]='4'
RNA_counts$genetype[rownames(RNA_counts)%in%unan_out_list[,1]]='3'

#######graph 2: proportion of lines expressed#####

RNA_counts$num_greater5=apply(RNA_counts[,1:29], 1, function(x) length(which((x)>=5)))
RNA_counts$num_greater1=apply(RNA_counts[,1:29], 1, function(x) length(which((x)>=1)))

RNA_counts$max=apply(RNA_counts[,1:29], 1, function(x) max(x))
RNA_counts_anscestral=RNA_counts[grep('fbg', rownames(RNA_counts), ignore.case = T),]
RNA_counts_dn=RNA_counts[rownames(RNA_counts)%in%final_screen_list[,1],]
RNA_counts_dn=RNA_counts_dn[!RNA_counts_dn$num_greater1==0,]

graph_lines_exp=ggplot(RNA_counts_dn, aes(x=num_greater1))+
  geom_bar(stat='count', color='black', size=0.25, width=0.75, alpha=0.5)+
  scale_fill_manual(values=c("#56B4E9"))+
  scale_x_continuous(expand = c(0, 0), limits = c(0,31), breaks=seq(0,30,5))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  ylab('gene count')+
  xlab('#lines expressed >1 TPM')

#######graph1:variance vs expression########################################################################
#options(scipen = 10000000000000)
graph_v_x_exp=ggplot(RNA_counts %>%
         arrange(genetype), aes(x=mean+1, y=var+1, color=genetype))+
  geom_point(size=1.5, alpha=1)+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x)) +
  ylab('variance + 1')+
  scale_x_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x)) +
  xlab('mean expression (TPM) + 1')+
  scale_color_manual(values=c( "grey50","#5DB1DD","#C75127","goldenrod"), 
                     labels = c("AN (-AG)",  "AG", "MSD","MOD"))+
  #facet_zoom(xlim = c(.01, .5),ylim=c(.01,1), zoom.size = 0.25, split=F,horizontal=F, shrink=T)+
  theme_classic()+
  #geom_segment(aes(x = 15, y = 3, xend = 1000, yend = 3), color='black')+
  #annotate("rect", xmin=1, xmax=2, ymin=1, ymax=2, alpha=0.63, fill="grey40") +
  theme(text = element_text(size=25))
graph_v_x_exp
######graph 3: chromosomes########################################################################

allgtf=read.delim(file='all.gtf', sep='\t', header=F)
#add de novo names

allgtf_dn=allgtf[grep('assembly',allgtf$V2 ),]
allgtf_dn=allgtf_dn[allgtf_dn$V3=='transcript',]
allgtf_dn$gene_id=sapply(strsplit(gsub('gene_id ', '', allgtf_dn$V9), ';'),`[`, 6)
allgtf_dn$gene_id=gsub(" ", "", allgtf_dn$gene_id)
allgtf_dn$transcript_id=sapply(strsplit(gsub('transcript_id ', '', allgtf_dn$V9), ';'),`[`, 9)

allgtf_unnan=allgtf_dn[allgtf_dn$gene_id%in%unan_out_list[,1],]
allgtf_dn=allgtf_dn[allgtf_dn$gene_id%in%final_screen_list[,1],]

allgtf_gene=allgtf[allgtf$V3=='gene',]

allgtf_dn1=allgtf_dn[!duplicated(allgtf_dn$gene_id),]
allgtf_unann1=allgtf_unnan[!duplicated(allgtf_unnan$gene_id),]
allgtf_gene$gene_id=sapply(strsplit(gsub('gene_id ', '', allgtf_gene$V9), ';'),`[`, 1)
allgtf_gene$gene_id=gsub(" ", "", allgtf_gene$gene_id)

allgtf_ag=allgtf_gene[allgtf_gene$gene_id%in%AG_specific_list$V1,]
allgtf_unnan=allgtf_unnan[!duplicated(allgtf_unnan$gene),]

chr_df1=as.data.frame(unique(allgtf_dn1$V1))
chr_df2=as.data.frame(unique(allgtf_dn1$V1))
chr_df3=as.data.frame(unique(allgtf_dn1$V1))
chr_df4=as.data.frame(unique(allgtf_dn1$V1))

colnames(chr_df1)[1]='chromosome'
colnames(chr_df2)[1]='chromosome'
colnames(chr_df3)[1]='chromosome'
colnames(chr_df4)[1]='chromosome'

chr_df1$proportion=lapply(1:nrow(chr_df1), function(x) length(which(allgtf_gene$V1==chr_df1$chromosome[x])))
chr_df1$proportion1=as.numeric(chr_df1$proportion)/nrow(allgtf_gene)
chr_df2$proportion=lapply(1:nrow(chr_df2), function(x) length(which(allgtf_dn1$V1==chr_df2$chromosome[x])))
chr_df2$proportion1=as.numeric(chr_df2$proportion)/nrow(allgtf_dn1)
chr_df3$proportion=lapply(1:nrow(chr_df3), function(x) length(which(allgtf_ag$V1==chr_df3$chromosome[x])))
chr_df3$proportion1=as.numeric(chr_df3$proportion)/nrow(allgtf_ag)
chr_df4$proportion=lapply(1:nrow(chr_df4), function(x) length(which(allgtf_unnan$V1==chr_df4$chromosome[x])))
chr_df4$proportion1=as.numeric(chr_df4$proportion)/nrow(allgtf_unnan)

chr_df1$type='1ancestral'
chr_df2$type='4MOD'
chr_df3$type='2AG'
chr_df4$type='3MSD'

chr_df=rbind(chr_df2,chr_df1,chr_df3, chr_df4)

#library(tidyr)
sum(unlist(chr_df[chr_df$type=="ancestral",'proportion']))
sum(unlist(chr_df[chr_df$type=="AG_specific",'proportion']))
sum(unlist(chr_df[chr_df$type=="denovo",'proportion']))
sum(unlist(chr_df[chr_df$type=="unann",'proportion']))
chr_df=chr_df[!chr_df$chromosome=='mitochondrion_genome',]

graph_chrom_counts=ggplot(chr_df, aes(x=chromosome, y=proportion1, fill=type))+
  geom_bar(position='dodge',stat='identity', aes(fill=type), color='black', size=0.25, width=0.75, alpha=0.8)+
  scale_fill_manual(values=c(  "#999999","#56B4E9", '#C75127',"#E69F00"), 
                     labels = c("AN (-AG)", "AG","MSD", 'MOD'))+
  #scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(limits = c(0,0.45), breaks=seq(0,.4,.05),expand = c(0, 0))+
  theme_classic()+
  ylab('gene proportion')+
  xlab("")+
  labs(fill=" ")+
  theme(text=element_text(size=15),
        legend.position = c(0.88,0.85))

#X chromosome only
graph_x_counts=ggplot(chr_df[chr_df$chromosome=='X',], aes(x=type, y=proportion1, fill=type))+
  geom_bar(position='dodge',stat='identity', aes(fill=type), color='black', size=0.25, width=0.75, alpha=0.85)+
  scale_fill_manual(values=c(  "#999999","#56B4E9", '#C75127',"#E69F00"), 
                    labels = c("AN", "AG","OU", 'YU'))+
  scale_x_discrete(labels=c("AN (-AG)", "AG","MSD", 'MOD'))+
  scale_y_continuous(limits = c(0,0.23), breaks=seq(0,.4,.05),expand = c(0, 0))+
  theme_classic()+
  ylab('proportion on X chromosome')+
  xlab("")+
  labs(fill=" ")+
  theme(text=element_text(size=25),
        legend.position = "none")+
  ggsignif::geom_signif(annotations = c(formatC("***"),
                        formatC("*")), 
                        y_position = c(.174,.2), xmin=c(1, 1), xmax=c(2, 3), textsize = 8)
graph_x_counts
library(rstatix)
chr_df_x=chr_df[chr_df$chromosome=='X',]
chr_df_x$proportion1
chr_df_x$proportion
chr_df_x$total<-(unlist(chr_df_x$proportion)/chr_df_x$proportion1)-unlist(chr_df_x$proportion)

chr_df_x$total

chr_df_fisher <- data.frame(notonX = chr_df_x$total, onX =unlist(chr_df_x$proportion))
row.names(chr_df_fisher)=c("MOD", "AN", "AG", "MSU")
pairwise_fisher_test(chr_df_fisher)

####distance to TSS########################################################################
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)

gr_allgenes = import.gff('all.gtf')

#reduce all genes to one range
gr_allgenes=unlist(reduce(split(gr_allgenes, gr_allgenes$gene_id)))

#populate df from granges objest
gene_TSS_df=data.frame(1:length(gr_allgenes))
gene_TSS_df$chr=as.character(seqnames(gr_allgenes))
gene_TSS_df$start1=(start(gr_allgenes))
gene_TSS_df$end1=(end(gr_allgenes))
gene_TSS_df$geneid=(names(gr_allgenes))

gene_TSS_df$strand=(as.character(strand(gr_allgenes)))

gene_TSS_df[which(gene_TSS_df$geneid=="FBgn0002781"),'strand']='-'
gene_TSS_df=gene_TSS_df[,-1]
gene_TSS_df$TSS=NA

gene_TSS_df=gene_TSS_df[gene_TSS_df$chr%in%c("2R", "3R", "3L", "2L", "X"),]

gene_TSS_df$TSS[which(gene_TSS_df$strand=="+")]=gene_TSS_df[which(gene_TSS_df$strand=="+"),'start1']
gene_TSS_df$TSS[which(gene_TSS_df$strand=="-")]=gene_TSS_df[which(gene_TSS_df$strand=="-"),'end1']
#read gene classes
gene_TSS_df$geneclass=NA 

######looking at the enriched, specific genes here
gene_TSS_df$geneclass='an'
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%final_screen_list[,1]]="dn"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%unan_out_list[,1]]="unanc"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%AG_specific_list[,1]]="anag"

EN_DN_list=read.delim(file ='EN_DN_list', sep='\t',header=F)
EN_UA_list=read.delim(file='EN_UA_list', sep='\t',header=F)
nonEN_DN_list=read.delim(file ='nonEN_DN_list', sep='\t',header=F)
nonEN_UA_list=read.delim(file='nonEN_UA_list', sep='\t',header=F)


gene_TSS_df$geneclass[gene_TSS_df$geneid%in%rownames(RNA_counts)[which(RNA_counts$genetype==1)]]='an'
gene_TSS_df=gene_TSS_df[!is.na(gene_TSS_df$geneclass),]

#distancing

gene_TSS_df$dist_to_an_plus=NA
gene_TSS_df$dist_to_an_minus=NA
gene_TSS_df$min_an=NA

###annotated-works
i=1
for(i in 1:nrow(gene_TSS_df)){
  temp_df=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                              gene_TSS_df$strand==gene_TSS_df$strand[i]&
                              gene_TSS_df$geneclass=='an'),]
  temp_df=temp_df[!(temp_df$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_an_plus[i]=min(abs(temp_df$TSS-gene_TSS_df$TSS[i]))
  
  temp_df1=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                              ((gene_TSS_df$strand)!=(gene_TSS_df$strand[i]))&
                              gene_TSS_df$geneclass=='an'),]
  temp_df1=temp_df1[!(temp_df1$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_an_minus[i]=min(abs(temp_df1$TSS-gene_TSS_df$TSS[i]))
  gene_TSS_df$min_an[i]=min(c(gene_TSS_df$dist_to_an_minus[i],gene_TSS_df$dist_to_an_plus[i]))
  #print(i)
}
####

gene_TSS_df$dist_to_anag_plus=NA
gene_TSS_df$dist_to_anag_minus=NA
gene_TSS_df$min_anag=NA
for(i in 1:nrow(gene_TSS_df)){
  temp_df=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                              gene_TSS_df$strand==gene_TSS_df$strand[i]&
                              gene_TSS_df$geneclass=='anag'),]
  temp_df=temp_df[!(temp_df$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_anag_plus[i]=min(abs(temp_df$TSS-gene_TSS_df$TSS[i]))
  
  temp_df1=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                               ((gene_TSS_df$strand)!=(gene_TSS_df$strand[i]))&
                               gene_TSS_df$geneclass=='anag'),]
  temp_df1=temp_df1[!(temp_df1$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_anag_minus[i]=min(abs(temp_df1$TSS-gene_TSS_df$TSS[i]))
  gene_TSS_df$min_anag[i]=min(c(gene_TSS_df$dist_to_anag_minus[i],gene_TSS_df$dist_to_anag_plus[i]))
  #print(i)
}

####

gene_TSS_df$dist_to_ua_plus=NA
gene_TSS_df$dist_to_ua_minus=NA
gene_TSS_df$min_ua=NA

for(i in 1:nrow(gene_TSS_df)){
  temp_df=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                              gene_TSS_df$strand==gene_TSS_df$strand[i]&
                              gene_TSS_df$geneclass=='unanc'),]
  temp_df=temp_df[!(temp_df$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_ua_plus[i]=min(abs(temp_df$TSS-gene_TSS_df$TSS[i]))
  
  temp_df1=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                               ((gene_TSS_df$strand)!=(gene_TSS_df$strand[i]))&
                               gene_TSS_df$geneclass=='unanc'),]
  temp_df1=temp_df1[!(temp_df1$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_ua_minus[i]=min(abs(temp_df1$TSS-gene_TSS_df$TSS[i]))
  
  gene_TSS_df$min_ua[i]=min(c(gene_TSS_df$dist_to_ua_minus[i],gene_TSS_df$dist_to_ua_plus[i]))
  #print(i)
}


##################

gene_TSS_df$min_to_AG=NA
for(i in 1:nrow(gene_TSS_df)){
  gene_TSS_df$min_to_AG[i]=min(c(gene_TSS_df$dist_to_anag_minus[i], gene_TSS_df$dist_to_anag_plus[i]))
}
gene_TSS_df$group="Annotated"
gene_TSS_df[gene_TSS_df$geneid%in%unan_out_list$MSD_list1,'group']="UA"
gene_TSS_df[gene_TSS_df$geneid%in%final_screen_list[,1],'group']="DN"
unique(gene_TSS_df$geneclass)
gene_TSS_df$enrichment="BE"
gene_TSS_df[gene_TSS_df$geneclass%in%c('anag', 'ENdn', 'UAdn'), 'enrichment']='AG en'

gene_TSS_df1=gene_TSS_df[gene_TSS_df$geneid%in%rownames(RNA_counts)[RNA_counts$num_greater1>0],]


ggplot(gene_TSS_df1, aes(x=group, y=min_to_AG+1, fill=enrichment))+
  geom_boxplot(size=1)+
  scale_y_log10(limits=c(1,10000000))+
  theme_classic()+
  scale_fill_manual(values = c("red3", "white"), 
                    labels=c('AG/En+', 'All Other'))+
  scale_x_discrete(labels=c('Annotated', 'MOD', 'MSD'))+
  ylab('Nearest AG neighbor (bp)')+
  theme(text =element_text(size=25))

wilcox.test(gene_TSS_df1[gene_TSS_df1$group=='Annotated',][gene_TSS_df1$enrichment=='BE', 'min_to_AG'],
            gene_TSS_df1[gene_TSS_df1$group=='Annotated',][gene_TSS_df1$enrichment=='AG en', 'min_to_AG'])

wilcox.test(gene_TSS_df1[gene_TSS_df1$group=='DN',][gene_TSS_df1$enrichment=='BE', 'min_to_AG'],
            gene_TSS_df1[gene_TSS_df1$group=='DN',][gene_TSS_df1$enrichment=='AG en', 'min_to_AG'])

wilcox.test(gene_TSS_df1[gene_TSS_df1$group=='UA',][gene_TSS_df1$enrichment=='BE', 'min_to_AG'],
            gene_TSS_df1[gene_TSS_df1$group=='UA',][gene_TSS_df1$enrichment=='AG en', 'min_to_AG'])

wilcox.test(gene_TSS_df1[gene_TSS_df1$group=='UA','min_to_AG'],
            gene_TSS_df1[gene_TSS_df1$group=='DN','min_to_AG'])


gene_TSS_df$min_to_AN=NA
for(i in 1:nrow(gene_TSS_df)){
  gene_TSS_df$min_to_AN[i]=min(c(gene_TSS_df$dist_to_an_minus[i], gene_TSS_df$dist_to_an_plus[i]))
}

ggplot(gene_TSS_df, aes(x=group, y=min_to_AN+1, fill=enrichment))+
  geom_boxplot(size=1)+
  scale_y_log10(limits=c(1,10000000))+
  theme_classic()+
  scale_fill_manual(values = c("red3", "white"), 
                    labels=c('AG/En+', 'All Other'))+
  scale_x_discrete(labels=c('Annotated', 'MOD', 'MSD'))+
  ylab('Nearest AG neighbor (bp)')+
  theme(text =element_text(size=25))
ggplot(gene_TSS_df, aes(x=geneclass, y=min_to_AN+1))+
  geom_boxplot()+
  scale_y_log10(limits=c(1,10000000))



ggplot(gene_TSS_df, aes(x=enrichment, fill=chr))+
  facet_wrap(~group, scales = 'free')+
  geom_bar(stat="count", position="dodge")+
  theme_minimal()




#################



###test10/10/23
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%final_screen_list[,1]]="ann"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%final_screen_list[,1]]="dn"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%unan_out_list[,1]]="unanc"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%unan_out_list[,1]]="unanc"


gene_TSS_df$dist_to_dn_plus=NA
gene_TSS_df$dist_to_dn_minus=NA
gene_TSS_df$min_dn=NA

for(i in 1:nrow(gene_TSS_df)){
  temp_df=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                              gene_TSS_df$strand==gene_TSS_df$strand[i]&
                              gene_TSS_df$geneclass=='dn'),]
  temp_df=temp_df[!(temp_df$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_dn_plus[i]=min(abs(temp_df$TSS-gene_TSS_df$TSS[i]))
  
  temp_df1=gene_TSS_df[which(gene_TSS_df$chr==gene_TSS_df$chr[i]&
                               ((gene_TSS_df$strand)!=(gene_TSS_df$strand[i]))&
                               gene_TSS_df$geneclass=='dn'),]
  temp_df1=temp_df1[!(temp_df1$geneid%in%gene_TSS_df$geneid[i]),]
  gene_TSS_df$dist_to_dn_minus[i]=min(abs(temp_df1$TSS-gene_TSS_df$TSS[i]))
  
  gene_TSS_df$min_dn[i]=min(c(gene_TSS_df$dist_to_dn_minus[i],gene_TSS_df$dist_to_dn_plus[i]))
  #print(i)
}

ggplot(gene_TSS_df, aes(x=geneclass, y=dist_to_dn_minus+1))+
  geom_boxplot()+
  scale_y_log10(limits=c(1,10000000))

ggplot(gene_TSS_df, aes(x=geneclass, y=dist_to_dn_plus+1))+
  geom_boxplot()+
  scale_y_log10(limits=c(1,10000000))

ggplot(gene_TSS_df, aes(x=geneclass, y=min_dist+1))+
  geom_boxplot()+
  scale_y_log10(limits=c(1,10000000))

a=gene_TSS_df[,c(7:ncol(gene_TSS_df))]
gene_TSS_df$group=NULL
gene_TSS_df$enrichment=NULL

gene_TSS_df1=gene_TSS_df[,c(7:ncol(gene_TSS_df))] %>%
  pivot_longer(!geneclass, names_to = "gene_class", values_to = "count")

gene_TSS_df1$direction='same orientation'
gene_TSS_df1$direction[grep("minus", gene_TSS_df1$gene_class)]='opposite orientation'
gene_TSS_df1$direction[grep("min_", gene_TSS_df1$gene_class)]='minimum'


gene_TSS_df1$gene_class=gsub("_plus", "", gene_TSS_df1$gene_class)
gene_TSS_df1$gene_class=gsub("_minus", "", gene_TSS_df1$gene_class)


all_distance=ggplot(gene_TSS_df1, aes(fill=geneclass, y=log(as.numeric(count)+1), x=gene_class))+
  geom_boxplot(outlier.size = 0.1)+
  #scale_y_log10(limits=c(1,10000000))+
  facet_grid(cols = vars(direction))+
  labs(fill="2nd gene class")+
  scale_x_discrete(labels=c('An; Exp+/- Sp-', 'An; Exp+ Sp+', "Un; DN"))+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00" ,"#C75127"), 
                    labels=c('An; Exp+/- Sp-', 'An; Exp+ Sp+', "Un; DN", "Un; Anc"))+
  xlab('1st gene class')+
  ylab('distance from TSS to nearest TSS + 1 (bp)')+
  theme_bw()



gene_TSS_df1_plot=gene_TSS_df1[gene_TSS_df1$direction=="minimum",]
unique(gene_TSS_df1_plot$gene_class)
gene_TSS_df1_plot=gene_TSS_df1_plot[!gene_TSS_df1_plot$gene_class=='min_an', ]

gene_TSS_df1_plot=droplevels(gene_TSS_df1_plot)
unique(gene_TSS_df1_plot$gene_class)
gene_TSS_df1_plot$gene_class[1]

gene_TSS_df1_plot[gene_TSS_df1_plot$gene_class%in%c("min_anag","min_ua","min_dn"),]
nearest_neighbor=ggplot(gene_TSS_df1_plot[gene_TSS_df1_plot$gene_class%in%c("min_anag","min_ua","min_dn"),], aes(fill=geneclass, y=(count+1), x=gene_class))+
  geom_boxplot(outlier.size = 0.1)+
  #scale_y_log10(limits=c(1,10000000))+
  #facet_grid(cols = vars(direction))+
  labs(fill=" ")+
  scale_x_discrete(labels=c('AG', "MOD","MSD"))+
  scale_y_log10()+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00" ,"#C75127"), #grey or purple? #999999#7A65A5
                    labels=c('AN (-AG)', 'AG', "MOD", "MSD"))+
  xlab(' ')+
  ylab('TSS distance + 1')+
  theme_bw()+
  theme(text=element_text(size=25))
options(scipen = 999999999)
nearest_neighbor

#distance by expression controls
gene_TSS_df_AG=gene_TSS_df[gene_TSS_df$geneid%in%AG_specific_list[,],]
gene_TSS_df_AG$exp=RNA_counts[row.names(RNA_counts)%in%gene_TSS_df_AG$geneid,'mean']

gene_TSS_df_AG$minanag=pmax(gene_TSS_df_AG$dist_to_anag_plus, gene_TSS_df_AG$dist_to_anag_minus)
AG_exVdistance_graph=ggplot(gene_TSS_df_AG, aes(x=min_to_AG, y=exp))+
  geom_point()+
  scale_x_log10(limits=c(1,10000000), breaks=c(1,100,10000,1000000))+
  scale_y_log10(limits=c(0.1,10000000), breaks=c(1,100,10000,1000000))+
  xlab('closest AG gene')+
  ylab('mean expression (TPM)')+
  geom_smooth(method = "lm", se=FALSE, color="black", linetype='dashed') +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  theme_classic()+
  theme(text=element_text(size=25))+
  ggtitle('AG genes')

#MOD
gene_TSS_df_MOD=gene_TSS_df[gene_TSS_df$geneid%in%final_screen_list[,1],]
gene_TSS_df_MOD$exp=RNA_counts[row.names(RNA_counts)%in%gene_TSS_df_MOD$geneid,'mean']
MOD_expVdistance_graph=ggplot(gene_TSS_df_MOD, aes(x=min_to_AG, y=exp))+
  geom_point()+
  scale_x_log10(limits=c(1,10000000), breaks=c(1,100,10000,1000000))+
  scale_y_log10(limits=c(0.01,100), breaks=c(.1,1,10,100))+
  xlab('closest AG gene')+
  ylab('mean expression (TPM)')+
  geom_smooth(method = "lm", se=FALSE, color="black", linetype='dashed') +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  theme_classic()+
  theme(text=element_text(size=25))+
  ggtitle('MOD genes')

gene_TSS_df_MSD=gene_TSS_df[gene_TSS_df$geneid%in%unan_out_list[,1],]
gene_TSS_df_MSD$exp=RNA_counts[row.names(RNA_counts)%in%gene_TSS_df_MSD$geneid,'mean']
MSD_expVdistance_graph=ggplot(gene_TSS_df_MSD, aes(x=min_to_AG, y=exp))+
  geom_point()+
  scale_x_log10(limits=c(1,10000000), breaks=c(1,100,10000,1000000))+
  scale_y_log10(limits=c(0.01,100), breaks=c(.1,1,10,100))+
  xlab('closest AG gene')+
  ylab('mean expression (TPM)')+
  geom_smooth(method = "lm", se=FALSE, color="black", linetype='dashed') +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  theme_classic()+
  theme(text=element_text(size=25))+
  ggtitle('MSD genes')


#regular genes
ggplot(gene_TSS_df_MSD, aes(x=min_to_AN, y=exp))+
  geom_point()+
  scale_x_log10(limits=c(1,10000000), breaks=c(1,100,10000,1000000))+
  scale_y_log10(limits=c(0.01,100), breaks=c(.1,1,10,100))+
  xlab('closest annotated gene')+
  ylab('mean expression (TPM)')+
  geom_smooth(method = "lm", se=FALSE, color="black", linetype='dashed') +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  theme_classic()+
  theme(text=element_text(size=25))+
  ggtitle('MSD genes')

ggplot(gene_TSS_df_AG, aes(x=min_to_AN+1, y=exp))+
  geom_point()+
  scale_x_log10(limits=c(1,10000000), breaks=c(1,100,10000,1000000))+
  scale_y_log10(limits=c(0.1,10000000), breaks=c(1,100,10000,1000000))+
  xlab('closest annotated gene')+
  ylab('mean expression (TPM)')+
  geom_smooth(method = "lm", se=FALSE, color="black", linetype='dashed') +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  theme_classic()+
  theme(text=element_text(size=25))+
  ggtitle('AG genes')

ggplot(gene_TSS_df_MOD, aes(x=min_to_AN, y=exp))+
  geom_point()+
  scale_x_log10(limits=c(1,10000000), breaks=c(1,100,10000,1000000))+
  scale_y_log10(limits=c(0.01,100), breaks=c(.1,1,10,100))+
  xlab('closest annotated gene')+
  ylab('mean expression (TPM)')+
  geom_smooth(method = "lm", se=FALSE, color="black", linetype='dashed') +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  theme_classic()+
  theme(text=element_text(size=25))+
  ggtitle('MOD genes')

#not significant correlation between exp and nearest neighbor
cor.test(log(gene_TSS_df_MOD$min_to_AG), log(gene_TSS_df_MOD$exp), method='pearson')
cor.test(log(gene_TSS_df_MSD$min_to_AG), log(gene_TSS_df_MSD$exp), method='pearson')
cor.test(log(gene_TSS_df_AG$min_to_AG), log(gene_TSS_df_AG$exp), method='pearson')

gene_TSS_df_tests=gene_TSS_df1[gene_TSS_df1$direction=='opposite orientation',]
gene_TSS_df_tests=as.data.frame(gene_TSS_df_tests)

anova_test=(aov(log(count+1)~geneclass+gene_class, data = gene_TSS_df_tests ))
TukeyHSD(anova_test)
pairwise.wilcox.test(gene_TSS_df_tests$count, gene_TSS_df_tests$geneclass, p.adjust.method = p.adjust.methods,
                     paired = FALSE)

#######length of de novo genes########################################################################
#from https://www.biostars.org/p/83901/
# First, import the GTF-file that you have also used as input for htseq-count

dn_unfiltered_gtf=read.delim(file='all_unfilteredtaco_1tpm.gtf', sep='\t',
           quote = "", header=F)
dmel_anc_gtf=read.delim(file='6.41_all.gtf', sep='\t',
                        quote = "", header=F)

txdb <- makeTxDbFromGFF("all.gtf",format="gtf")

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
gene_length_df = as.data.frame(sum(width(reduce(exons.list.per.gene))))

colnames(gene_length_df)='length'
gene_length_df$type='1_ancestral'

gene_length_df[row.names(gene_length_df)%in%final_screen_list[,1],'type']='4_MOD'
gene_length_df[row.names(gene_length_df)%in%AG_specific_list[,1],'type']='2_AG'
gene_length_df[row.names(gene_length_df)%in%unan_out_list[,1],'type']='3_MSD'



#significance values for exon length
pairwise.wilcox.test(gene_length_df$length, gene_length_df$type)

exonlength_no_sig_graph=ggplot(gene_length_df, aes(y=length, x=type, fill=type))+
  geom_boxplot(aes(fill=type), alpha=0.8,outlier.alpha = 1, width=0.6)+
  scale_fill_manual(values=c("#56B4E9", "#999999","#C75127", "#E69F00" ))+
  scale_x_discrete(labels=c("AN(-AG)", "AG", "MSD", "MOD"))+
  scale_colour_manual()+
  #scale_y_log10(limits=c(20,120000))+
  theme_classic()+
  ylab('Total exon length (bp)')+
  scale_x_discrete(labels = c("AN", "AG","MSD", "MOD"))+
  theme(text = element_text(size=25),
        legend.position = "none")

exonlength_sig_graph=ggplot(gene_length_df, aes(y=length, x=type, fill=type))+
  geom_boxplot(aes(fill=type), alpha=0.8,outlier.alpha = 1, width=0.6)+
  scale_fill_manual(values=c("#999999", "#56B4E9","#C75127", "#E69F00" ))+
  scale_x_discrete(labels=c("AN(-AG)", "AG", "MSD", "MOD"))+
  scale_colour_manual()+
  #scale_y_continuous(limits=c(2.5,17.4), breaks=c(3,5,7,9,11,13,15))+
  scale_y_log10(limits=c(10, 200000), breaks =c(10,100,1000,10000,100000))+
  theme_classic()+
  ylab('Total exon length (bp)')+
  ggsignif::geom_signif(annotations = c(formatC("bcd"),
                              formatC("ad"),
                              formatC("ad"),
                              formatC("abc")), 
              y_position = c(5.2), xmin=c(1, 2, 3, 4), xmax=c(1, 2, 3, 4), textsize = 6)+  theme(text = element_text(size=25),
        legend.position = "none")


exonlength_sig_graph

#############exon counts####################
exonscounts_df=data.frame(table(data.frame(exons.list.per.gene)$group_name))
rownames(exonscounts_df)=exonscounts_df[,1]
exonscounts_df$Var1=NULL
exonscounts_df$exon_count=exonscounts_df$Freq
exonscounts_df$Freq=NULL

exonscounts_df$type1='ancestral'
exonscounts_df[row.names(exonscounts_df)%in%final_screen_list[,1],'type1']='denovo'
exonscounts_df[row.names(exonscounts_df)%in%unan_out_list[,1],'type1']='unannotated'
exonscounts_df[row.names(exonscounts_df)%in%AG_specific_list[,1],'type1']='AG'


exon_dn=as.data.frame(table(exonscounts_df[exonscounts_df$type1=='denovo','exon_count']))
exon_AG=as.data.frame(table(exonscounts_df[exonscounts_df$type1=='AG','exon_count']))
exon_anc=as.data.frame(table(exonscounts_df[exonscounts_df$type1=='ancestral','exon_count']))
exon_unann=as.data.frame(table(exonscounts_df[exonscounts_df$type1=='unannotated','exon_count']))

exon_AG$Freq=exon_AG$Freq/length(which(exonscounts_df$type1=='AG'))
exon_dn$Freq=exon_dn$Freq/length(which(exonscounts_df$type1=='denovo'))
exon_anc$Freq=exon_anc$Freq/length(which(exonscounts_df$type1=='ancestral'))
exon_unann$Freq=exon_unann$Freq/length(which(exonscounts_df$type1=='unannotated'))
exon_AG$type1='AG'
exon_dn$type1='dn'
exon_anc$type1='anc'
exon_unann$type1="unann"

exonscounts_df1=exonscounts_df %>% count(exon_count, type1)

#wilcoxon tests for exons
str(exonscounts_df1$n)
exonscounts_df[exonscounts_df$type1=="AG",1]

#AG > MOD ***
wilcox.test(exonscounts_df[exonscounts_df$type1=="AG",1], exonscounts_df[exonscounts_df$type1=="denovo",1])
#AG = MSD NS
wilcox.test(exonscounts_df[exonscounts_df$type1=="AG",1], exonscounts_df[exonscounts_df$type1=="unannotated",1])
#MOD  MSD *
wilcox.test(exonscounts_df[exonscounts_df$type1=="denovo",1], exonscounts_df[exonscounts_df$type1=="unannotated",1])

#making a greater than 15 
ancestral_16=sum(exonscounts_df1[which(exonscounts_df1$type1=='ancestral'&exonscounts_df1$exon_count>7),'n'])
AG_16=sum(exonscounts_df1[which(exonscounts_df1$type1=='AG'&exonscounts_df1$exon_count>7),'n'])
denovo_16=sum(exonscounts_df1[which(exonscounts_df1$type1=='denovo'&exonscounts_df1$exon_count>7),'n'])
unann_16=sum(exonscounts_df1[which(exonscounts_df1$type1=='unannotated'&exonscounts_df1$exon_count>7),'n'])


exonscounts_df1=exonscounts_df1[-which(exonscounts_df1$type1=='ancestral'&exonscounts_df1$exon_count>7),]
exonscounts_df1=exonscounts_df1[-which(exonscounts_df1$type1=='AG'&exonscounts_df1$exon_count>7),]
exonscounts_df1=exonscounts_df1[-which(exonscounts_df1$type1=='denovo'&exonscounts_df1$exon_count>7),]
exonscounts_df1=exonscounts_df1[-which(exonscounts_df1$type1=='unannotated'&exonscounts_df1$exon_count>7),]


exonscounts_df1=rbind(exonscounts_df1, c(8,'ancestral',ancestral_16))
exonscounts_df1=rbind(exonscounts_df1, c(8,'AG',AG_16))
exonscounts_df1=rbind(exonscounts_df1, c(8,'denovo',denovo_16))
exonscounts_df1=rbind(exonscounts_df1, c(8,'unannotated',denovo_16))


exonscounts_df1$n=as.numeric(exonscounts_df1$n)
exonscounts_df1$exon_count=as.numeric(exonscounts_df1$exon_count)

exonscounts_df1[which(exonscounts_df1$type1=='AG'),'n']=
  exonscounts_df1[which(exonscounts_df1$type1=='AG'),'n']/sum(exonscounts_df1[which(exonscounts_df1$type1=='AG'),'n'])
exonscounts_df1[which(exonscounts_df1$type1=='ancestral'),'n']=
  exonscounts_df1[which(exonscounts_df1$type1=='ancestral'),'n']/sum(exonscounts_df1[which(exonscounts_df1$type1=='ancestral'),'n'])
exonscounts_df1[which(exonscounts_df1$type1=='denovo'),'n']=
  exonscounts_df1[which(exonscounts_df1$type1=='denovo'),'n']/sum(exonscounts_df1[which(exonscounts_df1$type1=='denovo'),'n'])
exonscounts_df1[which(exonscounts_df1$type1=='unannotated'),'n']=
  exonscounts_df1[which(exonscounts_df1$type1=='unannotated'),'n']/sum(exonscounts_df1[which(exonscounts_df1$type1=='unannotated'),'n'])


#re ordering and renaming plot
exonscounts_df1$exon_count[exonscounts_df1$exon_count=='8']='8+'
exonscounts_df1$type1=gsub("ancestral", "1ancestral", exonscounts_df1$type1)
exonscounts_df1$type1
graph_exon_count<-ggplot(exonscounts_df1, aes(y=n, x=exon_count, fill=type1))+
  geom_bar(position = position_dodge(preserve = "single"),stat='identity', aes(fill=type1), color='black',
           size=0.25, width=0.75, alpha=0.85)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))+
  xlab('Number of exons')+
  ylab("proportion")+
  #scale_x_continuous(limits = c(0.1,8.5), breaks=seq(0,8,1), expand = c(0, 0))+ #is a factor now with "8+"
  scale_y_continuous(limits = c(0,0.62), breaks=seq(0,.8,0.1), expand = c(0, 0))+
  scale_fill_manual(values=c("#999999", "#56B4E9", "#E69F00", '#C75127'),labels = c("AN (-AG)", "AG", "MOD", 'MSD'))+
  theme(legend.position = c(0.88, 0.83),
        text=element_text(size=25))






###################Simulans rarity and unannotated rarity###################

unann_list=read.delim(file='tables/MOD_list.txt', sep='\t', header = F)

unann_list_nogenome=read.delim(file='tables/MSD_list.txt', sep='\t', header = T)

simfixed_genelist=read.delim(file='tables/genelist_simfixed', sep='\t', header = F)
sim1genelist=read.delim(file='tables/gene_list_sim1', sep='\t', header = F)
sim2genelist=read.delim(file='tables/gene_list_sim2', sep='\t', header = F)
sim3genelist=read.delim(file='tables/gene_list_sim3', sep='\t', header = F)

RNA_counts_unann=RNA_counts[rownames(RNA_counts)%in%unan_out_list[,1],]

RNA_counts_unann$genetype='Unannotated'
names(RNA_counts_dn)
names(RNA_counts_unann)[37]='diptest'


#
unan_vs_dn_df=RNA_counts_unann
unan_vs_dn_df$sim_fixed=NA

unan_vs_dn_df[(rownames(unan_vs_dn_df)[unan_vs_dn_df$genetype=='Unannotated']%in%simfixed_genelist[,1]),'sim_fixed']='yes'
unan_vs_dn_df[(rownames(unan_vs_dn_df)[unan_vs_dn_df$genetype=='Unannotated']%in%sim1genelist[,1]),'sim_fixed']='1_line'
unan_vs_dn_df[(rownames(unan_vs_dn_df)[unan_vs_dn_df$genetype=='Unannotated']%in%sim2genelist[,1]),'sim_fixed']='2_line'
unan_vs_dn_df[(rownames(unan_vs_dn_df)[unan_vs_dn_df$genetype=='Unannotated']%in%sim3genelist[,1]),'sim_fixed']='3_line'

unan_vs_dn_df=unan_vs_dn_df[-which(is.na(unan_vs_dn_df$sim_fixed)==T),]

##simulans specific
unann_sim_df=unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated',]
unann_sim_df
unann_sim_df$num_greater1=apply(unann_sim_df[,c(1:29)], 1, function(x) length(which(x>1)))



unann_sim_df$sim_fixed=gsub("3_line", 3, unann_sim_df$sim_fixed)
unann_sim_df$sim_fixed=gsub("2_line", 2, unann_sim_df$sim_fixed)
unann_sim_df$sim_fixed=gsub("1_line", 1, unann_sim_df$sim_fixed)

pairwise.wilcox.test(unann_sim_df$num_greater1, as.numeric(unann_sim_df$sim_fixed))
#library(ggsignif)


ggplot(unann_sim_df, aes(y=num_greater1, x=sim_fixed))+
  geom_boxplot(fill= 'grey')+
  theme_classic()+
  ylab("# mel DGRP (>1 TPM)")+
  xlab("# sim (BLAST match)")+
  scale_x_discrete(labels=c('1/3','2/3','3/3'))+
  scale_y_continuous(limits=c(0,35))+
  geom_signif(annotations = c(formatC("*"),
                              formatC("**"),
                              formatC("NS")), 
              y_position = c(29, 33.5, 31), xmin=c(1, 1, 2), xmax=c(2, 3, 3), textsize = 6)+
  theme(text = element_text(size=25))

length(which(unann_sim_df$sim_fixed==3))/nrow(unann_sim_df)
  
nrow(unann_sim_df)
##EB6100","#002C60","#31446D","#686B76","#C1B36F", "#F9E14B"
unan_vs_dn_df


RNA_counts_dn$num_greater1=apply(RNA_counts_dn[,c(1:29)], 1, function(x) length(which(x>1)))
RNA_counts_unann$num_greater1=apply(RNA_counts_unann[,c(1:29)], 1, function(x) length(which(x>1)))


unan_vs_dn_df=rbind(RNA_counts_unann[,c(1:36)], RNA_counts_dn[,c(1:36)])
unan_vs_dn_df$num_greater1=apply(unan_vs_dn_df[,c(1:29)], 1, function(x) length(which(x>1)))

#de novo 
MOD_prop_graph=ggplot(unan_vs_dn_df[unan_vs_dn_df$genetype=='4',], aes(x=as.numeric(num_greater1), fill=genetype))+
  geom_histogram(color= 'black',  
                 position='stack', alpha=1, 
                 #position='identity', alpha=0.75,
                 bins=30)+
  theme_classic()+
  guides(fill=guide_legend(" "))+
  scale_fill_manual(values=c('goldenrod'), labels=c('de novo'))+
  scale_y_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5), limits=c(0,58))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  xlab("#lines expressed >1 TPM")+
  ylab('MOD count')+
  theme(legend.position='none',
        text = element_text(size=25))+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               y = 20, yend = 15,arrow = arrow(length = unit(0.04, "npc")), color="goldenrod", size=1.5)
#MSD
MSD_prop_graph=ggplot(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated',], aes(x=as.numeric(num_greater1), fill=genetype))+
  geom_histogram(color= 'black',  
                 position='stack', alpha=.8, 
                 #position='identity', alpha=0.75,
                 bins=30)+
  theme_classic()+
  guides(fill=guide_legend(" "))+
  scale_fill_manual(values=c('#C75127'), labels=c('de novo'))+
  scale_y_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5), limits=c(0,58))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  xlab("#lines expressed > 1tpm")+
  ylab('MSD count')+
  theme(legend.position='none',
        text = element_text(size=25))+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               y = 20, yend = 15,arrow = arrow(length = unit(0.04, "npc")), color="#C75127", size=1.5)

#sig lower proportion of MOD genes
wilcox.test(unan_vs_dn_df[unan_vs_dn_df$genetype=="Unannotated","num_greater1"],
            unan_vs_dn_df[unan_vs_dn_df$genetype=="4","num_greater1"])
            

#code to make rarefaction graphs
#uses some very inelegant for loops
RNA_counts_dn1=RNA_counts_dn[,c(1:29)]

temp_df=RNA_counts_dn1
rarefaction=as.data.frame(c(1:29))
rarefaction$count=NA

#count MOD genes first
for(i in 1:28){
  tmp_col=which.max(apply(temp_df, 2, function(x) length(which(x>1))))
  tmp_dn_genes=which(temp_df[,tmp_col]>1)
  
  rarefaction$count[i]=length(tmp_dn_genes)
  
  temp_df=temp_df[-tmp_dn_genes,]
  temp_df=temp_df[,-tmp_col]
  print(i)
}

RNA_counts_DNunann=rbind(RNA_counts[rownames(RNA_counts)%in%unan_out_list[,1],],
                         RNA_counts[rownames(RNA_counts)%in%final_screen_list$x,])
RNA_counts_DNunann=RNA_counts_DNunann[,c(1:29)]

temp_df=RNA_counts_DNunann
rarefaction$unann_count=NA

#count MSD genes
for(i in 1:28){
  tmp_col=which.max(apply(temp_df, 2, function(x) length(which(x>1))))
  tmp_dn_genes=which(temp_df[,tmp_col]>1)
  
  rarefaction$unann_count[i]=length(tmp_dn_genes)
  
  temp_df=temp_df[-tmp_dn_genes,]
  temp_df=temp_df[,-tmp_col]
  #print(i)
}

#set first point to zero
rarefaction$count[29]=0
rarefaction$unann_count[29]=0

#put counts together and populate df from there
rarefaction1=as.data.frame(c(cumsum(rarefaction$count),cumsum(rarefaction$unann_count)))
names(rarefaction1)='counts'

rarefaction1$type=c(rep('denovo', times=29), rep('1unann', times=29))
rarefaction1$num=c(rep(c(1:29), times=2))

rarefaction_graph=ggplot(rarefaction1, aes(x=num, y=counts, color=type)) + 
  geom_line(size=1.1) + 
  geom_point(size=1.8)+
  theme_classic()+
  scale_y_continuous(expand=c(0,0), 
                     limits=c(0,310),
                     breaks=pretty_breaks(n=10))+
  scale_x_continuous(limits=c(0,29),
                     breaks=pretty_breaks(n=10))+
  xlab('DGRP lines')+
  ylab('unannotated genes')+
  scale_color_manual(values=c('grey30','goldenrod'),
                     label=c('MOD+MSD','MOD'))+
  theme(text = element_text(size=25),
        legend.title = element_blank())

#
RNA_counts_dn$exon=exonscounts_df[(row.names(exonscounts_df)%in%row.names(RNA_counts_dn)),1]
RNA_counts_unann$exon=exonscounts_df[exonscounts_df$type1=='unannotated',1]

my.formula <- y ~ x #for linreg lines

Nexon_Nlines_MOD_graph=RNA_counts_dn%>%add_count(exon,num_greater1)%>%ggplot(aes(y=exon, x=num_greater1))+
  geom_point(aes(size=n))+
  theme_classic()+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula, linetype='solid') +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  scale_x_continuous(limits = c(0,31))+
  scale_y_continuous(limits=c(0.1,12), breaks = scales::pretty_breaks(n = 5))+
  ylab('# exons')+
  xlab('#lines expressed >1 TPM')+
  scale_size_continuous(range=c(1,5))+
  theme(legend.position = c(0.85,0.85),
        legend.title = element_blank(),
        text = element_text(size=30))+
  scale_size(range=c(1.2,8),breaks=c(0,1,5,10,20),
             labels=c(">=0",">=1",">=5",">=10",">=20"),guide="legend")

Nexon_Nlines_MSD_graph=RNA_counts_unann%>%add_count(exon,num_greater1)%>%ggplot(aes(y=exon, x=num_greater1))+
  geom_point(aes(size=n))+
  theme_classic()+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula, linetype='solid') +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  scale_x_continuous(limits = c(0,31))+
  scale_y_continuous(limits=c(0.1,12), breaks = scales::pretty_breaks(n = 5))+
  ylab('# exons')+
  xlab('#lines expressed >1 TPM')+
  scale_size_continuous(range=c(1,5))+
  theme(legend.position = c(0.85,0.85),
        legend.title = element_blank(),
        text = element_text(size=30))+
  scale_size(range=c(1.2,8),breaks=c(0,1,5,10,20),
             labels=c(">=0",">=1",">=5",">=10",">=20"),guide="legend")

##############
#
RNA_counts_unann1=RNA_counts_unann
RNA_counts_unann1$diptest=NULL
RNA_counts_unann1$exon=NULL

RNA_counts_unann1$type='un'
RNA_counts_dn$type='dn'
dn_v_unn_df=rbind(RNA_counts_dn, RNA_counts_unann1)

############
RNA_counts_dn$exonlen=gene_length_df[(row.names(gene_length_df)%in%row.names(RNA_counts_dn)),1]
RNA_counts_unann$exonlen=gene_length_df[(row.names(gene_length_df)%in%row.names(RNA_counts_unann)),1]

#sig for MSD but not MOD
cor.test(RNA_counts_unann$exon, RNA_counts_unann$num_greater1)
cor.test(RNA_counts_dn$exon, RNA_counts_dn$num_greater1)

exonlen_nlines_MOD_plot=ggplot(RNA_counts_dn, aes(y=exonlen, x=num_greater1))+
  geom_point(size=1,position=position_jitter(w=0.15))+
  theme_classic()+
  #scale_x_log10()+
  #scale_y_log10(limits=c(200,9000))+
  scale_y_continuous(limits=c(0,7900), breaks = c(0,1000,2000,3000,4000,5000,6000,7000))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  ylab('gene length (exons)')+
  xlab('#lines expressed >1 TPM')+
  ggtitle("de novo")

exonlen_nlines_MSD_plot=ggplot(RNA_counts_unann, aes(y=exonlen, x=num_greater1))+
  geom_point(size=2,position=position_jitter(w=0.15))+
  theme_classic()+
  #scale_x_log10()+
  scale_y_log10(limits=c(200,10000))+
  scale_x_continuous(limits=c(0,30))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 4))+
  ylab('exon length')+
  xlab('#lines expressed >1 TPM')+
  theme(text=element_text(size=30))

#not significant for either
cor.test(RNA_counts_dn$num_greater1, RNA_counts_dn$exonlen)
cor.test(RNA_counts_unann$num_greater1, RNA_counts_unann$exonlen)

#

maxexp_Nlines_MSD_plot=ggplot(RNA_counts_unann, aes(y=max, x=num_greater1))+
  geom_point(size=2)+
  theme_classic()+
  scale_y_log10(limits=c(1,220),breaks = c(1,3,10,30,100))+
  scale_x_continuous(limits=c(0,31))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  ylab('max expression')+
  xlab('#lines expressed >1 TPM')+
  theme(text=element_text(size=30))

a=lm(RNA_counts_unann$max~RNA_counts_unann$num_greater1)
summary(a)
b=lm(RNA_counts_dn$max~RNA_counts_dn$num_greater1)
summary(b)

ggplot(RNA_counts_dn, aes(y=max, x=num_greater1))+
  geom_point(size=2)+
  theme_classic()+
  scale_y_log10(limits=c(1,220),breaks = c(1,3,10,30,100))+
  scale_x_continuous(limits=c(0,31))+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE, size=8)+
  ylab('max expression')+
  xlab('#lines expressed >1 TPM')+
  theme(text=element_text(size=30))

cor.test(log(RNA_counts_dn$max), RNA_counts_dn$num_greater1)
cor.test(log(RNA_counts_unann$max), RNA_counts_unann$num_greater1)

#
chr_count_df
pairwise_fisher_test(chr_count_df[,3:4], p.adjust.method='holm')$p.adj

###############
TPM_threshold_plot=RNA_counts
TPM_threshold_plot1=TPM_threshold_plot[row.names(TPM_threshold_plot)%in%
                                      final_screen_list[,1],]
TPM_threshold_plot1=TPM_threshold_plot1[ order(TPM_threshold_plot1$max), ]
TPM_threshold_plot1$df=1:nrow(TPM_threshold_plot1)
TPM_threshold_plot1$df=abs(TPM_threshold_plot1$df-nrow(TPM_threshold_plot1))
TPM_threshold_plot1$class='MOD'

TPM_threshold_plot2=TPM_threshold_plot[row.names(TPM_threshold_plot)%in%
                                         unan_out_list[,1],]
TPM_threshold_plot2=TPM_threshold_plot2[ order(TPM_threshold_plot2$max), ]
TPM_threshold_plot2$df=1:nrow(TPM_threshold_plot2)
TPM_threshold_plot2$df=abs(TPM_threshold_plot2$df-nrow(TPM_threshold_plot2))
TPM_threshold_plot2$class='MSD'

graph_TPM_threshold=rbind(TPM_threshold_plot1,TPM_threshold_plot2)
colnames(TPM_threshold_plot)[37]

colnames(TPM_threshold_plot)[37]="counts1"

ggplot(TPM_threshold_plot, aes(y=counts1, x=max))+
         geom_line()


ggplot(graph_TPM_threshold, aes(y=counts1, x=max, color=class))+
  geom_line()+
  scale_x_log10(breaks = c(1,2,5,10,20,50,100))+
  scale_y_continuous(limits=c(0,150), breaks = scales::pretty_breaks(n = 10))+
  xlab('min TPM threshold')+
  ylab('#genes identified')+
  theme_classic()+
  scale_color_manual(values=c('goldenrod', "#C75127"), labels=c("MOD", "MSU"))+
  theme(text = element_text(size=25))

#Li's paper: candidate dn genes
Zhao_DN=c('FBgn0053664','FBgn0053669','FBgn0053668','FBgn0053667','FBgn0053666','FBgn0053665'
,'FBgn0260871','FBgn0085361','FBgn0040784','FBgn0260867','FBgn0264344','FBgn0051909')

Zhao_DN_exp=RNA_counts[row.names(RNA_counts)%in%Zhao_DN,]
write.table(Zhao_DN_exp, file='/Users/loganblair/Desktop/Zhao_DN.txt', sep='\t',
            quote=F)
apply(Zhao_DN_exp[,1:29], 1, max)

ggplot(Zhao_DN_exp, aes(x=max))+
  geom_histogram()

#filtering pipeline: potential DN AG genes
candidate_AG_DN=c('FBgn0053244','FBgn0053245','FBgn0259973',
'FBgn0261006','FBgn0263090','FBgn0263095','FBgn0265424',
'FBgn0265739','FBgn0265798','FBgn0265842','FBgn0265900',
'FBgn0267179','FBgn0267215','FBgn0267268','FBgn0267926')

AG_DN_final=c("FBgn0265424", "FBgn0265798", "FBgn0265842", "FBgn0267215")
  
AG_DN_candidates_final=RNA_counts[row.names(RNA_counts)%in%AG_DN_final,]
AG_DN_candidates_final=AG_DN_candidates_final[,1:29]
AG_DN_candidates_final$gene=row.names(AG_DN_candidates_final)
AG_DN_candidates_final=reshape2::melt(AG_DN_candidates_final, id.vars = "gene")

ggplot(AG_DN_candidates_final, aes(x=gene, y=value))+
  geom_boxplot()+
  theme_classic()+
  xlab(" ")+
  ylab("TPM")+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

### ORF binning
#write.table(RNA_counts_dn, file='MOD_withstats',sep='\t',
#            row.names = T, quote=F)
#write.table(RNA_counts_unann, file='MSD_withstats',sep='\t',
#            row.names = T, quote=F)

ORF_df = read.delim(file = "tables/MU_MSU_ORF_combined.csv", sep=',')
ORF_df$TU = str_split_fixed(ORF_df$ID, '_', 2) [,1]
dn_unfiltered_gtf1=dn_unfiltered_gtf[dn_unfiltered_gtf$V3=='transcript',]
dn_unfiltered_gtf1$TU = str_split_fixed(dn_unfiltered_gtf1$V9, ';', 10) [,9]
dn_unfiltered_gtf1$gene = str_split_fixed(dn_unfiltered_gtf1$V9, ';', 9) [,6]

dn_unfiltered_gtf1$TU = gsub("\"", "", dn_unfiltered_gtf1$TU)
dn_unfiltered_gtf1$gene = gsub("\"", "", dn_unfiltered_gtf1$gene)
dn_unfiltered_gtf1$TU = gsub(" transcript_id ", "", dn_unfiltered_gtf1$TU)
dn_unfiltered_gtf1$gene = gsub(" gene_id ", "", dn_unfiltered_gtf1$gene)
ORF_df$TU[1]

ORF_df$GID=dn_unfiltered_gtf1[match(ORF_df$TU, dn_unfiltered_gtf1$TU),11]
ORF_df=ORF_df[ORF_df$ORF>119,]
RNA_counts_dn$ORF="N"

RNA_counts_dn[row.names(RNA_counts_dn)%in%ORF_df$GID,"ORF"]="Y"

length(which(RNA_counts_dn$ORF=='Y'))
length(which(RNA_counts_dn$ORF=='N'))
t.test(RNA_counts_dn[RNA_counts_dn$ORF=='Y','mean'],RNA_counts_dn[RNA_counts_dn$ORF=='N','mean'])
t.test(RNA_counts_dn[RNA_counts_dn$ORF=='Y','max'],RNA_counts_dn[RNA_counts_dn$ORF=='N','max'])
t.test(RNA_counts_dn[RNA_counts_dn$ORF=='Y','exon'],RNA_counts_dn[RNA_counts_dn$ORF=='N','exon'])
t.test(RNA_counts_dn[RNA_counts_dn$ORF=='Y','num_greater1'],RNA_counts_dn[RNA_counts_dn$ORF=='N','num_greater1'])
t.test(RNA_counts_dn[RNA_counts_dn$ORF=='Y','exonlen'],RNA_counts_dn[RNA_counts_dn$ORF=='N','exonlen'])

colnames(RNA_counts_dn)
ggplot(RNA_counts_dn, aes(x=ORF, y=mean))+
  geom_boxplot()

#shortening ORF df for D sim analysis
MU_ORF_test <- ORF_df[ORF_df$class=="MU",]
MU_ORF_test1<-MU_ORF_test[!duplicated(MU_ORF_test$GID),]
MU_ORF_test1$chrom<-NA
MU_ORF_test1$start<-NA
MU_ORF_test1$end<-NA
MU_ORF_test1$strand<-NA
MU_ORF_test1$coord<-NA

for(i in 1:nrow(MU_ORF_test1)){
  
  GID_i <- MU_ORF_test1[i,'GID']
  MU_ORF_test1$start[i] <- gene_TSS_df_MOD[gene_TSS_df_MOD$geneid==GID_i, 'start1']
  MU_ORF_test1$end[i] <- gene_TSS_df_MOD[gene_TSS_df_MOD$geneid==GID_i, 'end1']
  MU_ORF_test1$chrom[i] <- gene_TSS_df_MOD[gene_TSS_df_MOD$geneid==GID_i, 'chr']
  MU_ORF_test1$strand[i] <- gene_TSS_df_MOD[gene_TSS_df_MOD$geneid==GID_i, 'strand']
}

MU_ORF_test1$coord <- paste0(MU_ORF_test1$chrom, ":", MU_ORF_test1$start, "-", MU_ORF_test1$end)
write.table(MU_ORF_test1, file = "MU_ORF_coords.tsv", sep='\t', quote=F)

#shortening MSU ORF df for D sim analysis
MSU_ORF_test <- ORF_df[ORF_df$class=="MSU",]
MSU_ORF_test1<-MSU_ORF_test[!duplicated(MSU_ORF_test$GID),]
MSU_ORF_test1$chrom<-NA
MSU_ORF_test1$start<-NA
MSU_ORF_test1$end<-NA
MSU_ORF_test1$strand<-NA
MSU_ORF_test1$coord<-NA

for(i in 1:nrow(MSU_ORF_test1)){
  print(i)
  GID_i <- MSU_ORF_test1[i,'GID']
  MSU_ORF_test1$start[i] <- gene_TSS_df_MSD[gene_TSS_df_MSD$geneid==GID_i, 'start1']
  MSU_ORF_test1$end[i] <- gene_TSS_df_MSD[gene_TSS_df_MSD$geneid==GID_i, 'end1']
  MSU_ORF_test1$chrom[i] <- gene_TSS_df_MSD[gene_TSS_df_MSD$geneid==GID_i, 'chr']
  MSU_ORF_test1$strand[i] <- gene_TSS_df_MSD[gene_TSS_df_MSD$geneid==GID_i, 'strand']
}

MSU_ORF_test1$coord <- paste0(MSU_ORF_test1$chrom, ":", MSU_ORF_test1$start, "-", MSU_ORF_test1$end)
write.csv(MSU_ORF_test1, file = "MSU_ORF_coords.csv", quote=F)

#which has longer ORFs?
ORF_consv_df <- read.delim(file='/Users/loganblair/Desktop/seq/github/ORFs-MOD-MSD.csv', na.strings = "", check.names = FALSE, header = T, sep=",")
colnames(ORF_consv_df)[8] <- "ORFlength"
ggplot(ORF_consv_df, aes(y=ORFlength, fill=class, x=conserved))+
  geom_boxplot()+
  scale_fill_manual(values = c("goldenrod", "#C75127"))+
  theme_classic()+
  ylab("ORF length (bp)")+
  xlab("conserved in D. sim")+
  labs(fill="class (ns)")+
  theme(text = element_text(size=20))+
  scale_y_continuous(limits = c(110,850))+
  ggsignif::geom_signif(annotations = c(formatC("*")), 
              y_position = c(770), xmin=c(1), xmax=c(2), textsize = 8)


dataframe1<- data.frame(ORFlength = ORF_consv_df$ORFlength, conserved= ORF_consv_df$conserved,
                        ORFclass = ORF_consv_df$class)
conORFlen_lm<- lm(ORFlength ~ conserved * ORFclass, data = dataframe1)
summary(aov(conORFlen_lm))
aov(conORFlen_lm)
conORFlen_lm
mylogit <- glm(formula = ORF$conserved ~ ORF$ORFlength*ORF$class, family = "binomial" )
summary(mylogit)
myanova <- anova(ORF$ORFlength, )
library(tidyverse)
dat <- ORF %>%
  dplyr::select(ORFlength, conserved, class)
aov(data = dat, formula = ORFlength ~ conserved * class)
TukeyHSD(aov(data = dat, formula = ORFlength ~ conserved + class))

#conserved vs not fisher test
fisher_conserved <- data.frame(not_conserved = c(42,69), conserved =c(9,20))
row.names(fisher_conserved) <- c("MOD", "MSD")
fisher.test(fisher_conserved)

#figure S4 + 5

RNA_counts_S4=read.delim(file='kallisto_expression.tsv',
                      sep='\t', stringsAsFactors = F, header=T)
RNA_counts_S4 <- RNA_counts_S4[row.names(RNA_counts_S4)%in%final_screen_list$x, ]
colnames(RNA_counts_S4)[30] <- "ed10"
colnames(RNA_counts_S4)[31] <- "iso1"

RNA_counts_S5<-data.frame(line = colnames(RNA_counts_S4), value = apply(RNA_counts_S4, 2, sum))
RNA_counts_S5$DGRP <-T
RNA_counts_S5$DGRP[30] <- F
RNA_counts_S5$DGRP[31] <- F

order(RNA_counts_S5$value, decreasing = T)
RNA_counts_S5 <- RNA_counts_S5[order(RNA_counts_S5$value, decreasing = T),]
RNA_counts_S5$line


RNA_counts_S5$orders<-order(RNA_counts_S5$value, decreasing = T)
ggplot(RNA_counts_S5, aes(x=as.factor(orders), y=value, fill=DGRP))+
         geom_bar(stat='identity', color='black')+
  scale_y_continuous(expand=c(0,0), limits = c(0,160))+
  scale_x_discrete(labels =c(RNA_counts_S5$line), guide = guide_axis(angle = 90))+
  xlab("")+
  ylab("sum de novo gene expression (TPM)")+
  theme(text=element_text(size=20),
        legend.position = 'none')

RNA_counts_S4_A<-data.frame(line = colnames(RNA_counts_S4), value = apply(RNA_counts_S4, 2, function(x) length(which(x>1))))
RNA_counts_S4_A$DGRP <-T
RNA_counts_S4_A$DGRP[30] <- F
RNA_counts_S4_A$DGRP[31] <- F

order(RNA_counts_S4_A$value, decreasing = T)
RNA_counts_S4_A <- RNA_counts_S4_A[order(RNA_counts_S4_A$value, decreasing = T),]

RNA_counts_S4_A$orders<-order(RNA_counts_S4_A$value, decreasing = T)
ggplot(RNA_counts_S4_A, aes(x=as.factor(orders), y=value, fill=DGRP))+
  geom_bar(stat='identity', color='black')+
  scale_y_continuous(expand=c(0,0), limits = c(0,26))+
  scale_x_discrete(labels =c(RNA_counts_S4_A$line), guide = guide_axis(angle = 90))+
  xlab("")+
  ylab("#de novo genes expressed >1 TPM")+
  theme(text=element_text(size=20),
        legend.position = 'none')





# graphs
rarefaction_graph
graph_lines_exp
graph_exon_count
graph_tau_enrichment
exonlen_nlines_MSD_plot
graph_x_counts
graph_TPM_threshold
nearest_neighbor
exonlength_sig_graph
MSD_prop_graph
MOD_prop_graph
graph_v_x_exp

#