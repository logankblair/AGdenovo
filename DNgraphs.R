library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggpmisc)
library(scales)
options(scipen=999999)

#load RNA count table from kallisto
RNA_counts=read.delim(file='/Users/loganblair/Desktop/seq/count_tables/lengthscaledcounts',
                      sep='\t', stringsAsFactors = F, header=T)
rownames(RNA_counts)=gsub(' ','', rownames(RNA_counts))

#updated - AG specific list after filtering orphan genes within the AG biased annotated gene list
AG_specific_list=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/AG_list_not_ORFan',
                            sep='\t', header=F)

#de novo and AG specific genes
RNA_counts=read.delim(file='/Users/loganblair/Desktop/seq/count_tables/lengthscaledcounts',
                      sep='\t', stringsAsFactors = F, header=T)
RNA_counts1=RNA_counts
RNA_counts=RNA_counts[,c(1:29)]#no iso1 or ed10

MU_list=read.delim(file='/Users/loganblair/Desktop/seq/final_list',
                                        sep='\t',stringsAsFactors = F, header=T)
MSU_list=read.table(file='/Users/loganblair/Desktop/seq/gene_tab_files/unnan_list_nogenome', sep='\t')

#RNA stats
RNA_counts$mean=apply(RNA_counts[,1:29], 1, mean)
RNA_counts$var=apply(RNA_counts[,1:29], 1, var)
RNA_counts$sd=apply(RNA_counts[,1:29], 1, sd)
RNA_counts$max=apply(RNA_counts[,1:29], 1, function(x) max(x))


#order I want to be graphed
RNA_counts$genetype='1'
RNA_counts$genetype[rownames(RNA_counts)%in%AG_specific_list[,1]]='2'
RNA_counts$genetype[rownames(RNA_counts)%in%MU_list[,1]]='4'
RNA_counts$genetype[rownames(RNA_counts)%in%MSU_list[,1]]='3'

#######proportion of lines expressed#####
RNA_counts$num_greater1=apply(RNA_counts[,1:29], 1, function(x) length(which((x)>=1)))
RNA_counts_anscestral=RNA_counts[grep('fbg', rownames(RNA_counts), ignore.case = T),]
RNA_counts_dn=RNA_counts[rownames(RNA_counts)%in%MU_list[,1],]

lines_exp_graph=ggplot(RNA_counts_dn, aes(x=num_greater1))+
  geom_bar(stat='count', color='black', size=0.25, width=0.75, alpha=0.5)+
  scale_fill_manual(values=c("#56B4E9"))+
  scale_x_continuous(expand = c(0, 0), limits = c(0,31), breaks=seq(0,30,5))+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  ylab('gene count')+
  xlab('#lines expressed >1 TPM')
lines_exp_graph

########################################################################
#######graph1:variance vs expression########################################################################
########################################################################

v_x_exp_plot=ggplot(RNA_counts %>%
                      arrange(genetype), aes(x=mean+1, y=var+1, color=genetype))+
  geom_point(size=1.2, alpha=1)+
  scale_y_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x)) +
  ylab('variance + 1')+
  scale_x_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x)) +
  xlab('mean expression + 1 (TPM)')+
  scale_color_manual(values=c( "grey50","#5DB1DD", "#C75127", "goldenrod"), 
                     labels = c("AN (-AG)",  "AG", "MSU","MU"))+
  #facet_zoom(xlim = c(.01, .5),ylim=c(.01,1), zoom.size = 0.25, split=F,horizontal=F, shrink=T)+
  theme_classic()+
  labs(color="gene class")+
  theme(legend.position = c(0.85,0.3),
        text=element_text(size=25))+
  guides(colour = guide_legend(override.aes = list(size=3.5)))
v_x_exp_plot
options(scipen=9999999999999)

########################################################################
######graph 3: chromosomes########################################################################
########################################################################

#getting location info from gtf file
allgtf=read.delim(file='/Users/loganblair/Desktop/seq/gtf_files/dm6_dn2022_unann_gtf', sep='\t', header=F)

#grab genenames, subset MU 
allgtf_MU=allgtf[grep('assembly',allgtf$V2 ),]
allgtf_MU=allgtf_MU[allgtf_MU$V3=='transcript',]
allgtf_MU$gene_id=sapply(strsplit(gsub('gene_id ', '', allgtf_MU$V9), ';'),`[`, 6)
allgtf_MU$gene_id=gsub(" ", "", allgtf_MU$gene_id)
allgtf_MU$transcript_id=sapply(strsplit(gsub('transcript_id ', '', allgtf_MU$V9), ';'),`[`, 9)

allgtf_MSU=allgtf_MU[allgtf_MU$gene_id%in%MSU_list[,1],]
allgtf_MU=allgtf_MU[allgtf_MU$gene_id%in%MU_list[,1],]

transcript_idx1=data.frame(allgtf_MU$transcript_id,allgtf_MU$gene_id)

allgtf_gene=allgtf[allgtf$V3=='gene',]

allgtf_MU1=allgtf_MU[!duplicated(allgtf_MU$gene_id),]
allgtf_unann1=allgtf_MSU[!duplicated(allgtf_MSU$gene_id),]

allgtf_gene$gene_id=sapply(strsplit(gsub('gene_id ', '', allgtf_gene$V9), ';'),`[`, 1)
allgtf_gene$gene_id=gsub(" ", "", allgtf_gene$gene_id)

allgtf_ag=allgtf_gene[allgtf_gene$gene_id%in%AG_specific_list$V1,]
#
allgtf_MSU=allgtf_MSU[!duplicated(allgtf_MSU$gene),]


chr_df1=as.data.frame(unique(allgtf_MU1$V1))
chr_df2=as.data.frame(unique(allgtf_MU1$V1))
chr_df3=as.data.frame(unique(allgtf_MU1$V1))
chr_df4=as.data.frame(unique(allgtf_MU1$V1))

colnames(chr_df1)[1]='chromosome'
colnames(chr_df2)[1]='chromosome'
colnames(chr_df3)[1]='chromosome'
colnames(chr_df4)[1]='chromosome'

chr_df1$proportion=lapply(1:nrow(chr_df1), function(x) length(which(allgtf_gene$V1==chr_df1$chromosome[x])))
chr_df1$proportion1=as.numeric(chr_df1$proportion)/nrow(allgtf_gene)
chr_df2$proportion=lapply(1:nrow(chr_df2), function(x) length(which(allgtf_MU1$V1==chr_df2$chromosome[x])))
chr_df2$proportion1=as.numeric(chr_df2$proportion)/nrow(allgtf_MU1)
chr_df3$proportion=lapply(1:nrow(chr_df3), function(x) length(which(allgtf_ag$V1==chr_df3$chromosome[x])))
chr_df3$proportion1=as.numeric(chr_df3$proportion)/nrow(allgtf_ag)

chr_df4$proportion=lapply(1:nrow(chr_df4), function(x) length(which(allgtf_MSU$V1==chr_df4$chromosome[x])))
chr_df4$proportion1=as.numeric(chr_df4$proportion)/nrow(allgtf_MSU)

chr_df1$type='1ancestral'
chr_df2$type='4denovo'
chr_df3$type='2AG_specific'
chr_df4$type='3unann'
chr_df=rbind(chr_df2,chr_df1,chr_df3, chr_df4)

chr_df=chr_df[!chr_df$chromosome=='mitochondrion_genome',]

ggplot(chr_df, aes(x=chromosome, y=proportion1, fill=type))+
  geom_bar(position='dodge',stat='identity', aes(fill=type), color='black', size=0.25, width=0.75, alpha=0.8)+
  scale_fill_manual(values=c(  "#999999","#56B4E9", '#C75127',"#E69F00"), 
                     labels = c("AN (-AG)", "AG","MSU", 'MU'))+
  #scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(limits = c(0,0.45), breaks=seq(0,.4,.05),expand = c(0, 0))+
  theme_classic()+
  ylab('gene proportion')+
  xlab("")+
  labs(fill=" ")+
  theme(text=element_text(size=15),
        legend.position = c(0.88,0.85))
  #geom_text(aes(x= chromosome,label=proportion),position=position_dodge(width=0.9), vjust=-1)


#X chromosome only
ggplot(chr_df[chr_df$chromosome=='X',], aes(x=type, y=proportion1, fill=type))+
  geom_bar(position='dodge',stat='identity', aes(fill=type), color='black', size=0.25, width=0.75, alpha=0.85)+
  scale_fill_manual(values=c(  "#999999","#56B4E9", '#C75127',"#E69F00"), 
                    labels = c("AN", "AG","OU", 'YU'))+
  scale_x_discrete(labels=c("AN", "AG","MSU", 'MU'))+
  scale_y_continuous(limits = c(0,0.23), breaks=seq(0,.4,.05),expand = c(0, 0))+
  theme_classic()+
  ylab('proportion on X chromosome')+
  xlab("")+
  labs(fill=" ")+
  theme(text=element_text(size=25),
        legend.position = "none")

########################################################################
####distance to TSS########################################################################
########################################################################
library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)

gr_allgenes = import.gff('/Users/loganblair/Desktop/seq/gtf_files/allgenes_MU_MSU_combined.gtf')


#magic to reduce all genes to one range
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

######loking at the enriched, specific genes here
gene_TSS_df$geneclass='an'
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%MU_list[,1]]="dn"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%MSU_list[,1]]="unanc"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%AG_specific_list[,1]]="anag"

EN_DN_list=read.delim(file ='/Users/loganblair/Desktop/seq/r_scripts/EN_DN_list', sep='\t',header=F)
EN_UA_list=read.delim(file='/Users/loganblair/Desktop/seq/r_scripts/EN_UA_list', sep='\t',header=F)
nonEN_DN_list=read.delim(file ='/Users/loganblair/Desktop/seq/r_scripts/nonEN_DN_list', sep='\t',header=F)
nonEN_UA_list=read.delim(file='/Users/loganblair/Desktop/seq/r_scripts/nonEN_UA_list', sep='\t',header=F)

#########
#########
#########
#########
#gene_TSS_df$geneclass[gene_TSS_df$geneid%in%EN_DN_list[,1]]="ENdn"
#gene_TSS_df$geneclass[gene_TSS_df$geneid%in%nonEN_DN_list[,1]]="nonENdn"
#gene_TSS_df$geneclass[gene_TSS_df$geneid%in%EN_UA_list[,1]]="UAdn"
#gene_TSS_df$geneclass[gene_TSS_df$geneid%in%nonEN_UA_list[,1]]="nonUAdn"
#########
#########
#########
#########
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
gene_TSS_df[gene_TSS_df$geneid%in%MSU_list$MU_list1,'group']="UA"
gene_TSS_df[gene_TSS_df$geneid%in%MU_list[,1],'group']="DN"
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
  scale_x_discrete(labels=c('Annotated', 'MU', 'MSU'))+
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
  scale_x_discrete(labels=c('Annotated', 'MU', 'MSU'))+
  ylab('Nearest AG neighbor (bp)')+
  theme(text =element_text(size=25))
ggplot(gene_TSS_df, aes(x=geneclass, y=min_to_AN+1))+
  geom_boxplot()+
  scale_y_log10(limits=c(1,10000000))



ggplot(gene_TSS_df, aes(x=enrichment, fill=chr))+
  facet_col(~group, scales = 'free')+
  geom_bar(stat="count", position="dodge")+
  theme_minimal()




#################

###test10/10/23
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%MU_list[,1]]="ann"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%MU_list[,1]]="dn"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%MSU_list[,1]]="unanc"
gene_TSS_df$geneclass[gene_TSS_df$geneid%in%MSU_list[,1]]="unanc"


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


ggplot(gene_TSS_df1, aes(fill=geneclass, y=log(as.numeric(count)+1), x=gene_class))+
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
ggplot(gene_TSS_df1_plot[gene_TSS_df1_plot$gene_class%in%c("min_anag","min_ua","min_dn"),], aes(fill=geneclass, y=(count+1), x=gene_class))+
  geom_boxplot(outlier.size = 0.1)+
  #scale_y_log10(limits=c(1,10000000))+
  #facet_grid(cols = vars(direction))+
  labs(fill=" ")+
  scale_x_discrete(labels=c('AG', "MU","MSU"))+
  scale_y_log10()+
  scale_fill_manual(values=c("#7A65A5", "#56B4E9", "#E69F00" ,"#C75127"), #grey or purple? #999999#7A65A5
                    labels=c('AN', 'AG', "MU", "MSU"))+
  xlab(' ')+
  ylab('TSS distance + 1')+
  theme_bw()+
  theme(text=element_text(size=25))



gene_TSS_df_tests=gene_TSS_df1[gene_TSS_df1$direction=='opposite orientation',]
gene_TSS_df_tests=as.data.frame(gene_TSS_df_tests)

anova_test=(aov(log(count+1)~geneclass+gene_class, data = gene_TSS_df_tests ))
TukeyHSD(anova_test)
pairwise.wilcox.test(gene_TSS_df_tests$count, gene_TSS_df_tests$geneclass, p.adjust.method = p.adjust.methods,
                     paired = FALSE)

#library('dunn.test')
gene_TSS_df_tests$count~gene_TSS_df_tests$geneclass+gene_TSS_df_tests$gene_class
kruskal.test(count ~ factor(geneclass) + factor(gene_class), data = gene_TSS_df_tests)
??kruskal.test
kruskal.test()

#filter only to main autosomes?

#######length of de novo genes########################################################################

#from https://www.biostars.org/p/83901/
# First, import the GTF-file that you have also used as input for htseq-count

dn_unfiltered_gtf=read.delim(file='/Users/loganblair/Desktop/seq/gtf_files/gtf_21-9-29/gtf_files/1tpm_29asmbl/refcom_29_1tpm/unfiltered_dn_1tpm.gtf', sep='\t',
           quote = "", header=F)
dmel_anc_gtf=read.delim(file='/Users/loganblair/Desktop/seq/gtf_files/gtf_21-9-29/dmel-all-r6.41.gtf', sep='\t',
                        quote = "", header=F)
all_gtf_2=rbind(dn_unfiltered_gtf,dmel_anc_gtf)
write.table(all_gtf_2,file="/Users/loganblair/Desktop/seq/gtf_files/dm6_dn2022gtf", sep='\t', quote=F, row.names = F, col.names = F)

txdb <- makeTxDbFromGFF("/Users/loganblair/Desktop/seq/gtf_files/dm6_dn2022gtf",format="gtf")
txdb
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
gene_length_df = as.data.frame(sum(width(reduce(exons.list.per.gene))))

colnames(gene_length_df)='length'
gene_length_df$type='1_ancestral'

gene_length_df[row.names(gene_length_df)%in%MU_list[,1],'type']='4_MU'
gene_length_df[row.names(gene_length_df)%in%AG_specific_list[,1],'type']='2_AG'
gene_length_df[row.names(gene_length_df)%in%MSU_list[,1],'type']='3_MSU'


ggplot(gene_length_df, aes(y=length, x=type))+
  geom_violin(aes(fill=type), alpha=0.5,outlier.alpha = 1, width=0.6)+
  scale_fill_manual(values=c(  "#999999","#E69F00","#56B4E9"), 
                    labels = c("anscestral", "AG specific","denovo"))+
  scale_colour_manual()+
  scale_y_log10()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 0,
        axis.title.x = element_blank())+
  ylab('Total exon length (bp)')+
  scale_x_discrete(labels = c("anscestral", "AG specific","denovo"))


ggplot(gene_length_df, aes(y=log(length, base = 10), x=type))+
  geom_boxplot(aes(fill=type),alpha=0.8, outlier.alpha = 1, width=0.6)+
  scale_fill_manual(values=c(  "grey50","#5DB1DD","#C75127","goldenrod"), 
                    labels = c("AN (-AG)", "AN AG","MU", "MSU"))+
  scale_colour_manual()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 0,
        axis.title.x = element_blank())+
  ylab('Total exon length (bp)')+
  scale_x_discrete(labels = c("AN (-AG)", "AG","MSU", "MU"))+
  theme(text=element_text(size=25))+
  geom_signif(annotations = c(formatC("***"),
                              formatC("***"),
                              formatC("***"),
                              formatC("NS"),
                              formatC("***"),
                              formatC("***")), 
              y_position = c(5.7, 5.8, 5.9, 5.4, 5.5, 5.1), xmin=c(1.05, 1, .95, 2, 2, 3), xmax=c(2.05, 3, 4.05, 3, 4, 4), textsize = 6, 
              tip_length = 0.03)+
  scale_y_continuous(limits=c(1,6), breaks=c(1,2,3,4,5), labels=c(10,100,1000,10000,100000))
  
10^3
gene_length_df %>% pairwise_wilcox_test(length~type)

#####
ggplot(gene_length_df, aes(y=length, x=type))+
  geom_violin(aes(fill=type), alpha=0.5,outlier.alpha = 1, width=0.6)+
  scale_fill_manual(values=c( "grey50","#5DB1DD","#C75127","goldenrod"), 
                    labels = c("anscestral", "AG specific","denovo", "unannotated"))+
  scale_colour_manual()+
  scale_y_log10()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = 0,
        axis.title.x = element_blank())+
  ylab('Total exon length (bp)')+
  scale_x_discrete(labels = c("anscestral", "AG specific","denovo", "unannotated"))

#library(ggsignif)

unique(gene_length_df$type)

ggplot(gene_length_df, aes(y=length, x=type, fill=type))+
  geom_boxplot(aes(fill=type), alpha=0.8,outlier.alpha = 1, width=0.6)+
  scale_fill_manual(values=c("#56B4E9", "#999999","#C75127", "#E69F00" ))+
  scale_x_discrete(labels=c("AN(-AG)", "AG", "MSU", "MU"))+
  scale_colour_manual()+
  scale_y_log10(limits=c(20,120000))+
  theme_classic()+
  ylab('Total exon length (bp)')+
  #scale_x_discrete(labels = c("AN", "AG","MSU", "MU"))+
  theme(text = element_text(size=25))
  #geom_signif(y_position = c(5.5, 5.7, 5.5), xmin = c(1.025, 2.025, 3.025), xmax = c(1.975, 2.975, 3.975),
  #  annotation = c("***", "***", "***"), tip_length =0)+
  #geom_signif(y_position = c(6, 6.3), xmin = c(1, 2), xmax = c(3, 4),
  #            annotation = c("NS", "***"), tip_length=0)+
  #geom_signif(y_position = c(6.6), xmin = c(1), xmax = c(4),
  #            annotation = c("***"), tip_length=0)
  
  










#############exon counts########
#graph 6 or something

exonscounts_df=data.frame(table(data.frame(exons.list.per.gene)$group_name))
rownames(exonscounts_df)=exonscounts_df[,1]
exonscounts_df$Var1=NULL
exonscounts_df$exon_count=exonscounts_df$Freq
exonscounts_df$Freq=NULL

exonscounts_df$type1='ancestral'
exonscounts_df[row.names(exonscounts_df)%in%MU_list[,1],'type1']='denovo'
exonscounts_df[row.names(exonscounts_df)%in%AG_specific_list[,1],'type1']='AG'

exon_dn=as.data.frame(table(exonscounts_df[exonscounts_df$type1=='denovo','exon_count']))
exon_AG=as.data.frame(table(exonscounts_df[exonscounts_df$type1=='AG','exon_count']))
exon_anc=as.data.frame(table(exonscounts_df[exonscounts_df$type1=='ancestral','exon_count']))
exon_AG$Freq=exon_AG$Freq/length(which(exonscounts_df$type1=='AG'))
exon_dn$Freq=exon_dn$Freq/length(which(exonscounts_df$type1=='denovo'))
exon_anc$Freq=exon_anc$Freq/length(which(exonscounts_df$type1=='ancestral'))
exon_AG$type1='AG'
exon_dn$type1='dn'
exon_anc$type1='anc'

exonscounts_df1=exonscounts_df %>% count(exon_count, type1)

#making a greater than 15 
AAancestral_16=sum(exonscounts_df1[which(exonscounts_df1$type1=='ancestral'&exonscounts_df1$exon_count>7),'n'])
AG_16=sum(exonscounts_df1[which(exonscounts_df1$type1=='AG'&exonscounts_df1$exon_count>7),'n'])
denovo_16=sum(exonscounts_df1[which(exonscounts_df1$type1=='denovo'&exonscounts_df1$exon_count>7),'n'])

exonscounts_df1=exonscounts_df1[-which(exonscounts_df1$type1=='ancestral'&exonscounts_df1$exon_count>7),]
exonscounts_df1=exonscounts_df1[-which(exonscounts_df1$type1=='AG'&exonscounts_df1$exon_count>7),]
exonscounts_df1=exonscounts_df1[-which(exonscounts_df1$type1=='denovo'&exonscounts_df1$exon_count>7),]


exonscounts_df1=rbind(exonscounts_df1, c(8,'ancestral',AAancestral_16))
exonscounts_df1=rbind(exonscounts_df1, c(8,'AG',AG_16))
exonscounts_df1=rbind(exonscounts_df1, c(8,'denovo',denovo_16))


exonscounts_df1$n=as.numeric(exonscounts_df1$n)
exonscounts_df1$exon_count=as.numeric(exonscounts_df1$exon_count)

exonscounts_df1[which(exonscounts_df1$type1=='AG'),'n']=
  exonscounts_df1[which(exonscounts_df1$type1=='AG'),'n']/sum(exonscounts_df1[which(exonscounts_df1$type1=='AG'),'n'])
exonscounts_df1[which(exonscounts_df1$type1=='ancestral'),'n']=
  exonscounts_df1[which(exonscounts_df1$type1=='ancestral'),'n']/sum(exonscounts_df1[which(exonscounts_df1$type1=='AAancestral'),'n'])
exonscounts_df1[which(exonscounts_df1$type1=='denovo'),'n']=
  exonscounts_df1[which(exonscounts_df1$type1=='denovo'),'n']/sum(exonscounts_df1[which(exonscounts_df1$type1=='denovo'),'n'])


             
ggplot(exonscounts_df1, aes(y=n, x=exon_count, fill=type1))+
  geom_bar(position = position_dodge(preserve = "single"),stat='identity', aes(fill=type1), color='black', size=0.25, width=0.75, alpha=0.5)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
          legend.key = element_rect(fill = "white"))+
  xlab('Number of Exons')+
  ylab("proportion")+
  scale_x_continuous(limits = c(0.1,8.5), breaks=seq(0,8,1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0,0.62), breaks=seq(0,.8,0.1), expand = c(0, 0))+
  scale_fill_manual(values=c("#999999","#E69F00", "#56B4E9"), 
                    labels = c("anscestral (n=17570)", "AG specific (n=537)", "denovo (n=122)"))
  










#############exon counts####################
exonscounts_df=data.frame(table(data.frame(exons.list.per.gene)$group_name))
rownames(exonscounts_df)=exonscounts_df[,1]
exonscounts_df$Var1=NULL
exonscounts_df$exon_count=exonscounts_df$Freq
exonscounts_df$Freq=NULL

exonscounts_df$type1='ancestral'
exonscounts_df[row.names(exonscounts_df)%in%MU_list[,1],'type1']='denovo'
exonscounts_df[row.names(exonscounts_df)%in%MSU_list[,1],'type1']='unannotated'
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
ggplot(exonscounts_df1, aes(y=n, x=exon_count, fill=type1))+
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
  scale_fill_manual(values=c("#999999", "#56B4E9",'#C75127', "#E69F00"), 
                    labels = c("AN", "AG", "MSU", 'MU'))+
  theme(legend.position = c(0.88, 0.83),
        text=element_text(size=25))

median(exonscounts_df[which(exonscounts_df$type1=='ancestral'),1])
median(exonscounts_df[which(exonscounts_df$type1=='AG'),1])
median(exonscounts_df[which(exonscounts_df$type1=='unannotated'),1])
median(exonscounts_df[which(exonscounts_df$type1=='denovo'),1])

length(which(exonscounts_df$type1=='unannotated'))


####
#graph 6 but with unannotated column
#############wilcox rank sum tests####################



##Number per chromosome - used pearson chi square online
#denovo not sig different from either

#############dist to ancestral genes
#
test <- wilcox.test(as.numeric(TSS_graph[TSS_graph$type=='Aancestral','dist_to_TSS']),
                    as.numeric(TSS_graph[TSS_graph$type=='denovo','dist_to_TSS']))
test #de novo and AG both sig farther than Ancest, but not different from each other

#dn vs AG 0.02055
#dn vs anc NS 0.5047
#AG vs anc 0.00002077

test <- wilcox.test(as.numeric(TSS_graph[TSS_graph$type=='Aancestral','dist_to_TSS_AG']),
                    as.numeric(TSS_graph[TSS_graph$type=='denovo','dist_to_TSS_AG']))
#dn vs anc 0.02819
#others obviously very sig




###########EXONLENGTHEXONLENGTHEXONLENGTHEXONLENGTHEXONLENGTH#############
test <- wilcox.test(as.numeric(gene_length_df[gene_length_df$type=='AG','length']),
                    as.numeric(gene_length_df[gene_length_df$type=='denovo','length']))
test #de novo and AG both sig longer than Ancest, but not different from each other

test=pairwise.wilcox.test(gene_length_df$length, gene_length_df$type)
test


#########EXONS##########
#####for exons all pairwise are really sig

pairwise.wilcox.test(exonscounts_df$exon_count, exonscounts_df$type1)

test[1,2]=wilcox.test(as.numeric((exonscounts_df[exonscounts_df$type=='AAancestral','exon_count'])),
                    as.numeric(log(exonscounts_df[exonscounts_df$type=='AG','exon_count'])))$p.value
test[1,3]=wilcox.test(as.numeric(log(exonscounts_df[exonscounts_df$type=='AAancestral','exon_count'])),
                    as.numeric(log(exonscounts_df[exonscounts_df$type=='denovo','exon_count'])))$p.value
test[1,4]=wilcox.test(as.numeric(log(exonscounts_df[exonscounts_df$type=='AAancestral','exon_count'])),
                    as.numeric(log(exonscounts_df[exonscounts_df$type=='unannotated','exon_count'])))$p.value
test[2,3]=wilcox.test(as.numeric(log(exonscounts_df[exonscounts_df$type=='AG','exon_count'])),
                    as.numeric(log(exonscounts_df[exonscounts_df$type=='denovo','exon_count'])))
test[2,4]=wilcox.test(as.numeric(log(exonscounts_df[exonscounts_df$type=='AG','exon_count'])),
                      as.numeric(log(exonscounts_df[exonscounts_df$type=='denovo','exon_count'])))
test <- wilcox.test(as.numeric(log(exonscounts_df[exonscounts_df$type=='AG','exon_count'])),
                    as.numeric(log(exonscounts_df[exonscounts_df$type=='denovo','exon_count'])))
test

#A-AG 0.00000001791
#A-dn 0.000000008985
#AG-dn 0.0001076


##########Expression graphs##########


label_fill <- function(orig, .offset=1, .mod=2, .fill=""){
  ## replace
  ii <- as.logical(
    ## offset==0 keeps first
    (1:length(orig)-1+.offset) %% .mod
  )
  orig[ii] <- .fill
  orig
}

every_third <- scale_x_continuous(
  breaks=my_breaks, 
  labels=label_fill(my_breaks, .mod=1)
)

ggplot(RNA_counts_dn, aes(x=num_greater2))+
  geom_bar(stat='count', fill="#56B4E9", alpha=0.5, color='black', size=0.15)+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_rect(fill = "white"))+
  scale_x_continuous(
    breaks=1:28, labels=label_fill, expand = c(0.02, 0)
  )+
  scale_y_continuous(limits = c(0,12), breaks=seq(0,12,1), expand = c(0, 0))+
  xlab("# lines with expression > 2 cpm")+
  ylab("# de novo genes")

#

#library('diptest')

RNA_counts[,1:29]
RNA_counts$dip_pvalues=apply(RNA_counts[,1:29], 1, function(x) dip.test(x)$p.value)
RNA_counts_dn$diptest=apply(RNA_counts_dn[,3:9], 1, function(x) dip.test(x)$p.value)
FBGNid
which(RNA_counts[rownames(RNA_counts)%in%MU_list[,1],]$dip_pvalues<0.05)
sort(RNA_counts_dn$diptest)


#distribution of de novo genes expressed
RNA_counts_dn$max
ggplot(RNA_counts_dn, aes(x=max))+
  geom_histogram(bins=60)+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  xlab("max TPM")

ggplot(RNA_counts_dn, aes(x=max, y=num_greater1))+
  geom_point()+
  theme_classic()+
  ylab("# lines expressed greater than 1 TPM")+
  xlab("max TPM expression")+
  scale_x_log10()


###################Simulans rarity and unannotated rarity###################

unann_list=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/outgroup_unannotated', sep='\t', header = F)

unann_list_nogenome=read.delim(file='/Users/loganblair/Desktop/seq/final_UA_list', sep='\t', header = T)

simfixed_genelist=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/genelist_simfixed', sep='\t', header = F)
sim1genelist=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/gene_list_sim1', sep='\t', header = F)
sim2genelist=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/gene_list_sim2', sep='\t', header = F)
sim3genelist=read.delim(file='/Users/loganblair/Desktop/seq/gene_tab_files/gene_list_sim3', sep='\t', header = F)

RNA_counts_unann=RNA_counts[rownames(RNA_counts)%in%MSU_list[,1],]

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
RNA_counts_unann[1,]
RNA_counts_dn[1,]

RNA_counts_dn$num_greater1=apply(RNA_counts_dn[,c(1:29)], 1, function(x) length(which(x>1)))
RNA_counts_unann$num_greater1=apply(RNA_counts_unann[,c(1:29)], 1, function(x) length(which(x>1)))


unan_vs_dn_df=rbind(RNA_counts_unann[,c(1:36)], RNA_counts_dn[,c(1:36)])
unan_vs_dn_df$num_greater1=apply(unan_vs_dn_df[,c(1:29)], 1, function(x) length(which(x>1)))
unan_vs_dn_df[1,]
unan_vs_dn_df$num_greater1
ggplot(unan_vs_dn_df, aes(x=as.numeric(num_greater1), fill=genetype))+
  geom_histogram(color= 'black',  
                 position='stack', alpha=1,
                 #position='identity', alpha=0.75,
                 bins=30)+
  theme_classic()+
  guides(fill=guide_legend(" "))+
  scale_fill_manual(values=c("goldenrod", 'grey40'), labels=c('Un; DN', 'Un; Anc'))+
  scale_y_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5), limits=c(0,100))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  xlab("#lines expressed >1 TPM")+
  ylab('count')+
  theme(legend.position=c(.8,.75))+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               y = 40, yend = 35,arrow = arrow(length = unit(0.03, "npc")), color="goldenrod", size=1)+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               y = 40, yend = 35,arrow = arrow(length = unit(0.03, "npc")), color="grey40", size=1)
  
MSUprop= mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1'])
MSUprop/29
MUprop=mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1'])
MUprop/29

#de novo 
ggplot(unan_vs_dn_df[unan_vs_dn_df$genetype=='4',], aes(x=as.numeric(num_greater1), fill=genetype))+
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
  ylab('MU count')+
  theme(legend.position='none',
        text = element_text(size=25))+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               y = 20, yend = 15,arrow = arrow(length = unit(0.04, "npc")), color="goldenrod", size=1.5)
#MSU
ggplot(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated',], aes(x=as.numeric(num_greater1), fill=genetype))+
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
  ylab('MSU count')+
  theme(legend.position='none',
        text = element_text(size=25))+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               y = 20, yend = 15,arrow = arrow(length = unit(0.04, "npc")), color="#C75127", size=1.5)



wilcox.test(unan_vs_dn_df[unan_vs_dn_df$genetype=="Unannotated","num_greater1"],
            unan_vs_dn_df[unan_vs_dn_df$genetype=="4","num_greater1"])
            






unan_vs_dn_df_counts=unan_vs_dn_df%>% count(genetype,num_greater1) %>%
  tidyr::complete(genetype,num_greater1)

ggplot(unan_vs_dn_df_counts, aes(y=n, x=num_greater1, fill=genetype))+
  geom_bar(color= 'black',stat='identity',
                 position='dodge', alpha=0.8, width=0.7)+
  theme_classic()+
  guides(fill=guide_legend("unannotated type"))+
  scale_fill_manual(values=c("goldenrod", 'grey40'), labels=c('mel only (n=122)', 'mel + yak/sim (n=147)'))+
  scale_y_continuous(expand=c(0,0), limits = c(0,55))+
  xlab("#lines expressed > 1tpm")+
  ylab("count")+
  theme(legend.position=c(.8,.75))+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']),
               y = 40, yend = 35,arrow = arrow(length = unit(0.03, "npc")), color="goldenrod", size=1)+
  geom_segment(x = mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               xend =  mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1']),
               y = 40, yend = 35,arrow = arrow(length = unit(0.03, "npc")), color="grey40", size=1)

4.806723/29
mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1'])
length(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']==1)
length(which(unan_vs_dn_df[unan_vs_dn_df$genetype=='4','num_greater1']==1))
mean(unan_vs_dn_df[unan_vs_dn_df$genetype=='Unannotated','num_greater1'])
     
##

#############rarefaction#############
RNA_counts_dn1=RNA_counts_dn[,c(1:29)]

temp_df=RNA_counts_dn1
rarefaction=as.data.frame(c(1:29))
rarefaction$count=NA
for(i in 1:28){
  tmp_col=which.max(apply(temp_df, 2, function(x) length(which(x>1))))
  tmp_dn_genes=which(temp_df[,tmp_col]>1)
  
  rarefaction$count[i]=length(tmp_dn_genes)
  
  temp_df=temp_df[-tmp_dn_genes,]
  temp_df=temp_df[,-tmp_col]
  print(i)
}

RNA_counts_DNunann=rbind(RNA_counts[rownames(RNA_counts)%in%MSU_list[,1],],
                         RNA_counts[rownames(RNA_counts)%in%MU_list$x,])

RNA_counts_DNunann=RNA_counts_DNunann[,c(1:29)]

temp_df=RNA_counts_DNunann
rarefaction$unann_count=NA
for(i in 1:28){
  tmp_col=which.max(apply(temp_df, 2, function(x) length(which(x>1))))
  tmp_dn_genes=which(temp_df[,tmp_col]>1)
  
  rarefaction$unann_count[i]=length(tmp_dn_genes)
  
  temp_df=temp_df[-tmp_dn_genes,]
  temp_df=temp_df[,-tmp_col]
  print(i)
}

rarefaction$count[29]=0
rarefaction$unann_count[29]=0

rarefaction1=as.data.frame(c(cumsum(rarefaction$count),cumsum(rarefaction$unann_count)))
names(rarefaction1)='counts'

rarefaction1$type=c(rep('denovo', times=29), rep('1unann', times=29))
rarefaction1$num=c(rep(c(1:29), times=2))

ggplot(rarefaction1, aes(x=num, y=counts, color=type)) + 
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
                     label=c('MU+MSU','MU'))+
  theme(text = element_text(size=25),
        legend.title = element_blank())

#############common vs rare#############

RNA_counts_dn$exon=exonscounts_df[(row.names(exonscounts_df)%in%row.names(RNA_counts_dn)),1]
RNA_counts_unann$exon=exonscounts_df[exonscounts_df$type1=='unannotated',1]

my.formula <- y ~ x #for linreg lines

RNA_counts_dn%>%add_count(exon,num_greater1)%>%ggplot(aes(y=exon, x=num_greater1))+
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

RNA_counts_unann%>%add_count(exon,num_greater1)%>%ggplot(aes(y=exon, x=num_greater1))+
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

cor.test(RNA_counts_unann$exon, RNA_counts_unann$num_greater1)
cor.test(RNA_counts_dn$exon, RNA_counts_dn$num_greater1)


ggplot(RNA_counts_dn, aes(y=exonlen, x=num_greater1))+
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


ggplot(RNA_counts_unann, aes(y=exonlen, x=num_greater1))+
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

cor.test(RNA_counts_dn$num_greater1, RNA_counts_dn$exonlen)
cor.test(RNA_counts_unann$num_greater1, RNA_counts_unann$exonlen)

#
nrow(RNA_counts_dn)
ggplot(RNA_counts_unann, aes(y=max, x=num_greater1))+
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
b=a=lm(RNA_counts_dn$max~RNA_counts_dn$num_greater1)
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


#DN info table
info_table_MU=RNA_counts_dn
info_table_MU$genetype=NULL
info_table_MU$type=NULL
info_table_MU$exonlen=NULL
info_table_MU$exon=NULL

info_table_MU$transcript_length=gene_length_df[rownames(gene_length_df)%in%rownames(info_table_MU),'length']
info_table_MU$exon_number=exonscounts_df[rownames(exonscounts_df)%in%rownames(info_table_MU),'exon_count']
info_table_MU$start=gene_TSS_df[gene_TSS_df$geneid%in%rownames(info_table_MU),'start1']
info_table_MU$end=gene_TSS_df[gene_TSS_df$geneid%in%rownames(info_table_MU),'end1']
info_table_MU$chr=gene_TSS_df[gene_TSS_df$geneid%in%rownames(info_table_MU),'chr']
info_table_MU$strand=gene_TSS_df[gene_TSS_df$geneid%in%rownames(info_table_MU),'strand']
info_table_MU$ed10=RNA_counts1[rownames(RNA_counts1)%in%rownames(info_table_MU),30]
info_table_MU$iso1=RNA_counts1[rownames(RNA_counts1)%in%rownames(info_table_MU),31]
info_table_MU$gene_id=rownames(info_table_MU)
ncol(info_table_MU)
outgroup_exp=read.delim("/Users/loganblair/Desktop/seq/Rplots/paper1_v2/unannotated_outgroup_expression.tsv", sep='\t',
                        header=T)
info_table_MU$w501=outgroup_exp[rownames(outgroup_exp)%in%info_table_MU$gene_id,'w501']
info_table_MU$lara10=outgroup_exp[rownames(outgroup_exp)%in%info_table_MU$gene_id,'lara10']
info_table_MU$SZ116=outgroup_exp[rownames(outgroup_exp)%in%info_table_MU$gene_id,'X116']
info_table_MU$tai18=outgroup_exp[rownames(outgroup_exp)%in%info_table_MU$gene_id,'z_tai18']

info_table_MU=info_table_MU[,c(43,39,37,38,40,35,36,1:29,41:47,30:34)]
write.table(info_table_MU, file='/Users/loganblair/Desktop/seq/AG_denovo_R/info_table_MU_genes.tsv', sep='\t',
            row.names=F, col.names = T, quote=F)

#MSU
info_table_MSU=RNA_counts_unann

info_table_MSU$genetype=NULL
info_table_MSU$type=NULL
info_table_MSU$exonlen=NULL
info_table_MSU$exon=NULL

info_table_MSU$transcript_length=gene_length_df[rownames(gene_length_df)%in%rownames(info_table_MSU),'length']
info_table_MSU$exon_number=exonscounts_df[rownames(exonscounts_df)%in%rownames(info_table_MSU),'exon_count']

info_table_MSU$start=allgtf_MSU[allgtf_MSU$gene_id%in%rownames(info_table_MSU),4]
info_table_MSU$end=allgtf_MSU[allgtf_MSU$gene_id%in%rownames(info_table_MSU),5]
info_table_MSU$chr=allgtf_MSU[allgtf_MSU$gene_id%in%rownames(info_table_MSU),1]
info_table_MSU$strand=allgtf_MSU[allgtf_MSU$gene_id%in%rownames(info_table_MSU),7]
info_table_MSU$ed10=RNA_counts1[rownames(RNA_counts1)%in%rownames(info_table_MSU),30]
info_table_MSU$iso1=RNA_counts1[rownames(RNA_counts1)%in%rownames(info_table_MSU),31]
info_table_MSU$gene_id=rownames(info_table_MSU)

outgroup_exp=read.delim("/Users/loganblair/Desktop/seq/Rplots/paper1_v2/unannotated_outgroup_expression.tsv", sep='\t',
           header=T)
info_table_MSU$w501=outgroup_exp[rownames(outgroup_exp)%in%info_table_MSU$gene_id,'w501']
info_table_MSU$lara10=outgroup_exp[rownames(outgroup_exp)%in%info_table_MSU$gene_id,'lara10']
info_table_MSU$SZ116=outgroup_exp[rownames(outgroup_exp)%in%info_table_MSU$gene_id,'X116']
info_table_MSU$tai18=outgroup_exp[rownames(outgroup_exp)%in%info_table_MSU$gene_id,'z_tai18']

info_table_MSU=info_table_MSU[,c(43,39,37,38,40,35,36,1:29,41:47,30:34)]
info_table_MSU$gene_id.1=NULL
write.table(info_table_MSU, file='/Users/loganblair/Desktop/seq/AG_denovo_R//info_table_MSU_genes.tsv', sep='\t',
            row.names=F, col.names = T, quote=F)


