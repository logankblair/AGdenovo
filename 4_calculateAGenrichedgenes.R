setwd("~/Desktop/github1/")

options(scipen=999999)

#load RNA count table from kalisto
RNA_counts=read.delim(file='kallisto_expression.tsv',
                      sep='\t', stringsAsFactors = F, header=T)

#####tau vs exp######
#load table from Leader et al 2018, showing expression of genes across different tissues in D mel
tau_df=read.delim(file="tables/FlyAtlas_data_table.tsv", sep='\t', header=T)
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

#enrichment graph 
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
write.table(AG_specific_list, file='tables/AG_biased_gene_list.tsv',
            sep='\t',col.names=F,row.names=F)
