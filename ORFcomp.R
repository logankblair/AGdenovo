
MSU_ORF_list=read.delim(file="/Users/loganblair/Desktop/seq/ORF/Nov2023/MSU_cpat.longest.orf", header=F)
#remove even rows
MSU_ORF_list=MSU_ORF_list[-2*(1:(nrow(MSU_ORF_list)/2)),]
MSU_ORF_list=as.data.frame(MSU_ORF_list)
#reformat columns
MSU_ORF_list[,1]=gsub(">", "", MSU_ORF_list[,1])
MSU_ORF_list$transcript_id=str_split_fixed(MSU_ORF_list$MSU_ORF_list,pattern = '_', n=6)[,1]
MSU_ORF_list$length=str_split_fixed(MSU_ORF_list$MSU_ORF_list,pattern = '_', n=6)[,3]
MSU_ORF_list$length=str_split_fixed(MSU_ORF_list$length,pattern = ':', n=6)[,2]

#convert transcript to gene ID
big_gtf=read.delim(file="/Users/loganblair/Desktop/seq/gtf_files/gtf_21-9-29/gtf_files/1tpm_29asmbl/refcom_29_1tpm/unfiltered_dn_1tpm.gtf",
           header=F, sep='\t')
big_gtf=big_gtf[big_gtf$V3=='transcript',]
big_gtf$gene_id=str_split_fixed(big_gtf[,9], pattern = ';', n=9)[,6]
big_gtf$gene_id=gsub(" gene_id ", "", big_gtf$gene_id)
big_gtf$transcript_id=str_split_fixed(big_gtf[,9], pattern = ';', n=10)[,9]
big_gtf$transcript_id=gsub(" transcript_id ", "", big_gtf$transcript_id)

big_gtf[which(MSU_ORF_list$transcript_id[1]==big_gtf$transcript_id),'gene_id']
big_gtf[MSU_ORF_list$transcript_id[1]%in%big_gtf$transcript_id,'gene_id']
MSU_ORF_list$gene_id=NA
for(i in 1:nrow(MSU_ORF_list)){
  MSU_ORF_list$gene_id[i]=big_gtf[which(MSU_ORF_list$transcript_id[i]==big_gtf$transcript_id),'gene_id']
}
big_gtf[big_gtf$transcript_id%in%MSU_ORF_list$transcript_id, 'gene_id']
length(unique(MSU_ORF_list$gene_id))



MU_ORF_list=read.delim(file="/Users/loganblair/Desktop/seq/ORF/Nov2023/MU_cpat.longest.orf", header=F)
#remove even rows
MU_ORF_list=MU_ORF_list[-2*(1:(nrow(MU_ORF_list)/2)),]
MU_ORF_list=as.data.frame(MU_ORF_list)
#reformat columns
MU_ORF_list[,1]=gsub(">", "", MU_ORF_list[,1])
MU_ORF_list$transcript_id=str_split_fixed(MU_ORF_list$MU_ORF_list,pattern = '_', n=6)[,1]
MU_ORF_list$length=str_split_fixed(MU_ORF_list$MU_ORF_list,pattern = '_', n=6)[,3]
MU_ORF_list$length=str_split_fixed(MU_ORF_list$length,pattern = ':', n=6)[,2]

#convert transcript to gene ID
big_gtf=read.delim(file="/Users/loganblair/Desktop/seq/gtf_files/gtf_21-9-29/gtf_files/1tpm_29asmbl/refcom_29_1tpm/unfiltered_dn_1tpm.gtf",
                   header=F, sep='\t')
big_gtf=big_gtf[big_gtf$V3=='transcript',]
big_gtf$gene_id=str_split_fixed(big_gtf[,9], pattern = ';', n=9)[,6]
big_gtf$gene_id=gsub(" gene_id ", "", big_gtf$gene_id)
big_gtf$transcript_id=str_split_fixed(big_gtf[,9], pattern = ';', n=10)[,9]
big_gtf$transcript_id=gsub(" transcript_id ", "", big_gtf$transcript_id)

big_gtf[which(MU_ORF_list$transcript_id[1]==big_gtf$transcript_id),'gene_id']
big_gtf[MU_ORF_list$transcript_id[1]%in%big_gtf$transcript_id,'gene_id']
MU_ORF_list$gene_id=NA
for(i in 1:nrow(MU_ORF_list)){
  MU_ORF_list$gene_id[i]=big_gtf[which(MU_ORF_list$transcript_id[i]==big_gtf$transcript_id),'gene_id']
}
big_gtf[big_gtf$transcript_id%in%MU_ORF_list$transcript_id, 'gene_id']
length(unique(MU_ORF_list$gene_id))

MU_ORF_list$length=as.numeric(MU_ORF_list$length)
MU_ORF_list=MU_ORF_list[order(MU_ORF_list$length, decreasing = T),]
MU_ORF_list_longest=MU_ORF_list[!duplicated(MU_ORF_list$gene_id),]

MSU_ORF_list$length=as.numeric(MSU_ORF_list$length)
MSU_ORF_list=MSU_ORF_list[order(MSU_ORF_list$length, decreasing = T),]
MSU_ORF_list_longest=MSU_ORF_list[!duplicated(MSU_ORF_list$gene_id),]


MSU_ORF_list_longest$class="MSU"
MU_ORF_list_longest$class="MU"
MSU_ORF_list_longest$MSU_ORF_list=NULL
MU_ORF_list_longest$MU_ORF_list=NULL
ORF_df=rbind(MSU_ORF_list_longest, MU_ORF_list_longest)

a=ggplot(ORF_df, aes(y=length, x=class)+
  geom_boxplot(size=2)
a
t.test(ORF_df[ORF_df$class=="MU",'length'], ORF_df[ORF_df$class=="MSU",'length'])
