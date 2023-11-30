#Shell one liners used in analysis of "New insights into the dynamics of de novo gene origin from increased population sampling of Drosophila melanogaster accessory glands"
#first three lines are included in the snakefile but are also included them here to create a comprehensive list

#read trimming with trim_galore
"trim_galore {input} --paired --retain_unpaired -o trimmed

#mapping with star. Requires dmel star idx to be buit first
STAR  --genomeDir _dmel-all-chromosome-r6.41 --readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat --outFilterMatchNminOverLread  --chimSegmentMin 29 --quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx --twopassMode Basic --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./mapped/{wildcards.sample} --alignWindowsPerReadNmax 30000

#stringtie
stringtie {input_bam} -G /genome_files/dmel-all-r6.41.gtf --rf -A genetab/{wildcards.sample}_gene_abund.tab -o gtf/{wildcards.sample}.gtf

#merge 29 assemblies with taco
taco_run gtf_list.txt -o 1tpm_assembly.gtf --ref-genome-fasta dmel-all-chromosome-r6.41.fasta --filter-min-expr 1 --gtf-expr-attr TPM

#take merged assemblies (29) and merge with dmel reference, indicating unique transcripts not found in reference
taco_refcomp -o refcomp_final_1tpm -r dmel-all-chromosome-r6.41.gtf -t 1tpm_assembly.gtf

#combine output fasta with dmel reference
cat refcomp_final_1tpm.fasta dmel-all-transcript-r6.41.fasta > combined.fasta

#build kallisto index using fasta of combined unannotated sequences and annotated sequences
kallisto index -i kal_unan_idx -combined.fasta 

#measure expression (again) using kallisto
kallisto quant {sample}1.fq.gz {sample}2.fq.gz -i kal_unan_idx --rf-stranded -o kal_quant_{sample}

#run trinity for 3 D sim and 1 D yak RNA seq replicates for AG tissue
Trinity --seqType fq --left {sample}_1.fq --right {sample}_2.fq

#construct blast database of outgroup sequences
makeblastdb -in {outgroup_fasta} -input_type fasta -dbtype nucl -title {outgroup_fasta}.db -parse_seqids -out {outgroup_fasta}.db

#blast match to outgroup
blastn-db  {outgroup_database}  -query stringtie_1tpm_transcripts.fasta -out {outgroup_database}_.whole_body -outfmt 6 -evalue 1e-6