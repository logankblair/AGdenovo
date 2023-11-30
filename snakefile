#input is list of fastq RNA-seq libraries, output is line-specific bamfile and the stringtie transcriptome assemblies
#trimming

rule trim_galore_pe:
    input:
        ["{input}_1.fastq.gz", "{input}_2.fastq.gz"]
    output:
        "trimmed/{input}_1_val_1.fq.gz", "trimmed/{input}_2_val_2.fq.gz"
    shell:
        "trim_galore {input} --paired --retain_unpaired -o trimmed"

#mapping
rule star_map:
    input:
        fq1 = "trimmed/{sample}_1_val_1.fq.gz", fq2= "trimmed/{sample}_2_val_2.fq.gz"
    output:
        "star_{sample}_1_val_1.fq.gz_p.Aligned.sortedByCoord.out.bam"
    shell:
        "star_map.sh -t 2 -1 {input.fq1} -2 {input.fq2}  --idx idx/_dmel-all-chromosome-r6.41 --alignWindowsPerReadNmax 30000  --outSAMstrandField intronMotif  --outSAMattributes All  --outFileNamePrefix ./mapped/{wildcards.sample}"

#stringtie
rule stringtie:
    input:
        "star_{sample}_1_val_1.fq.gz_p.Aligned.sortedByCoord.out.bam"
    output:
        "gtf/{sample}.gtf"
    shell:
        "stringtie {input} -G /genome_files/dmel-all-r6.41.gtf --rf -A genetab/{wildcards.sample}_gene_abund.tab -o gtf/{wildcards.sample}.gtf"