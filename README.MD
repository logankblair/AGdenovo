This repo contains analysis for the manuscript "New insights into the dynamics of de novo gene origin"

Summary of analysis prior to R scripts:

1)  Map reads (STAR)
2)  Reference-guided transcriptome assembly for each DGRP line in Drosophila melanogaster (stringtie)
3)  Merge assemblies, both between DGRP lines and with the reference annotation (TACO)
4)  Use merged assemblies to recalculate gene expression (kalisto)
5)  Generate outgroup transcriptome assemblies to filter de novo genes from merged assemblies (Trinity). These outgroup transcriptomes are taken from Yang et al. (2018).
6)  Align merged assemblies to outgroup transcriptome assemblies + reference genome sequences (BLAST). Unfortunately, I was unable to include these BLAST outputs due to file size limits.

Order and purpose of R scripts:

1.  Generate table for expression of genes within D. melanogaster, parsing kalisto output. Used to determine whether de novo genes meet the minimum expression cutoff.
2.  Generate table for expression of genes outside D. melanogaster, parsing kalisto output. Used to assess whether genes are unique to D. melanogaster.
3.  Determine transcript presence in D. melanogaster and absence in related species. Parses BLAST information and uses expression tables to determine whether genes belong to MOD/MSD categories. Synteny analysis is incoorperated midway through script.
4.  Calculate tau values for annotated genes
5.  Analysis of bidirectional promoters.
6.  Generate figures for manuscript.
