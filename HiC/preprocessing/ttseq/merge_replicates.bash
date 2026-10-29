samtools merge -@ 20 auxinINTS11.bam data/ttseq_bams/INTS11_IAA_plus_R1_sorted.bam data/ttseq_bams/INTS11_IAA_plus_R2_sorted.bam data/ttseq_bams/INTS11_IAA_plus_R3_sorted.bam
samtools merge -@ 20 noauxinINTS11.bam data/ttseq_bams/INTS11_IAA_minus_R1_sorted.bam data/ttseq_bams/INTS11_IAA_minus_R2_sorted.bam data/ttseq_bams/INTS11_IAA_minus_R3_sorted.bam
