hisat2 -x ../../data/ref/yeast_index -U ../../data/raw/SRR1916152.fastq | samtools view -bS - | samtools sort -o EV_3.srt.bam
hisat2 -x ../../data/ref/yeast_index -U ../../data/raw/SRR1916153.fastq | samtools view -bS - | samtools sort -o EV_4.srt.bam
hisat2 -x ../../data/ref/yeast_index -U ../../data/raw/SRR1916154.fastq | samtools view -bS - | samtools sort -o DNMT3B_2.srt.bam
hisat2 -x ../../data/ref/yeast_index -U ../../data/raw/SRR1916155.fastq | samtools view -bS - | samtools sort -o DNMT3B_3.srt.bam
hisat2 -x ../../data/ref/yeast_index -U ../../data/raw/SRR1916156.fastq | samtools view -bS - | samtools sort -o DNMT3B_4.srt.bam
