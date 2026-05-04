
# try su adeno1 -> Default, nessun parametro toccato
java -jar -Xmx40g /home/user/nextflow/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar scanfastq -d /media/user/8Tb/adeno1_run1_fastq_pass/ -o adeno1_sicelore --bcEditDistance 1
BUILD=~/reference/refdata-gex-GRCm39-2024-A/fasta/genome.fa
minimap2 -ax splice -uf --sam-hit-only -t 30 $BUILD passed/*.fastq.gz | samtools view -bS -@ 20 - | samtools sort -m 2G -@ 20 -o passed.bam -&& samtools index passed.bam
