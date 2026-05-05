


# Sicelore pipeline for single cell nanopore data
# try su adeno1 -> Default, nessun parametro toccato
java -jar -Xmx40g /home/user/nextflow/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar scanfastq -d /media/user/8Tb/adeno1_run1_fastq_pass/ -o adeno1_sicelore --bcEditDistance 1
BUILD=~/reference/refdata-gex-GRCm39-2024-A/fasta/genome.fa
minimap2 -ax splice -uf --sam-hit-only -t 30 $BUILD passed/*.fastq.gz | samtools view -bS -@ 20 - | samtools sort -m 2G -@ 20 -o passed.bam -&& samtools index passed.bam

# Step 2
java -jar -Xmx40g /home/user/nextflow/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar assignumis --inFileNanopore passed.bam --outfile adeno1_sicelore_assigned.bam

# Step 3
java -jar -Xmx32g /home/user/nextflow/sicelore-2.1/Jar/Sicelore-2.1.jar SelectValidCellBarcode I=BarcodesAssigned.tsv O=ValidBarcodes.csv MINUMI=1 ED0ED1RATIO=1

## Formatting gtf to refFlat
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz
java /home/user/nextflow/sicelore-2.1/Jar/Sicelore-2.1.jar gtfToGenePred -genePredExt -geneNameAsName2 gencode.vM38.primary_assembly.annotation.gtf gencode.vM38.primary_assembly.annotation.refflat.txt
paste <(cut -f 12 gencode.vM38.primary_assembly.annotation.refflat.txt) <(cut -f 1-10 gencode.vM38.primary_assembly.annotation.refflat.txt) > gencode.vM38.refFlat

# Step 4a
java -jar -Xmx64g /home/user/nextflow/sicelore-2.1/Jar/Sicelore-2.1.jar IsoformMatrix I=passedParsed.bam GENETAG=GE UMITAG=U8 CELLTAG=BC REFFLAT=gencode.vM38.refFlat CSV=barcodes.csv DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false OUTDIR=. PREFIX=sicelore