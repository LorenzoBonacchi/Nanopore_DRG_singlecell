# Workstation Adeno1 Run
# here samtools seems to work
# ------------------------------------------ #
# step 1: barcode and UMI assignment
# ------------------------------------------ #
# ------------------------------------------ #

java -jar -Xmx40g ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar scanfastq -d ../adeno1_run1_fastq_pass/ -o /home/lab-user/data/adeno1_sicelore --bcEditDistance 1

# ------------------------------------------ #
# step 2: alignment to minimap2
# ------------------------------------------ #
# ------------------------------------------ #

BUILD=~/reference/refdata-gex-GRCm39-2024-A/fasta/genome.fa
minimap2 -ax splice -uf --sam-hit-only -t 40 $BUILD passed/*passed.fastq.gz | samtools view -bS -@ 20 - | samtools sort -m 2G -@ 20 -o adeno1_passed.bam -&& samtools index adeno1_passed.bam

# ------------------------------------------ #
# step 3_ umi assignment
# ------------------------------------------ #
# ------------------------------------------ #

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz
#java -jar /home/lab-user/tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar \
#gtfToGenePred \
#-genePredExt \
#-geneNameAsName2 \
#../../reference/gencode.vM38.primary_assembly.annotation.gtf \
#gencode.vM38.primary_assembly.annotation.refflat.txt
#
#paste <(cut -f 12 gencode.vM38.primary_assembly.annotation.refflat.txt) <(cut -f 1-10 gencode.vM38.primary_assembly.annotation.refflat.txt) > gencode.vM38.refFlat

## Alla fine ho usato il gtf e non il refFlat, sembra funzionare lo stesso

java -jar -Xmx40g ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar assignumis --inFileNanopore adeno1_passed.bam --outfile adeno1_sicelore_assigned.bam --ONTgene GE -a ../../reference/gencode.vM38.primary_assembly.annotation.gtf

# ------------------------------------------ #
# Step 4a
# ------------------------------------------ #
# ------------------------------------------ #
java -jar -Xmx32g  ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar SelectValidCellBarcode I=BarcodesAssigned.tsv O=ValidBarcodes.csv MINUMI=1 ED0ED1RATIO=1

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.primary_assembly.annotation.gtf.gz
        
gunzip gencode.vM38.primary_assembly.annotation.gtf.gz
gtfToGenePred -genePredExt -geneNameAsName2 gencode.vM38.primary_assembly.annotation.gtf gencode.vM38.primary_assembly.annotation.refflat.txt
paste <(cut -f 12 gencode.vM38.primary_assembly.annotation.refflat.txt) <(cut -f 1-10 gencode.vM38.primary_assembly.annotation.refflat.txt) > gencode.vM38.refFlat

java -jar -Xmx64g  ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar IsoformMatrix I=adeno1_sicelore_assigned.bam  GENETAG=GE UMITAG=U8 CELLTAG=BC REFFLAT=gencode.vM38.refFlat CSV=barcodes.csv DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false OUTDIR=. PREFIX=sicelore




## ------------------------------------------ #
## ------------------------------------------ #
## ------------------------------------------ #
## ------------------------------------------ #
## ------------------------------------------ #
## ------------------------------------------ #
# SHAM RUN 1
java -jar -Xmx40g ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar scanfastq -d ../adeno1_run1_fastq_pass/ -o /home/lab-user/data/sham1run1_sicelore --bcEditDistance 1
BUILD=~/reference/refdata-gex-GRCm39-2024-A/fasta/genome.fa
minimap2 -ax splice -uf --sam-hit-only -t 40 $BUILD passed/*passed.fastq.gz | samtools view -bS -@ 20 - | samtools sort -m 2G -@ 20 -o sham1run1_passed.bam -&& samtools index sham1run1_passed.bam
java -jar -Xmx40g ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar assignumis --inFileNanopore sham1run1_passed.bam --outfile sham1run1_sicelore_assigned.bam --ONTgene GE -a ../../reference/gencode.vM38.primary_assembly.annotation.gtf
java -jar -Xmx32g  ../../tools/sicelore-2.1/Jar/Sicelore-2.1.jar SelectValidCellBarcode I=BarcodesAssigned.tsv O=ValidBarcodes.csv MINUMI=10 ED0ED1RATIO=2
java -jar -Xmx64g  ../../tools/sicelore-2.1/Jar/Sicelore-2.1.jar  IsoformMatrix I=sham1run1_sicelore_assigned.bam  GENETAG=GE UMITAG=U8 CELLTAG=BC REFFLAT=gencode.vM38.refFlat CSV=barcodes.csv DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false OUTDIR=. PREFIX=sicelore


## ------------------------------------------ #
# ADENO SC4 RUN
java -jar -Xmx40g ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar scanfastq -d ../adeno1_run1_fastq_pass/ -o /home/lab-user/data/adeno_sc4_sicelore --bcEditDistance 1
BUILD=~/reference/refdata-gex-GRCm39-2024-A/fasta/genome.fa
minimap2 -ax splice -uf --sam-hit-only -t 40 $BUILD passed/*passed.fastq.gz | samtools view -bS -@ 20 - | samtools sort -m 2G -@ 20 -o adeno_sc4_passed.bam -&& samtools index adeno_sc4_passed.bam
java -jar -Xmx40g ../../tools/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar assignumis --inFileNanopore adeno_sc4_passed.bam --outfile adeno_sc4_sicelore_assigned.bam --ONTgene GE -a ../../reference/gencode.vM38.primary_assembly.annotation.gtf
java -jar -Xmx32g  ../../tools/sicelore-2.1/Jar/Sicelore-2.1.jar SelectValidCellBarcode I=BarcodesAssigned.tsv O=ValidBarcodes.csv MINUMI=10 ED0ED1RATIO=2
java -jar -Xmx64g  ../../tools/sicelore-2.1/Jar/Sicelore-2.1.jar  IsoformMatrix I=adeno_sc4_sicelore_assigned.bam  GENETAG=GE UMITAG=U8 CELLTAG=BC REFFLAT=gencode.vM38.refFlat CSV=barcodes.csv DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false OUTDIR=. PREFIX=sicelore

