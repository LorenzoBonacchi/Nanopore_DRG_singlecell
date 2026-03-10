
# from Nanopore merged fastq directory
zcat adeno.combined.fastq.gz > adeno.combined.fastq
zcat sham.combined.fastq.gz > sham.combined.fastq

# Sudo required until we can set up permissions for the nextflow output directories
sudo nextflow run epi2me-labs/wf-single-cell \
    --expected_cells 8000 \
    --fastq adeno.combined.fastq \
    --kit 3prime:v3 \
    --ref_genome_dir ~/reference/refdata-gex-GRCm39-2024-A \
    --out_dir adeno_wf \
    -resume \
    -process.executor local \
    -process.maxForks 12 \
    -process.cpus 4 \
    -process.memory 16.GB

sudo nextflow run epi2me-labs/wf-single-cell \
    --expected_cells 8000 \
    --fastq sham.combined.fastq \
    --kit 3prime:v3 \
    --ref_genome_dir ~/reference/refdata-gex-GRCm39-2024-A \
    --out_dir sham_wf \
    -resume \
    -process.executor local \
    -process.maxForks 12 \
    -process.cpus 4 \
    -process.memory 16.GB
