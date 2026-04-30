
# from Nanopore merged fastq directory
zcat adeno.combined.fastq.gz > adeno.combined.fastq
zcat sham.combined.fastq.gz > sham.combined.fastq

# Sudo required until we can set up permissions for the nextflow output directories
sudo nextflow run epi2me-labs/wf-single-cell \
    --expected_cells 10000 \
    --fastq adeno1_merged.fastq  \
    --kit 3prime:v3 \
    --ref_genome_dir ~/reference/refdata-gex-GRCm39-2024-A \
    --out_dir adeno1_merged_10000_wf \
    -resume \
    --threads 30 \
    --matrix_max_mito 20 \
    --mito_prefix mt- \
    -process.executor local \
    -process.maxForks 12 \
    -process.cpus 4 \
    -process.memory 16.GB

sudo nextflow run epi2me-labs/wf-single-cell \
    --expected_cells 5000 \
    --fastq adeno1_merged.fastq  \
    --kit 3prime:v3 \
    --ref_genome_dir ~/reference/refdata-gex-GRCm39-2024-A \
    --out_dir adeno1_merged_5000_wf \
    -resume \
    --threads 30 \
    --matrix_max_mito 20 \
    --mito_prefix mt- \
    -process.executor local \
    -process.maxForks 12 \
    -process.cpus 4 \
    -process.memory 16.GB

sudo nextflow run epi2me-labs/wf-single-cell \
    --expected_cells 5000 \
    --fastq adeno1_merged.fastq  \
    --kit 3prime:v3 \
    --ref_genome_dir ~/reference/refdata-gex-GRCm39-2024-A \
    --out_dir adeno1_merged_5000_wf \
    -resume \
    --threads 30 \
    --matrix_max_mito 20 \
    --mito_prefix mt- \
    -process.executor local \
    -process.maxForks 12 \
    -process.cpus 4 \
    -process.memory 16.GB

sudo nextflow run epi2me-labs/wf-single-cell \
    --expected_cells 5000 \
    --fastq adeno1_merged.fastq  \
    --kit 3prime:v3 \
    --ref_genome_dir ~/reference/refdata-gex-GRCm39-2024-A \
    --out_dir adeno1_merged_5000_wf \
    -resume \
    --threads 30 \
    --matrix_max_mito 20 \
    --mito_prefix mt- \
    -process.executor local \
    -process.maxForks 12 \
    -process.cpus 4 \
    -process.memory 16.GB




