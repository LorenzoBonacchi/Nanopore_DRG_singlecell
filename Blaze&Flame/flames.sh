# In blaze output dir
# with a copy of the reference genome and annotation files --> genome.fa and genes.gtf

docker run -it \
  -v $PWD:/data \
  ghcr.io/mritchielab/flames:20af1ce


library(FLAMES)
setwd("/data")
outdir <- "flames_out"
dir.create(outdir, showWarnings = FALSE)
config_file <- create_config(
  outdir = outdir,
  type = "sc_3end",
  threads = 30,
  do_barcode_demultiplex = FALSE,
  barcode_parameters.max_bc_editdistance = 3,
  isoform_parameters.min_sup_cnt = 10,
  multithread_isoform_identification = TRUE,
  oarfish_quantification = TRUE,
  additional_arguments.oarfish = c("--model-coverage")
)

sce <- sc_long_pipeline(
  annotation = "genes.gtf",
  fastq = "matched_reads.fastq.gz",
  genome_fa = "genome.fa",
  outdir = outdir,
  barcodes_file = "whitelist.csv",
  config_file = config_file
)

# ------------------------
# Started at 9:30am circa


