## A549 vs primary cells (HBE, ALI) - IL1B, GC, and combination

ssh alex.gao1@arc.ucalgary.ca

### setup ARC
```bash
salloc -c 4 --mem 32GB --time 05:00:00
module load R/4.4.1
module load kallisto/0.46.1

MULTIQC="/home/alex.gao1/tools/multiqc_latest.sif"
#v1.19

export PATH=/home/alex.gao1/tools/FastQC:$PATH
# FastQC v0.12.1

mkdir -p /work/newton_lab/ag_analysis/a549-hbe-ali-comparison
cd /work/newton_lab/ag_analysis/a549-hbe-ali-comparison
```


### import metafiles
```bash
scp ./data/meta_a549_il1b_bud.txt alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/a549-hbe-ali-comparison
scp ./data//meta_ali_il1b_bud.txt alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/a549-hbe-ali-comparison
scp ./data//meta_hbe_il1b_dex.txt alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/a549-hbe-ali-comparison

dos2unix ./*.txt
```

### run QC on raw fastq files

```bash
mkdir -p fastqc
cd fastqc

# build list of fastq files to process
awk '
FNR==1 {next}
{
  base=$6
  for(i=7;i<=10;i++){
    if($i ~ /\.fastq\.gz$/){
      print base $i
    }
  }
}
' ../meta_*.txt | sort -u > fastq_files.txt

# prepare slurm script
cat <<'EOF' > run_fastqc_array.slurm
#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --output=/work/newton_lab/ag_analysis/a549-hbe-ali-comparison/fastqc/logs/%A_%a.out
#SBATCH --error=/work/newton_lab/ag_analysis/a549-hbe-ali-comparison/fastqc/logs/%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1

export PATH=/home/alex.gao1/tools/FastQC:$PATH
# FastQC v0.12.1

cd /work/newton_lab/ag_analysis/a549-hbe-ali-comparison/fastqc/

FASTQ_LIST=fastq_files.txt
OUTDIR=fastqc_results

mkdir -p "${OUTDIR}"
mkdir -p logs

# get this task's FASTQ
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FASTQ_LIST}")

echo "Processing ${FASTQ}"

fastqc -t ${SLURM_CPUS_PER_TASK} -o "${OUTDIR}" "${FASTQ}"

EOF

# submit fastqc runs
N=$(wc -l < fastq_files.txt)
sbatch --array=1-${N}%40 run_fastqc_array.slurm

# run multiqc

apptainer exec --bind /work:/work "$MULTIQC" multiqc ../fastqc/fastqc_results -o multiqc
```


### setup kallisto jobs
```bash
mkdir -p kallisto
mkdir -p kallisto/logs

# run_kallisto_singleEnd.slurm takes the path of $1 fastq file and outputs to $2
# don't actually have to run this, kallisto files are already in MDSC519_project/kallisto/tc
cat <<'EOF' > kallisto/run_kallisto_singleEnd.slurm
#!/bin/bash

#SBATCH --job-name=kallisto
#SBATCH --output=/work/newton_lab/ag_analysis/a549-hbe-ali-comparison/kallisto/logs/%j_kallistoSingle.out
#SBATCH --error=/work/newton_lab/ag_analysis/a549-hbe-ali-comparison/kallisto/logs/%j_kallistoSingle.err
#SBATCH --time=03:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

module load kallisto/0.46.1
cd /work/newton_lab/ag_analysis/a549-hbe-ali-comparison/

index=/work/newton_lab/ag_analysis/ref_seq/hg38_p14_kallisto.idx
fastq1=$1
output=$2
mkdir -p "$(dirname "$output")"

kallisto quant \
  -i "$index" \
  -o "$output" \
  --single -l 180 -s 20 \
  --threads 4 \
  -b 100 \
  "$fastq1"
EOF

# run_kallisto_pairedEnd.slurm takes the path of $1, $2, $3, $4 fastq files and outputs to $5
cat <<'EOF' > kallisto/run_kallisto_pairedEnd.slurm
#!/bin/bash

#SBATCH --job-name=kallisto
#SBATCH --output=/work/newton_lab/ag_analysis/a549-hbe-ali-comparison/kallisto/logs/%j_kallistoPaired.out
#SBATCH --error=/work/newton_lab/ag_analysis/a549-hbe-ali-comparison/kallisto/logs/%j_kallistoPaired.err
#SBATCH --time=03:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

module load kallisto/0.46.1
cd /work/newton_lab/ag_analysis/a549-hbe-ali-comparison/

index=/work/newton_lab/ag_analysis/ref_seq/hg38_p14_kallisto.idx
fastq1=$1
fastq2=$2
fastq3=$3
fastq4=$4
output=$5
mkdir -p "$(dirname "$output")"

kallisto quant \
    -i "$index" \
    -o "$output" \
    --threads 4 \
    -b 100 \
    "$fastq1" "$fastq2" "$fastq3" "$fastq4"
EOF
```


### run kallisto jobs
```bash
# a549 jobs
# TODO
# technically don't need to rerun as the metafile error (fixed now) only affected the 12h timepoint

# ali jobs
awk -F'\t' 'NR>1 {
  cmd = "sbatch kallisto/run_kallisto_pairedEnd.slurm " \
        $6 $7 " " \
        $6 $8 " " \
        $6 $9 " " \
        $6 $10 " " \
        $11
  system(cmd)
}' meta_ali_il1b_bud.txt

# hbe jobs
awk -F'\t' 'NR>1 {
  cmd = "sbatch kallisto/run_kallisto_pairedEnd.slurm " \
        $6 $7 " " \
        $6 $8 " " \
        $6 $9 " " \
        $6 $10 " " \
        $11
  system(cmd)
}' meta_hbe_il1b_dex.txt
```


### run QC on kallisto output
```bash
apptainer exec --bind /work:/work "$MULTIQC" multiqc kallisto/logs/ -o kallisto_multiqc
```

```bash (powershell)
scp alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/a549-hbe-ali-comparison/multiqc/multiqc_report.html ./data
```


### modelling and DEA work
```bash
mkdir -p sleuth
R
```

```r setup
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/a549-hbe-ali-comparison")

library(dplyr)
library(readr)
library(sleuth)
library(stringr)

system("mkdir sleuth/objects")
system("ls sleuth")
```

```r prepare sleuth inputs
# transcript ids to gene names
t2g <- read_tsv("/work/newton_lab/ag_analysis/ref_seq/Homo_sapiens.GRCh38.p14.cdna.all_mart_export.txt") %>%
  select(target_id = `Transcript stable ID`, Gene = `Gene name`) %>%
  na.omit() %>% 
  distinct()

#filter requiring at least 5 reads in at least 20% of samples
new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

# full s2c for tpms
s2c <- rbind(
        read_tsv("meta_a549_il1b_bud.txt"),
        read_tsv("meta_ali_il1b_bud.txt"),
        read_tsv("meta_hbe_il1b_dex.txt")) %>%
    mutate(
        treatment = if_else(treatment == "Dex" | treatment == "Bud", "GC", treatment),
        treatment = if_else(treatment == "I+D" | treatment == "I+B", "combo", treatment)) %>%
    mutate(
        treatment = factor(treatment, levels = c("NS", "IL1B", "GC", "combo")),
        celltype = factor(celltype, levels = c("A549", "ALI", "HBE"))
    ) %>%
    select(sample, rep, celltype, treatment, time, path = kallisto_path)
    
# s2c at 6h for comparisons
s2c_6h <- filter(s2c, time == 6)

# s2c at 6h for each celltype
s2c_a549_6h <- filter(s2c, time == 6 & celltype == "A549")
s2c_ali_6h <- filter(s2c, time == 6 & celltype == "ALI")
s2c_hbe_6h <- filter(s2c, time == 6 & celltype == "HBE")

# only il1b at 6h. also specify donors
s2c_il1b_6h = s2c %>%
    filter(time == 6 & treatment %in% c("NS", "IL1B")) %>%
    mutate(
        donor = paste0(celltype, rep),
        donor = if_else(str_detect(celltype, "A549"), "A549donor", donor),
        .before = celltype
    )

# s2c at 6h for each celltype, IL1B only, specify donors
s2c_a549_il1b_6h <- s2c %>%
  filter(time == 6, celltype == "A549", treatment %in% c("NS", "IL1B")) %>%
  mutate(donor = "A549donor")

s2c_ali_il1b_6h <- s2c %>%
  filter(time == 6, celltype == "ALI", treatment %in% c("NS", "IL1B")) %>%
  mutate(donor = rep)

s2c_hbe_il1b_6h <- s2c %>%
  filter(time == 6, celltype == "HBE", treatment %in% c("NS", "IL1B")) %>%
  mutate(donor = rep)

# add more s2cs here
```

```r ~treatment model, each celltype at 6h
# a549
so_a549_6h = sleuth_prep(
    sample_to_covariates = s2c_a549_6h,
    full_model = ~treatment,
    target_mapping = t2g,
    gene_mode = TRUE, aggregation_column = "Gene",
    filter_fun = new_filter, num_cores = 4
)

so_a549_6h = sleuth_fit(so_a549_6h, ~treatment) %>%
    sleuth_wt("treatmentIL1B") %>%
    sleuth_wt("treatmentGC") %>%
    sleuth_wt("treatmentcombo")

# ali
so_ali_6h = sleuth_prep(
    sample_to_covariates = s2c_ali_6h,
    full_model = ~treatment,
    target_mapping = t2g,
    gene_mode = TRUE, aggregation_column = "Gene",
    filter_fun = new_filter, num_cores = 1
)

so_ali_6h = sleuth_fit(so_ali_6h, ~treatment) %>%
    sleuth_wt("treatmentIL1B") %>%
    sleuth_wt("treatmentGC") %>%
    sleuth_wt("treatmentcombo")

# hbe
so_hbe_6h = sleuth_prep(
    sample_to_covariates = s2c_hbe_6h,
    full_model = ~treatment,
    target_mapping = t2g,
    gene_mode = TRUE, aggregation_column = "Gene",
    filter_fun = new_filter, num_cores = 4
)

so_hbe_6h = sleuth_fit(so_hbe_6h, ~treatment) %>%
    sleuth_wt("treatmentIL1B") %>%
    sleuth_wt("treatmentGC") %>%
    sleuth_wt("treatmentcombo")

# save
saveRDS(so_a549_6h, "sleuth/objects/so_a549_6h.rds")
saveRDS(so_ali_6h, "sleuth/objects/so_ali_6h.rds")
saveRDS(so_hbe_6h, "sleuth/objects/so_hbe_6h.rds")

# pull data
a <- sleuth_results(so_a549_6h, "treatmentIL1B") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "A549", treatment = "IL1B", time = "6h", .before = 2)

b <- sleuth_results(so_a549_6h, "treatmentGC") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "A549", treatment = "GC", time = "6h", .before = 2)

c <- sleuth_results(so_a549_6h, "treatmentcombo") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "A549", treatment = "combo", time = "6h", .before = 2)

d <- sleuth_results(so_ali_6h, "treatmentIL1B") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "ALI", treatment = "IL1B", time = "6h", .before = 2)

e <- sleuth_results(so_ali_6h, "treatmentGC") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "ALI", treatment = "GC", time = "6h", .before = 2)

f <- sleuth_results(so_ali_6h, "treatmentcombo") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "ALI", treatment = "combo", time = "6h", .before = 2)

g <- sleuth_results(so_hbe_6h, "treatmentIL1B") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "HBE", treatment = "IL1B", time = "6h", .before = 2)

h <- sleuth_results(so_hbe_6h, "treatmentGC") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "HBE", treatment = "GC", time = "6h", .before = 2)

i <- sleuth_results(so_hbe_6h, "treatmentcombo") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "HBE", treatment = "combo", time = "6h", .before = 2)

results_individual = rbind(a,b,c,d,e,f,g,h,i) %>%
    mutate(treatment = factor(treatment, levels = c("IL1B", "GC", "combo"))) %>%
    arrange(Gene, celltype, treatment) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    group_by(Gene) %>%
    filter(!all(is.na(log2fold))) %>%
    ungroup() %>%
    mutate(
        log2fold = if_else(is.na(log2fold), 0, log2fold),
        FDR = if_else(is.na(FDR), 1, FDR))

write_tsv(results_individual, "sleuth/A549vsPrimary_univariateDEA.txt")
```

```r ~treatment+donor model, each celltype at 6h
so_a549_il1b_6h = sleuth_prep(
    sample_to_covariates = s2c_a549_il1b_6h,
    full_model = ~ treatment,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter_fun = new_filter,
    num_cores = 4
)

so_ali_il1b_6h <- sleuth_prep(
  sample_to_covariates = s2c_ali_il1b_6h,
  full_model = ~ donor + treatment,
  target_mapping = t2g,
  gene_mode = TRUE,
  aggregation_column = "Gene",
  filter_fun = new_filter,
  num_cores = 4
)

so_hbe_il1b_6h <- sleuth_prep(
  sample_to_covariates = s2c_hbe_il1b_6h,
  full_model = ~ donor + treatment,
  target_mapping = t2g,
  gene_mode = TRUE,
  aggregation_column = "Gene",
  filter_fun = new_filter,
  num_cores = 4
)

so_a549_il1b_6h = sleuth_fit(so_a549_il1b_6h, ~treatment) %>%
    sleuth_wt("treatmentIL1B")

so_ali_il1b_6h = sleuth_fit(so_ali_il1b_6h, ~donor + treatment) %>%
    sleuth_wt("treatmentIL1B")

so_hbe_il1b_6h = sleuth_fit(so_hbe_il1b_6h, ~donor + treatment) %>%
    sleuth_wt("treatmentIL1B")

a <- sleuth_results(so_a549_il1b_6h, "treatmentIL1B") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "A549", treatment = "IL1B", time = "6h", .before = 2)
b <- sleuth_results(so_ali_il1b_6h, "treatmentIL1B") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "ALI", treatment = "IL1B", time = "6h", .before = 2)
c <- sleuth_results(so_hbe_il1b_6h, "treatmentIL1B") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(celltype = "HBE", treatment = "IL1B", time = "6h", .before = 2)

results_focused = rbind(a,b,c) %>%
    arrange(Gene, celltype) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    group_by(Gene) %>%
    filter(!all(is.na(log2fold))) %>%
    ungroup() %>%
    mutate(
        log2fold = if_else(is.na(log2fold), 0, log2fold),
        FDR = if_else(is.na(FDR), 1, FDR))

write_tsv(results_focused, "sleuth/A549vsPrimary_DEA_donorCorrected_IL1Bfocused.txt")
```

```r all celltypes at 6h, pull TPMs
so_6h = sleuth_prep(
    sample_to_covariates = s2c_6h,
    full_model = ~treatment*celltype,
    target_mapping = t2g,
    gene_mode = TRUE, aggregation_column = "Gene",
    filter_fun = new_filter, num_cores = 4
    )

# # define the universe of differing responses
# so_6h = sleuth_fit(so_6h, ~celltype*treatment, "full")
# so_6h = sleuth_fit(so_6h, ~celltype+treatment, "reduced")
# so_6h = sleuth_lrt(so_6h, "reduced", "full")

# # run wald tests
# so_6h = sleuth_wt(so_6h, "celltypeALI:treatmentIL1B")
# so_6h = sleuth_wt(so_6h, "celltypeHBE:treatmentIL1B")
# so_6h = sleuth_wt(so_6h, "celltypeALI:treatmentGC")
# so_6h = sleuth_wt(so_6h, "celltypeHBE:treatmentGC")
# so_6h = sleuth_wt(so_6h, "celltypeALI:treatmentcombo")
# so_6h = sleuth_wt(so_6h, "celltypeHBE:treatmentcombo")

# save
saveRDS(so_6h, "sleuth/objects/so_6h.rds")

# # pull results
# a <- sleuth_results(so_6h, "reduced:full", test_type = "lrt") %>%
#     select(Gene = target_id, lrt_FDR = qval)
# b <- sleuth_results(so_6h, "celltypeALI:treatmentIL1B", test_type = "wt") %>%
#     mutate(b = b / log(2)) %>%
#     select(Gene = target_id, ALI_IL1B_diff = b, ALI_IL1B_FDR = qval)
# c <- sleuth_results(so_6h, "celltypeHBE:treatmentIL1B", test_type = "wt") %>%
#     mutate(b = b / log(2)) %>%
#     select(Gene = target_id, HBE_IL1B_diff = b, HBE_IL1B_FDR = qval)
# d <- sleuth_results(so_6h, "celltypeALI:treatmentGC", test_type = "wt") %>%
#     mutate(b = b / log(2)) %>%
#     select(Gene = target_id, ALI_GC_diff = b, ALI_GC_FDR = qval)
# e <- sleuth_results(so_6h, "celltypeHBE:treatmentGC", test_type = "wt") %>%
#     mutate(b = b / log(2)) %>%
#     select(Gene = target_id, HBE_GC_diff = b, HBE_GC_FDR = qval)
# f <- sleuth_results(so_6h, "celltypeALI:treatmentcombo", test_type = "wt") %>%
#     mutate(b = b / log(2)) %>%
#     select(Gene = target_id, ALI_combo_diff = b, ALI_combo_FDR = qval)
# g <- sleuth_results(so_6h, "celltypeHBE:treatmentcombo", test_type = "wt") %>%
#     mutate(b = b / log(2)) %>%
#     select(Gene = target_id, HBE_combo_diff = b, HBE_combo_FDR = qval)

# results_6h <- a %>%
#     left_join(b, by = "Gene") %>%
#     left_join(c, by = "Gene") %>%
#     left_join(d, by = "Gene") %>%
#     left_join(e, by = "Gene") %>%
#     left_join(f, by = "Gene") %>%
#     left_join(g, by = "Gene") %>%
#     arrange(Gene) %>%
#     filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
#     filter(!if_all(-Gene, is.na)) %>%
#     mutate(
#         across(ends_with("_diff"), ~ifelse(is.na(.), 0, .)),
#         across(ends_with("_FDR"),  ~ifelse(is.na(.), 1, .)),
#         lrt_FDR = ifelse(is.na(lrt_FDR), 1, lrt_FDR))

tpm_6h <- kallisto_table(so_6h, use_filtered = FALSE) %>%
    select(Gene = target_id, rep, celltype, treatment, time, tpm) %>%
    arrange(Gene, celltype, treatment) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    group_by(Gene) %>%
    filter(!all(tpm == 0)) %>%
    ungroup() %>%
    mutate(log2tpm = log2(tpm + 0.1)) %>%
    group_by(Gene, celltype) %>%
    mutate(log2fold = log2tpm - log2tpm[treatment == "NS"]) %>%
    ungroup() %>%
    mutate(fold = 2^log2fold) %>%
    select(Gene, rep, celltype, treatment, time, tpm, log2tpm, fold, log2fold)

# write_tsv(results_6h, "sleuth/A549vsPrimary_multivariateDEA.txt")
write_tsv(tpm_6h,     "sleuth/A549vsPrimary_tpm.txt")
```

```r celltype/treatment LRT for IL1B
so_il1b_6h = sleuth_prep(
    sample_to_covariates = s2c_il1b_6h,
    target_mapping = t2g,
    gene_mode = TRUE, aggregation_column = "Gene",
    filter_fun = new_filter, num_cores = 4
)

saveRDS(so_il1b_6h, "sleuth/objects/so_il1b_6h.rds")

# reduced model: baseline differences between celltypes, one global IL1B response
so_il1b_6h <- sleuth_fit(
    so_il1b_6h,
    ~ celltype + treatment,
    'reduced'
)

# full model: still baseline differences between celltypes, but IL1B response can differ by celltype
so_il1b_6h <- sleuth_fit(
    so_il1b_6h,
    ~ celltype + treatment + celltype:treatment,
    'full'
)

# LRT: does allowing IL1B responses to differ by cell type significantly improve fit compared to a single shared IL1B response?
so_il1b_6h <- sleuth_lrt(
    so_il1b_6h,
    "reduced",
    "full"
)

lrt_results = sleuth_results(so_il1b_6h, "reduced:full", test_type = "lrt") %>%
    rename(Gene = target_id) %>%
    rename(FDR = qval)

write_tsv(lrt_results, "sleuth/A549vsPrimary_lrt_IL1B_6h.txt")
```

```r donor./treatment LRT for IL1B
# so_il1b_6h = sleuth_prep(
#     sample_to_covariates = s2c_il1b_6h,
#     target_mapping = t2g,
#     gene_mode = TRUE, aggregation_column = "Gene",
#     filter_fun = new_filter, num_cores = 4
# )

# saveRDS(so_il1b_6h, "sleuth/objects/so_il1b_6h.rds")

# reduced model: baseline differences between donors, one global IL1B response
so_il1b_6h = sleuth_fit(
    so_il1b_6h,
    ~ donor + treatment,
    'reduced'
)

# full model: still baseline differences between donors, but IL1B response can differ by donor
so_il1b_6h <- sleuth_fit(
    so_il1b_6h,
    ~ donor + treatment + donor:treatment,
    'full'
)

# LRT: does allowing IL1B responses to differ by donor significantly improve fit compared to a single shared IL1B response?
so_il1b_6h <- sleuth_lrt(
    so_il1b_6h,
    "reduced",
    "full"
)

lrt2_results = sleuth_results(so_il1b_6h, "reduced:full", test_type = "lrt") %>%
    rename(Gene = target_id) %>%
    rename(FDR = qval)

write_tsv(lrt2_results, "sleuth/A549vsPrimary_lrt2_IL1B_6h.txt")
```