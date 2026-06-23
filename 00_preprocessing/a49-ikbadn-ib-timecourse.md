## A549 under naive/GFP/IKBADN +- IL1B/Bud/both

ssh alex.gao1@arc.ucalgary.ca

### setup ARC
```bash
salloc -c 8 --mem 128GB --time 05:00:00
module load R/4.4.1
module load kallisto/0.46.1

MULTIQC="/home/alex.gao1/tools/multiqc_latest.sif"
#v1.19

export PATH=/home/alex.gao1/tools/FastQC:$PATH
# FastQC v0.12.1

mkdir -p /work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse
cd /work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse
export WORK_DIR="/work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse"
```

### import metafiles
```bash
scp ./data/meta_a549_ikbadn_ib.txt alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse
dos2unix ./*.txt
```

### run QC on raw fastq files

```bash
mkdir -p $WORK_DIR/fastqc
cd $WORK_DIR/fastqc

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
#SBATCH --output=/work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse/fastqc/logs/%A_%a.out
#SBATCH --error=/work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse/fastqc/logs/%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1

export PATH=/home/alex.gao1/tools/FastQC:$PATH
# FastQC v0.12.1

cd /work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse/fastqc/

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
mkdir -p $WORK_DIR/fastqc
cd $WORK_DIR/fastqc
apptainer exec --bind /work:/work "$MUlsLTIQC" multiqc ../fastqc/fastqc_results -o multiqc
```


### setup kallisto jobs
```bash
mkdir -p $WORK_DIR/kallisto
mkdir -p $WORK_DIR/kallisto/logs
cd $WORK_DIR/kallisto

# run_kallisto_pairedEnd.slurm takes the path of $1, $2 fastq files and outputs to $3
cat <<'EOF' > run_kallisto_pairedEnd.slurm
#!/bin/bash

#SBATCH --job-name=kallisto
#SBATCH --output=/work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse/kallisto/logs/%j_kallistoPaired.out
#SBATCH --error=/work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse/kallisto/logs/%j_kallistoPaired.err
#SBATCH --time=03:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

module load kallisto/0.46.1
cd /work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse/

index=/work/newton_lab/ag_analysis/ref_seq/hg38_p14_kallisto.idx
fastq1=$1
fastq2=$2
output=$3
mkdir -p "$(dirname "$output")"

kallisto quant \
    -i "$index" \
    -o "$output" \
    --threads 4 \
    -b 100 \
    "$fastq1" "$fastq2"
EOF
```

### run kallisto jobs
```bash
awk -F'\t' 'NR>1 {
  cmd = "sbatch run_kallisto_pairedEnd.slurm " \
        $6 $7 " " \
        $6 $8 " " \
        $11
  system(cmd)
}' $WORK_DIR/meta_a549_ikbadn_ib.txt
```

### run QC on kallisto output
```bash
cd $WORK_DIR/kallisto
apptainer exec --bind /work:/work "$MULTIQC" multiqc $WORK_DIR/kallisto/logs/
```


### modelling and DEA work

```bash
mkdir -p $WORK_DIR/sleuth
mkdir -p $WORK_DIR/sleuth/objects
cd $WORK_DIR

R
```

```r setup
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/a549-ikbadn-ib-timecourse")

library(dplyr)
library(readr)
library(sleuth)
library(tidyr)
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
s2c = read_tsv("meta_a549_ikbadn_ib.txt") %>%
    mutate(
        treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "IB")),
        condition = factor(condition, levels = c("naive", "Ad-GFP", "Ad-IKBA"))
    ) %>%
    select(sample, rep, condition, treatment, time, path = kallisto_path)

# s2c at timepoints for multivariate stuff
s2c_1h = filter(s2c, time == 1)
s2c_2h = filter(s2c, time == 2)
s2c_6h = filter(s2c, time == 6)

# to do: individual timepoint x condition for pure treatment testing

```

```r prepare sleuth objects
# full sleuth object for pulling tpms 
so = sleuth_prep(
    sample_to_covariates = s2c,
    full_model = ~treatment*condition,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

# time point sleuth objects for multivariate analysis
so_1h = sleuth_prep(
    sample_to_covariates = s2c_1h,
    full_model = ~treatment*condition,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_2h = sleuth_prep(
    sample_to_covariates = s2c_2h,
    full_model = ~treatment*condition,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

so_6h = sleuth_prep(
    sample_to_covariates = s2c_6h,
    full_model = ~treatment*condition,
    target_mapping = t2g,
    gene_mode = TRUE,
    aggregation_column = "Gene",
    filter = new_filter,
    num_cores = 4
)

# single time point and condition objects for univariate analysis
univariate_sleuth_objs = list()

for(t in c(1, 2, 6)){
    for(c in c("naive", "Ad-GFP", "Ad-IKBA")){

        s = s2c %>%
            filter(time == t & condition == c)

        print(paste0("Processing ", name, "..."))

        name = paste0("so_", c, "_", t, "h")

        o = sleuth_prep(
            sample_to_covariates = s,
            full_model = ~treatment,
            target_mapping = t2g,
            gene_mode = TRUE,
            aggregation_column = "Gene",
            filter = new_filter,
            num_cores = 4
        )

        univariate_sleuth_objs[[name]] = o

    }
}
```


```r pull tpms
tpm = kallisto_table(so, use_filtered = FALSE)

tpm_clean = tpm %>%
    select(Gene = target_id, rep, condition, treatment, time, tpm) %>%
    arrange(Gene, condition, treatment) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    group_by(Gene) %>%
    filter(!all(tpm == 0)) %>%
    ungroup() %>%
    mutate(log2tpm = log2(tpm + 0.1)) %>%
    group_by(Gene, time, condition) %>%
    mutate(log2fold = log2tpm - log2tpm[treatment == "NS"]) %>%
    ungroup() %>%
    mutate(fold = 2^log2fold) %>%
    arrange(Gene, time, condition, treatment) %>%
    select(Gene, rep, condition, treatment, time, tpm, log2tpm, fold, log2fold)

write_tsv(tpm_clean, "sleuth/a549-ikbadn-ib-timecourse_tpm.txt")
```


```r multivariate modelling
run_multivar_sleuth_tests <- function(so) {
  
  so %>%
            # test for condition-treatment interactions
    sleuth_fit(~condition*treatment, "full") %>%
    sleuth_fit(~condition+treatment, "reduced") %>%
    sleuth_lrt("reduced", "full") %>%
            # test pairs
    sleuth_wt("conditionAd-GFP:treatmentIL1B") %>%
    sleuth_wt("conditionAd-GFP:treatmentBud") %>%
    sleuth_wt("conditionAd-GFP:treatmentIB") %>%
    sleuth_wt("conditionAd-IKBA:treatmentIL1B") %>%
    sleuth_wt("conditionAd-IKBA:treatmentBud") %>%
    sleuth_wt("conditionAd-IKBA:treatmentIB")
}

so_1h <- run_multivar_sleuth_tests(so_1h)
so_2h <- run_multivar_sleuth_tests(so_2h)
so_6h <- run_multivar_sleuth_tests(so_6h)
```

```r pull multivariate results

get_multivar_results_table <- function(so, timepoint){

    a = sleuth_results(so, "reduced:full", test_type = "lrt") %>%
        select(Gene = target_id, lrt_FDR = qval)

    b = sleuth_results(so, "conditionAd-GFP:treatmentIL1B", test_type = "wt") %>% 
        mutate(b = b/log(2)) %>%
        select(Gene = target_id, GFP_IL1B_diff = b, GFP_IL1B_FDR = qval)

    c = sleuth_results(so, "conditionAd-GFP:treatmentBud", test_type = "wt") %>% 
        mutate(b = b/log(2)) %>%
        select(Gene = target_id, GFP_Bud_diff = b, GFP_Bud_FDR = qval)

    d = sleuth_results(so, "conditionAd-GFP:treatmentIB", test_type = "wt") %>% 
        mutate(b = b/log(2)) %>%
        select(Gene = target_id, GFP_IB_diff = b, GFP_IB_FDR = qval)

    e = sleuth_results(so, "conditionAd-IKBA:treatmentIL1B", test_type = "wt") %>% 
        mutate(b = b/log(2)) %>%
        select(Gene = target_id, IKBA_IL1B_diff = b, IKBA_IL1B_FDR = qval)

    f = sleuth_results(so, "conditionAd-IKBA:treatmentBud", test_type = "wt") %>% 
        mutate(b = b/log(2)) %>%
        select(Gene = target_id, IKBA_Bud_diff = b, IKBA_Bud_FDR = qval)

    g = sleuth_results(so, "conditionAd-IKBA:treatmentIB", test_type = "wt") %>% 
        mutate(b = b/log(2)) %>%
        select(Gene = target_id, IKBA_IB_diff = b, IKBA_IB_FDR = qval)

    results <- a %>%
        left_join(b, by = "Gene") %>%
        left_join(c, by = "Gene") %>%
        left_join(d, by = "Gene") %>%
        left_join(e, by = "Gene") %>%
        left_join(f, by = "Gene") %>%
        left_join(g, by = "Gene") %>%
        arrange(Gene) %>%
        filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
        filter(!if_all(-Gene, is.na)) %>%
        mutate(
            across(ends_with("_diff"), ~ifelse(is.na(.), 0, .)),
            across(ends_with("_FDR"),  ~ifelse(is.na(.), 1, .)),
            lrt_FDR = ifelse(is.na(lrt_FDR), 1, lrt_FDR),
            time = timepoint
        )

    results
}

results_1h <- get_multivar_results_table(so_1h, "1h")
results_2h <- get_multivar_results_table(so_2h, "2h")
results_6h <- get_multivar_results_table(so_6h, "6h")

multivar = rbind(results_1h, results_2h, results_6h) %>%
                select(Gene, time, everything()) %>%
                arrange(Gene, time) %>%
                as_tibble()

write_tsv(multivar, "sleuth/a549-ikbadn-ib-timecourse_multivarDEA.txt")
```


```r univariate modelling

for (nm in names(univariate_sleuth_objs)){
    
    print(paste0("Processing ", nm, "..."))

    univariate_sleuth_objs[[nm]] = univariate_sleuth_objs[[nm]] %>%
        sleuth_fit(., ~treatment) %>%
        sleuth_wt(., "treatmentIL1B") %>%
        sleuth_wt(., "treatmentBud") %>%
        sleuth_wt(., "treatmentIB")
}

```

```r pull univariate results

dea = tibble()

for (nm in names(univariate_sleuth_objs)){
    
    print(paste0("Processing ", nm, "..."))

    a = sleuth_results(univariate_sleuth_objs[[nm]], "treatmentIL1B") %>%
        mutate(treatment = "IL1B")
    b = sleuth_results(univariate_sleuth_objs[[nm]], "treatmentBud") %>%
        mutate(treatment = "Bud")
    c = sleuth_results(univariate_sleuth_objs[[nm]], "treatmentIB") %>%
        mutate(treatment = "IB")
        select(Gene = target_id, treatment, log2fold = b, FDR = qval) %>%
        mutate(group = nm)

    dea = rbind(dea, df) %>% as_tibble()

}

dea_clean = dea %>%
    separate(group, into = c("so", "condition", "time"), sep = "_") %>%
    select(Gene, condition, treatment, time, log2fold, FDR) %>%
    mutate(
        time = gsub("h", "", time),
        time = factor(time, levels = c(1, 2, 6)),
        treatment = factor(treatment, levels = c("IL1B", "Bud", "IB")),
        condition = factor(condition, levels = c("naive", "Ad-GFP", "Ad-IKBA"))
        ) %>%
    arrange(Gene, time, condition, treatment) %>%
    filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
    group_by(Gene) %>%
    filter(!all(is.na(log2fold))) %>%
    ungroup() %>%
    mutate(
        log2fold = if_else(is.na(log2fold), 0, log2fold),
        FDR = if_else(is.na(FDR), 1, FDR)
        )

write_tsv(dea_clean, "sleuth/a549-ikbadn-ib-timecourse_DEA.txt")


```

```r saving and retrieving

# use as needed

saveRDS(so, "sleuth/objects/so.rds")
saveRDS(so_1h, "sleuth/objects/so_1h.rds")
saveRDS(so_2h, "sleuth/objects/so_2h.rds")
saveRDS(so_6h, "sleuth/objects/so_6h.rds")
saveRDS(univariate_sleuth_objs, "sleuth/objects/univariate_sleuth_objs.rds")

so = readRDS("sleuth/objects/so.rds")
so_1h = readRDS("sleuth/objects/so_1h.rds")
so_2h = readRDS("sleuth/objects/so_2h.rds")
so_6h = readRDS("sleuth/objects/so_6h.rds")
univariate_sleuth_objs = readRDS("sleuth/objects/univariate_sleuth_objs.rds")
```


