ssh alex.gao1@arc.ucalgary.ca

## ARC SETUP
```bash SETUP TOOLS
salloc -c 4 --mem 64GB --time 05:00:00

MULTIQC="/home/alex.gao1/tools/multiqc_latest.sif"
apptainer exec --bind /work:/work "$MULTIQC" multiqc --version  #v1.19

module load kallisto
kallisto --version  #0.46.1

module load R   #4.4.1

mkdir -p /work/newton_lab/ag_analysis/A549_IL1B_Bud_tc
cd /work/newton_lab/ag_analysis/A549_IL1B_Bud_tc
```

## INHERIT FILES FROM PREVIOUS ANALYSIS
```bash
# symlink to MDSC519 project
ln -s /work/newton_lab/ag_analysis/MDSC519_project/kallisto/tc kallisto

# copy in meta file
cp /work/newton_lab/ag_analysis/A549_vs_primary/meta_a549_il1b_bud.txt .
```

## SLEUTH - univariate analysis by timepoint
```bash
mkdir -p ./sleuth
mkdir -p ./sleuth/objects
cd ./sleuth
R
```

```R PREP SLEUTH INPUTS
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/A549_IL1B_Bud_tc/sleuth")

library(dplyr)
library(readr)
library(sleuth)

# transcript ids to gene names
t2g <- read_tsv("../../ref_seq/Homo_sapiens.GRCh38.p14.cdna.all_mart_export.txt") %>%
  select(target_id = `Transcript stable ID`, Gene = `Gene name`) %>%
  na.omit() %>% 
  distinct()

#filter requiring at least 5 reads in at least 20% of samples
new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

# meta files
s2c <- read_tsv("../meta_a549_il1b_bud.txt") %>%
  select(sample, rep, treatment, time, path = kallisto_path) %>%
  mutate(treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "I+B")))

s2c_1h  = filter(s2c, time == 1)
s2c_2h  = filter(s2c, time == 2)
s2c_6h  = filter(s2c, time == 6)
s2c_12h = filter(s2c, time == 12)
s2c_24h = filter(s2c, time == 24)
```

```R CREATE SLEUTH OBJECTS
# function to process timecourse DEA groups
run_sleuth <- function(myS2C){
  # build sleuth object
  out <- sleuth_prep(sample_to_covariates = myS2C,
                     full_model = ~treatment,
                     target_mapping = t2g,
                     gene_mode = TRUE, 
                     aggregation_column = "Gene",
                     filter_fun = new_filter,
                     num_cores = 4)
  # fit model
  out <- sleuth_fit(out, ~treatment)
  # run wald tests
  out <- sleuth_wt(out, "treatmentIL1B")
  out <- sleuth_wt(out, "treatmentBud")
  out <- sleuth_wt(out, "treatmentI+B")
  gc()
  out
}

# build sleuth objects
so_1h  <- run_sleuth(s2c_1h)
saveRDS(so_1h, "objects/so_1h.rds")
so_2h  <- run_sleuth(s2c_2h)
saveRDS(so_2h, "objects/so_2h.rds")
so_6h  <- run_sleuth(s2c_6h)
saveRDS(so_6h, "objects/so_6h.rds")
so_12h <- run_sleuth(s2c_12h)
saveRDS(so_12h, "objects/so_12h.rds")
so_24h <- run_sleuth(s2c_24h)
saveRDS(so_24h, "objects/so_24h.rds")

# so_1h  = readRDS("objects/so_1h.rds")
# so_2h  = readRDS("objects/so_2h.rds")
# so_6h  = readRDS("objects/so_6h.rds")
# so_12h = readRDS("objects/so_12h.rds")
# so_24h = readRDS("objects/so_24h.rds")

so_all <- sleuth_prep(sample_to_covariates = s2c,
                      full_model = ~treatment,
                      target_mapping = t2g,
                      gene_mode = TRUE,
                      aggregation_column = "Gene",
                      filter_fun = new_filter,
                      num_cores = 4)
saveRDS(so_all, "objects/so_all.rds")

# so_all = readRDS("objects/so_all.rds")
```

```R ASSEMBLE AND CLEAN DATA
#function to pull results
pull_sleuth_results <- function(obj, treatment, time){
  
  out <- sleuth_results(obj, paste0("treatment", treatment), "wt") %>%
    select(Gene = target_id, log2fold = b, FDR = qval) %>%
    mutate(log2fold = log2fold / log(2)) %>%
    mutate(treatment = treatment, time = time, .before = 2)
  out
}

#assemble results
a <- pull_sleuth_results(so_1h, "IL1B", 1)
b <- pull_sleuth_results(so_1h, "Bud", 1)
c <- pull_sleuth_results(so_1h, "I+B", 1)

d <- pull_sleuth_results(so_2h, "IL1B", 2)
e <- pull_sleuth_results(so_2h, "Bud", 2)
f <- pull_sleuth_results(so_2h, "I+B", 2)

g <- pull_sleuth_results(so_6h, "IL1B", 6)
h <- pull_sleuth_results(so_6h, "Bud", 6)
i <- pull_sleuth_results(so_6h, "I+B", 6)

j <- pull_sleuth_results(so_12h, "IL1B", 12)
k <- pull_sleuth_results(so_12h, "Bud", 12)
l <- pull_sleuth_results(so_12h, "I+B", 12)

m <- pull_sleuth_results(so_24h, "IL1B", 24)
n <- pull_sleuth_results(so_24h, "Bud", 24)
o <- pull_sleuth_results(so_24h, "I+B", 24)

dea_results <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) %>%
  mutate(
    treatment = factor(treatment, levels = c("IL1B", "Bud", "I+B")),
    time      = factor(time,      levels = c(1, 2, 6, 12, 24))
    ) %>%
  arrange(Gene, treatment, time) %>%
  filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
  group_by(Gene) %>%
  filter(!all(is.na(log2fold))) %>%
  ungroup() %>%
  mutate(
    log2fold = if_else(is.na(log2fold), 0, log2fold),
    FDR = if_else(is.na(FDR), 1, FDR)
    )

rm(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)

write_tsv(dea_results, "A549_IL1B_Bud_tc_univariateDEA.txt")

#assemble tc tpm table
tpm <- kallisto_table(so_all, use_filtered = FALSE) %>%
  select(Gene = target_id, rep, treatment, time, tpm) %>%
  mutate(
    treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "I+B")),
    time      = factor(time, levels = c(1, 2, 6, 12, 24))
  ) %>%
  filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
  group_by(Gene, rep, time) %>%
  filter(any(treatment == "NS")) %>%       #R3_12h_I+B is missing
  mutate(log2tpm = log2(tpm + 0.1)) %>%
  mutate(log2fold = log2tpm - log2tpm[treatment == "NS"]) %>%
  ungroup() %>%
  mutate(fold = 2^log2fold) %>%
  select(Gene, rep, treatment, time, tpm, log2tpm, fold, log2fold)

write_tsv(tpm, "A549_IL1B_Bud_tc_tpm.txt")
```

## SLEUTH - IL1B+Bud vs IL1B
```bash
mkdir -p ./sleuth
mkdir -p ./sleuth/objects
cd ./sleuth
R
```

```R
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/A549_IL1B_Bud_tc/sleuth")

library(dplyr)
library(readr)
library(sleuth)

# t2g mapping
t2g <- read_tsv("../../ref_seq/Homo_sapiens.GRCh38.p14.cdna.all_mart_export.txt") %>%
  select(target_id = `Transcript stable ID`, Gene = `Gene name`) %>%
  na.omit() %>% 
  distinct()

# metafile, only IL1B and IB samples
s2c <- read_tsv("../meta_a549_il1b_bud.txt") %>%
  select(sample, rep, treatment, time, path = kallisto_path) %>%
  filter(treatment %in% c("IL1B", "I+B")) %>%
  mutate(treatment = factor(treatment, levels = c("IL1B", "I+B")))

# filter requiring at least 5 reads in at least 20% of samples
new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

# make DEA groups
s2c_1h <- filter(s2c, time == 1)
s2c_2h <- filter(s2c, time == 2)
s2c_6h <- filter(s2c, time == 6)
s2c_12h <- filter(s2c, time == 12)
s2c_24h <- filter(s2c, time == 24)

# create sleuth objects
so_1h = sleuth_prep(sample_to_covariates = s2c_1h,
                    full_model = ~treatment,
                    target_mapping = t2g,
                    gene_mode = TRUE, aggregation_column = "Gene",
                    filter_fun = new_filter, num_cores = 4)

so_2h = sleuth_prep(sample_to_covariates = s2c_2h,
                    full_model = ~treatment,
                    target_mapping = t2g,
                    gene_mode = TRUE, aggregation_column = "Gene",
                    filter_fun = new_filter, num_cores = 4)

so_6h = sleuth_prep(sample_to_covariates = s2c_6h,
                    full_model = ~treatment,
                    target_mapping = t2g,
                    gene_mode = TRUE, aggregation_column = "Gene",
                    filter_fun = new_filter, num_cores = 4)

so_12h = sleuth_prep(sample_to_covariates = s2c_12h,
                    full_model = ~treatment,
                    target_mapping = t2g,
                    gene_mode = TRUE, aggregation_column = "Gene",
                    filter_fun = new_filter, num_cores = 4)

so_24h = sleuth_prep(sample_to_covariates = s2c_24h,
                    full_model = ~treatment,
                    target_mapping = t2g,
                    gene_mode = TRUE, aggregation_column = "Gene",
                    filter_fun = new_filter, num_cores = 4)

# run tests
so_1h = so_1h %>% sleuth_fit(., ~treatment) %>% sleuth_wt(., "treatmentI+B")
so_2h = so_2h %>% sleuth_fit(., ~treatment) %>% sleuth_wt(., "treatmentI+B")
so_6h = so_6h %>% sleuth_fit(., ~treatment) %>% sleuth_wt(., "treatmentI+B")
so_12h = so_12h %>% sleuth_fit(., ~treatment) %>% sleuth_wt(., "treatmentI+B")
so_24h = so_24h %>% sleuth_fit(., ~treatment) %>% sleuth_wt(., "treatmentI+B")

# pull tables
a = sleuth_results(so_1h, "treatmentI+B", "wt") %>%
   select(Gene = target_id, fold = b, FDR = qval) %>%
   mutate(fold = fold / log(2)) %>%
   mutate(treatment = "I+B", time = 1, .before = 2)

b = sleuth_results(so_2h, "treatmentI+B", "wt") %>%
   select(Gene = target_id, fold = b, FDR = qval) %>%
   mutate(fold = fold / log(2)) %>%
   mutate(treatment = "I+B", time = 2, .before = 2)

c = sleuth_results(so_6h, "treatmentI+B", "wt") %>%
   select(Gene = target_id, fold = b, FDR = qval) %>%
   mutate(fold = fold / log(2)) %>%
   mutate(treatment = "I+B", time = 6, .before = 2)

d = sleuth_results(so_12h, "treatmentI+B", "wt") %>%
   select(Gene = target_id, fold = b, FDR = qval) %>%
   mutate(fold = fold / log(2)) %>%
   mutate(treatment = "I+B", time = 12, .before = 2)

e = sleuth_results(so_24h, "treatmentI+B", "wt") %>%
   select(Gene = target_id, fold = b, FDR = qval) %>%
   mutate(fold = fold / log(2)) %>%
   mutate(treatment = "I+B", time = 24, .before = 2)

dea = rbind(a,b,c,d,e) %>%
   group_by(Gene) %>%
   filter(!all(is.na(fold))) %>%
   filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
   ungroup() %>%
   arrange(Gene, time) %>%
   select(Gene, time, IB_diff = fold, IB_FDR = FDR) %>%
   mutate(
    IB_diff = if_else(is.na(IB_diff), 0, IB_diff),
    IB_FDR = if_else(is.na(IB_FDR), 1, IB_FDR)
    )

write_tsv(dea, "A549_IL1B_Bud_tc_IBdiff.txt")

# save sleuth objects
saveRDS(so_1h, "objects/so_IBdiff_1h.rds")
saveRDS(so_2h, "objects/so_IBdiff_2h.rds")
saveRDS(so_6h, "objects/so_IBdiff_6h.rds")
saveRDS(so_12h, "objects/so_IBdiff_12h.rds")
saveRDS(so_24h, "objects/so_IBdiff_24h.rds")

q()
n
```


## SLEUTH - transcript level analysis
```bash
R
```

```R
.libPaths("/home/alex.gao1/R")
setwd("/work/newton_lab/ag_analysis/MDSC519_project/sleuth")

library(dplyr)
library(readr)
library(sleuth)

# t2g mapping
t2g <- read_tsv("/work/newton_lab/ag_analysis/ref_seq/Homo_sapiens.GRCh38.p14.cdna.all_mart_export.txt") %>%
  select(target_id = `Transcript stable ID`, Gene = `Gene name`) %>%
  na.omit() %>% 
  distinct()

# s2c
tc_s2c <- read_tsv("tc_meta.txt") %>%
  mutate(treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "IL1B + Bud")))

# filter requiring at least 5 reads in at least 20% of samples
new_filter <- function(row, min_reads = 5, min_prop = 0.2){mean(row >= min_reads) >= min_prop}

# normalize reads
so_tc_transcripts <- sleuth_prep(sample_to_covariates = tc_s2c,
                                target_mapping = t2g,
                                filter_fun = new_filter,
                                num_cores = 4)

saveRDS(so_tc_transcripts, "tc_objects/so_tc_transcripts.rds")
so_tc_transcripts = readRDS("tc_objects/so_tc_transcripts.rds")

# assemble tc tpm table
tc_tpm_transcripts <- kallisto_table(so_tc_transcripts, use_filtered = FALSE) %>%
                      select(transcript = target_id, rep, treatment, time, tpm) %>%
                      mutate(target_id = sub("\\.\\d+$", "", transcript)) %>%
                      left_join(t2g, by = "target_id") %>%
                      filter(grepl("^[A-Za-z0-9]+$", Gene)) %>%
                      select(Gene, transcript, rep, treatment, time, tpm) %>%
                      arrange(Gene, transcript, rep, treatment, time, tpm)

write_tsv(tc_tpm_transcripts, "tc_data/tc_transcripts_tpm.txt")

q()
n
```

## EXPORTS
cd C:/Users/alexg/"OneDrive - University of Calgary"/"seqProject NFKB_GR"/"IL Bud timecourse RNAseq"/"branch2_temporal_programs"/data/raw

scp alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/A549_IL1B_Bud_tc/sleuth/"A549_IL1B_Bud_tc_tpm.txt" ./

scp alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/A549_IL1B_Bud_tc/sleuth/"A549_IL1B_Bud_tc_univariateDEA.txt" ./

scp alex.gao1@arc.ucalgary.ca:/work/newton_lab/ag_analysis/A549_IL1B_Bud_tc/sleuth/"A549_IL1B_Bud_tc_IBdiff.txt" ./