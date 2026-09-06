library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

if(!file.exists("../data/peaks/dba.rds")){
  
  # ===== raw data =====
  
  dat_raw <- data.frame()
  
  for (f in list.files("../data/diffbind")) {
    
    if (str_detect(f, "_EDGER")){next}
    
    info <- str_split(f, "_")[[1]]
    targ <- info[3]
    cont <- info[4]
    
    if (cont %in% c("DC", "DB")){next}
    
    df <- read_tsv(paste0("../data/diffbind/", f))
    
    missing <- setdiff(c("Conc_A","Conc_B","Conc_C","Conc_D"), names(df))
    
    df_clean <- df %>%
      mutate(!!!setNames(rep(list(NA_real_), length(missing)), missing)) %>%
      relocate(Conc_A, Conc_B, Conc_C, Conc_D) %>%
      mutate(target = targ, contrast = cont)
    
    dat_raw <- rbind(dat_raw, df_clean)
    
  }
  
  # ===== cleanup =====
  
  dba <- dat_raw %>%
    select(
      chrom = seqnames, start, end, width,
      target, contrast,
      log2conc_NS = Conc_A,
      log2conc_IL1B = Conc_B,
      log2conc_Bud = Conc_C,
      log2conc_IB = Conc_D,
      log2fold_diffbind = Fold,
      pval = p.value, FDR
    ) %>%
    mutate(treatment = case_when(
      contrast == "BA" ~ "IL1B",
      contrast == "CA" ~ "Bud",
      contrast == "DA" ~ "IB"
    ),
    .before = log2conc_NS) %>%
    
    # manually calculate log2fold from log2conc
    mutate(log2fold_conc = case_when(
      treatment == "IL1B" ~ log2conc_IL1B - log2conc_NS,
      treatment == "Bud"  ~ log2conc_Bud  - log2conc_NS,
      treatment == "IB"   ~ log2conc_IB   - log2conc_NS,
    ), .before = log2fold_diffbind) %>%
    
    # assign unique id's
    mutate(browserQuery = paste0(chrom, ":", start, "-", end), .before = width) %>%
    arrange(browserQuery) %>%
    group_by(target) %>%
    mutate(
      peak_id = paste0(
        target,
        "peak", 
        match(browserQuery, unique(browserQuery))), 
      .before = 1
    ) %>%
    ungroup() %>%
    
    # identify max conc
    group_by(peak_id, treatment) %>%
    mutate(maxConc = max(
      unique(na.omit(log2conc_NS)),
      unique(na.omit(log2conc_IL1B)),
      unique(na.omit(log2conc_Bud)),
      unique(na.omit(log2conc_IB))
    ), .before = log2conc_NS) %>%
    ungroup() %>%
    
    # assign peak significance
    mutate(sig = case_when(
      log2fold_conc >= 1 & FDR <= 0.05 & maxConc >= 4 ~ "up",
      log2fold_conc <=-1 & FDR <= 0.05 & maxConc >= 4 ~ "dn",
      TRUE ~ "ns"
    )) %>%
    
    # calculate il1b+bud differential
    group_by(peak_id, target) %>%
    mutate(IB_diff = log2fold_conc[treatment == "IB"] - log2fold_conc) %>%
    ungroup() %>%
    
    # assign trends
    mutate(trend = case_when(
      IB_diff <= -0.5 ~ "decreasing",
      IB_diff >= 0.5  ~ "increasing",
      TRUE ~ "unchanged"
    ))
  
  # simplify trends: RELA peaks take IB vs IL1B trend, GR takes IB vs Bud
  dba_trends <- dba %>%
    filter(
      treatment == "IL1B" & target == "RELA" |
        treatment == "Bud" & target == "GR"
    ) %>%
    select(peak_id, chrom, start, end, 
           browserQuery, target, treatment, sig, trend)
  
  # strip out old trend column and merge in relevant trends
  dba <- dba %>%
    select(-trend) %>%
    left_join(
      dba_trends,
      by = join_by(
        peak_id, chrom, start, end, browserQuery, target, treatment, sig
      )
    ) %>%
    group_by(peak_id, target) %>%
    mutate(trend = unique(na.omit(trend))) %>%
    ungroup()
  
  
  # ===== establish genomic ranges =====
  
  gene_ranges = genes(
    EnsDb.Hsapiens.v86,
    filter = GeneBiotypeFilter("protein_coding")
  )
  
  peaks_df = dba %>%
    select(peak_id, chrom, start, end) %>%
    mutate(chrom = gsub("chr", "", chrom)) %>%
    distinct()
  
  peaks_ranges = makeGRangesFromDataFrame(
    peaks_df,
    seqnames.field     = "chrom",
    start.field        = "start",
    end.field          = "end",
    keep.extra.columns = TRUE
  )
  
  # ===== find genes roughly within 30 kb =====
  
  hits = findOverlaps(
    peaks_ranges,
    gene_ranges,
    maxgap = 30000,
    ignore.strand = TRUE
  )
  
  # ===== calculate distance from TSS to peak center =====
  tss_ranges = promoters(
    gene_ranges,
    upstream = 0,
    downstream = 1
  )
  
  peak_centers = resize(peaks_ranges, width = 1, fix = "center")
  
  distances = distance(
    peak_centers[queryHits(hits)],
    tss_ranges[subjectHits(hits)]
  )
  
  # ===== clean and refine =====
  
  peak_gene_map = tibble(
    peak_id = peaks_ranges$peak_id[queryHits(hits)],
    Gene = gene_ranges$gene_name[subjectHits(hits)],
    distance_to_TSS = distances
  ) %>%
    filter(grepl("^[A-Z0-9]+$", Gene)) %>%
    filter(distance_to_TSS <= 30000)
  
  peak_gene_summary = peak_gene_map %>%
    group_by(peak_id) %>%
    summarize(genes = paste(unique(Gene), collapse = ","))
  
  peaks_to_genes = list(
    "long_data" = peak_gene_map,
    "gene_summary" = peak_gene_summary
  )
  
  
  # ===== write files =====
  
  saveRDS(dba, "../data/peaks/dba.rds")
  saveRDS(peaks_to_genes, "../data/peaks/peaks_to_genes.rds")
  
  cleanup_objects()
  
} else{
  
  dba = readRDS("../data/peaks/dba.rds")
  peaks_to_genes = readRDS("../data/peaks/peaks_to_genes.rds")
  
}

