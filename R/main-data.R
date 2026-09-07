tc_dea <- read_tsv("../data/raw/A549_IL1B_Bud_tc_univariateDEA.txt") %>%
  mutate(sig = case_when(
    log2fold <= -1 & FDR <= 0.05 ~ "dn",
    log2fold >=  1 & FDR <= 0.05 ~ "up",
    TRUE ~ "ns"
  )) %>%
  group_by(Gene, treatment) %>%
  mutate(
    peak = time[which.max(abs(log2fold))],
    peak_sig = sig[which(time == peak)]) %>%
  ungroup() %>%
  mutate(
    time = factor(time, levels = c(1, 2, 6, 12, 24)),
    peak = factor(peak, levels = c(1, 2, 6, 12, 24)),
    sig = factor(sig, levels = c("up", "dn", "ns")),
    peak_sig = factor(peak_sig, levels = c("up", "dn", "ns")),
    treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "I+B"))
  )

tc_tpm <- read_tsv("../data/raw/A549_IL1B_Bud_tc_tpm.txt") %>%
  mutate(
    time = factor(time, levels = c(1, 2, 6, 12, 24)),
    treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "I+B"))
  )

ikba_dea <- read_tsv("../data/raw/a549-ikbadn-ib-timecourse_DEA.txt") %>%
  mutate(sig = case_when(
    log2fold <= -1 & FDR <= 0.05 ~ "dn",
    log2fold >=  1 & FDR <= 0.05 ~ "up",
    TRUE ~ "ns"
  )) %>%
  group_by(Gene, treatment, condition) %>%
  mutate(
    peak = time[which.max(abs(log2fold))],
    peak_sig = sig[which(time == peak)]) %>%
  ungroup() %>%
  mutate(
    treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "IB")),
    condition = factor(condition, levels = c("naive", "Ad-GFP", "Ad-IKBA")),
    time = factor(time, levels = c(1, 2, 6))
  )

ikba_tpm <- read_tsv("../data/raw/a549-ikbadn-ib-timecourse_tpm.txt") %>%
  mutate(
    treatment = factor(treatment, levels = c("NS", "IL1B", "Bud", "IB")),
    condition = factor(condition, levels = c("naive", "Ad-GFP", "Ad-IKBA")),
    time = factor(time, levels = c(1, 2, 6))
  )

primary_dea = read_tsv("../data/raw/A549vsPrimary_univariateDEA.txt") %>%
  mutate(
    sig = case_when(
      log2fold >= 1 & FDR <= 0.05 ~ "up",
      log2fold <=-1 & FDR <= 0.05 ~ "dn",
      TRUE ~ "ns"
    )
  ) %>%
  mutate(
    celltype = factor(celltype, levels = c("A549", "ALI", "HBE")),
    treatment = factor(treatment, levels = c("NS", "IL1B", "GC", "combo"))
  )

hbe_tc_dea = read_tsv("../data/raw/hbe_timecourse_dea.txt") %>%
  mutate(
    sig = case_when(
      log2fold >= 1 & FDR <= 0.05 ~ "up",
      log2fold <=-1 & FDR <= 0.05 ~ "dn",
      TRUE ~ "ns"
    )
  ) %>%
  mutate(
    celltype = factor(celltype, levels = c("A549", "ALI", "HBE")),
    treatment = factor(treatment, levels = c("NS", "IL1B", "GC", "combo")),
    time = factor(time, levels = c(2, 6, 24))
  )

primary_tpm = read_tsv("../data/raw/A549vsPrimary_tpm.txt") %>%
  mutate(
    celltype = factor(celltype, levels = c("A549", "ALI", "HBE")),
    treatment = factor(treatment, levels = c("NS", "IL1B", "GC", "combo"))
  )