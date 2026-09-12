library(dorothea)
library(tidyverse)

regulons <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C")) %>%
  mutate(likelihood = case_when(
    confidence == "A" ~ 1.00,
    confidence == "B" ~ 0.75,
    confidence == "C" ~ 0.50
  ))

regulons_summary <- regulons %>%
  group_by(tf) %>%
  summarize(
    max_confidence = if_else("A" %in% confidence, "A", if_else("B" %in% confidence, "B", "C")),
    n_targets      = length(unique(target)),
    target_list    = list(target)
  )

rg <- regulons %>%
  split(.$tf) %>%
  map(~list(
    tfmode     = setNames(.x$mor, .x$target),
    likelihood = setNames(.x$likelihood, .x$target)
  ))

dorothea_regulons = list(
  "regulons" = regulons,
  "summary" = regulons_summary,
  "viper_input" = rg
)

rm(regulons, regulons_summary, rg)