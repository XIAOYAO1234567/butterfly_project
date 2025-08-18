library(dplyr)
library(readr)
library(GenomicRanges)

theta_file <- "theta.10kb.pestPG"
coverage_dir <- path.expand("~/Desktop/coverage_output_chrwise_10kb")
outdir <- "angsd_output_10kb"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cov_q <- 0.05
by_chrom <- TRUE
snp_count_file <- NULL
snp_min <- 20

theta_df <- read_tsv(theta_file, comment = "#", col_names = FALSE)
colnames(theta_df)[1:14] <- c("Region","CHROM","WinCenter","tW","tP","tF","tH","tL","Tajima","tajima2","tajima3","misc1","misc2","nSites")

theta_df <- theta_df %>%
  mutate(
    WinCenter = as.numeric(WinCenter),
    tP = as.numeric(tP),
    nSites = as.numeric(nSites),
    Tajima = as.numeric(Tajima),
    pi = tP / nSites,
    tajimaD = Tajima,
    BIN_START = WinCenter - 5000 + 1,
    BIN_END = WinCenter + 5000,
    PI = pi
  ) %>%
  filter(!is.na(PI), !is.na(tajimaD), !is.na(nSites)) %>%
  mutate(REGION = paste0(CHROM, ":", BIN_START, "-", BIN_END))

coverage_files <- list.files(coverage_dir, pattern = "\\.bed$", full.names = TRUE, recursive = TRUE)
stopifnot(length(coverage_files) > 0)

coverage_all <- coverage_files %>%
  lapply(read_tsv, col_names = c("CHROM","START","END","COV"), show_col_types = FALSE) %>%
  bind_rows() %>%
  mutate(BIN_START = as.numeric(START) + 1, BIN_END = as.numeric(END)) %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(COV = mean(as.numeric(COV), na.rm = TRUE), .groups = "drop")

theta_cov <- theta_df %>% left_join(coverage_all, by = c("CHROM","BIN_START","BIN_END"))

pi_low  <- quantile(theta_cov$PI,      0.05, na.rm = TRUE)
pi_high <- quantile(theta_cov$PI,      0.95, na.rm = TRUE)
td_low  <- quantile(theta_cov$tajimaD, 0.05, na.rm = TRUE)
td_high <- quantile(theta_cov$tajimaD, 0.95, na.rm = TRUE)

theta_cov <- theta_cov %>%
  mutate(
    flag_pi = case_when(PI > pi_high ~ "pi_high", PI < pi_low ~ "pi_low", TRUE ~ "none"),
    flag_td = case_when(tajimaD > td_high ~ "td_high", tajimaD < td_low ~ "td_low", TRUE ~ "none"),
    overlap_flag = case_when(
      flag_pi == "pi_high" & flag_td == "td_high" ~ "both_high",
      flag_pi == "pi_low"  & flag_td == "td_low"  ~ "both_low",
      TRUE ~ "none"
    )
  )

if (by_chrom) {
  theta_cov <- theta_cov %>%
    group_by(CHROM) %>%
    mutate(
      cov_low_thr  = quantile(COV, cov_q, na.rm = TRUE),
      cov_high_thr = quantile(COV, 1 - cov_q, na.rm = TRUE),
      COV_Z        = (COV - mean(COV, na.rm = TRUE)) / sd(COV, na.rm = TRUE),
      LOW_COV      = if_else(!is.na(COV) & COV < cov_low_thr, TRUE, FALSE)
    ) %>%
    ungroup()
} else {
  cov_low_thr_global  <- quantile(theta_cov$COV, cov_q, na.rm = TRUE)
  cov_high_thr_global <- quantile(theta_cov$COV, 1 - cov_q, na.rm = TRUE)
  mu_cov <- mean(theta_cov$COV, na.rm = TRUE)
  sd_cov <- sd(theta_cov$COV, na.rm = TRUE)
  theta_cov <- theta_cov %>%
    mutate(
      cov_low_thr  = cov_low_thr_global,
      cov_high_thr = cov_high_thr_global,
      COV_Z        = (COV - mu_cov) / sd_cov,
      LOW_COV      = if_else(!is.na(COV) & COV < cov_low_thr_global, TRUE, FALSE)
    )
}

if (!is.null(snp_count_file)) {
  snp_df <- read_tsv(snp_count_file, show_col_types = FALSE)
  names_lower <- tolower(names(snp_df))
  if (all(c("chrom","start","end","snp_count") %in% names_lower)) {
    snp_df <- snp_df %>% rename_with(~tolower(.x)) %>%
      mutate(BIN_START = as.numeric(start) + 1, BIN_END = as.numeric(end)) %>%
      select(CHROM = chrom, BIN_START, BIN_END, SNP_COUNT = snp_count)
  } else if (all(c("chrom","bin_start","bin_end","snp_count") %in% names_lower)) {
    snp_df <- snp_df %>% rename_with(~tolower(.x)) %>%
      select(CHROM = chrom, BIN_START = bin_start, BIN_END = bin_end, SNP_COUNT = snp_count)
  } else {
    stop("Unrecognized SNP file columns.")
  }
  theta_cov <- theta_cov %>%
    left_join(snp_df, by = c("CHROM","BIN_START","BIN_END")) %>%
    mutate(LOW_SNP = if_else(!is.na(SNP_COUNT) & SNP_COUNT < snp_min, TRUE, FALSE))
} else {
  theta_cov <- theta_cov %>% mutate(SNP_COUNT = NA_integer_, LOW_SNP = NA)
}

theta_cov <- theta_cov %>%
  mutate(label = case_when(
    overlap_flag == "both_low"                         ~ "sweep",
    flag_td == "td_high" | overlap_flag == "both_high" ~ "balancing",
    TRUE                                               ~ "neutral"
  ))

window_table <- theta_cov %>%
  transmute(
    CHROM, BIN_START, BIN_END, REGION,
    COV, cov_low_thr, cov_high_thr, COV_Z, LOW_COV,
    SNP_COUNT, LOW_SNP,
    PI, tajimaD, nSites,
    flag_pi, flag_td, overlap_flag,
    label
  )

out_windows <- file.path(outdir, "window_flags_10kb_labeled.tsv")
write_tsv(window_table, out_windows)

blk_windows <- window_table %>% filter(LOW_COV %in% TRUE | LOW_SNP %in% TRUE)

if (nrow(blk_windows) > 0) {
  gr_blk <- with(blk_windows, GRanges(seqnames = CHROM, ranges = IRanges(start = BIN_START, end = BIN_END)))
  gr_blk_merged <- reduce(gr_blk)
  blacklist_df <- data.frame(
    CHROM = as.character(seqnames(gr_blk_merged)),
    START = start(gr_blk_merged),
    END   = end(gr_blk_merged)
  ) %>%
    rowwise() %>%
    mutate(
      N_WINDOWS = sum(window_table$CHROM == CHROM & window_table$BIN_START >= START & window_table$BIN_END <= END),
      N_LOW_COV = sum(window_table$CHROM == CHROM & window_table$BIN_START >= START & window_table$BIN_END <= END & window_table$LOW_COV %in% TRUE, na.rm = TRUE),
      N_LOW_SNP = sum(window_table$CHROM == CHROM & window_table$BIN_START >= START & window_table$BIN_END <= END & window_table$LOW_SNP %in% TRUE, na.rm = TRUE),
      REASON = case_when(
        N_LOW_COV > 0 & N_LOW_SNP > 0 ~ "low_coverage+low_snp",
        N_LOW_COV > 0                 ~ "low_coverage",
        N_LOW_SNP > 0                 ~ "low_snp",
        TRUE                          ~ "unknown"
      ),
      SPAN_KB = round((END - START + 1) / 1000, 1)
    ) %>%
    ungroup() %>%
    arrange(CHROM, START)
  out_blacklist <- file.path(outdir, "blacklist_regions_10kb.tsv")
  write_tsv(blacklist_df, out_blacklist)
} else {
  out_blacklist <- file.path(outdir, "blacklist_regions_10kb.tsv")
  write_tsv(tibble(CHROM=character(), START=integer(), END=integer(), N_WINDOWS=integer(), N_LOW_COV=integer(), N_LOW_SNP=integer(), REASON=character(), SPAN_KB=double()), out_blacklist)
}

cat(
  "Done.\n",
  "Table #1 (windows + labels): ", out_windows, "\n",
  "Table #2 (blacklist):        ", out_blacklist, "\n"
)
