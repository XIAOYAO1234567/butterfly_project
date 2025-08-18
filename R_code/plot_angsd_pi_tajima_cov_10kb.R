library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(gridExtra)

theta_file <- "theta.10kb.pestPG"
coverage_dir <- path.expand("~/Desktop/coverage_output_chrwise_10kb")
outdir <- "angsd_pi_tajima_plots_10kb"
dir.create(outdir, showWarnings = FALSE)

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
  filter(!is.na(PI), !is.na(tajimaD), !is.na(nSites), nSites > 0)

coverage_files <- list.files(coverage_dir, pattern = "\\.bed$", full.names = TRUE, recursive = TRUE)

coverage_all <- coverage_files %>%
  lapply(function(f) {
    df <- read_tsv(f, col_names = FALSE, show_col_types = FALSE)
    tibble(CHROM = df[[1]], START = as.integer(df[[2]]), END = as.integer(df[[3]]), COV = suppressWarnings(as.numeric(df[[7]])))
  }) %>%
  bind_rows() %>%
  mutate(BIN_START = START + 1) %>%
  group_by(CHROM, BIN_START) %>%
  summarise(COV = mean(COV, na.rm = TRUE), .groups = "drop")

theta_cov <- theta_df %>% left_join(coverage_all, by = c("CHROM","BIN_START"))

pi_low  <- quantile(theta_cov$PI, 0.05, na.rm = TRUE)
pi_high <- quantile(theta_cov$PI, 0.95, na.rm = TRUE)
td_low  <- quantile(theta_cov$tajimaD, 0.05, na.rm = TRUE)
td_high <- quantile(theta_cov$tajimaD, 0.95, na.rm = TRUE)
cov_low <- quantile(theta_cov$COV, 0.05, na.rm = TRUE)

theta_cov <- theta_cov %>%
  mutate(
    flag_pi = case_when(PI > pi_high ~ "pi_high", PI < pi_low ~ "pi_low", TRUE ~ "none"),
    flag_td = case_when(tajimaD > td_high ~ "td_high", tajimaD < td_low ~ "td_low", TRUE ~ "none"),
    overlap_flag = case_when(flag_pi == "pi_high" & flag_td == "td_high" ~ "both_high",
                             flag_pi == "pi_low"  & flag_td == "td_low"  ~ "both_low",
                             TRUE ~ "none")
  )

chroms <- unique(theta_cov$CHROM)
pdf(file.path(outdir, "angsd_pi_tajima_plots_10kb.pdf"), width = 10, height = 12)

print(
  ggplot(theta_cov, aes(x = PI, y = tajimaD)) +
    geom_point(color = "grey60", size = 1, alpha = 0.5) +
    geom_point(data = filter(theta_cov, overlap_flag == "both_high"), color = "red", size = 2, shape = 21, fill = "gold") +
    geom_point(data = filter(theta_cov, overlap_flag == "both_low"),  color = "blue", size = 2, shape = 21, fill = "deepskyblue") +
    labs(title = "Genome-wide pi vs. Tajima's D", x = expression("pi"), y = "Tajima's D") +
    theme_minimal()
)

grid.arrange(
  ggplot(theta_cov, aes(x = PI)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = c(pi_low, pi_high), linetype = "dashed", color = c("blue","red")) +
    labs(title = "Histogram of pi", x = expression("pi"), y = "Count") +
    theme_minimal(),
  ggplot(theta_cov, aes(x = tajimaD)) +
    geom_histogram(bins = 50, fill = "tomato", color = "black") +
    geom_vline(xintercept = c(td_low, td_high), linetype = "dashed", color = c("blue","red")) +
    labs(title = "Histogram of Tajima's D", x = "Tajima's D", y = "Count") +
    theme_minimal(),
  ncol = 1
)

chr_offsets <- theta_cov %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(BIN_END, na.rm = TRUE), .groups = "drop") %>%
  mutate(offset = cumsum(dplyr::lag(chr_len, default = 0)))

theta_cov <- theta_cov %>%
  left_join(chr_offsets, by = "CHROM") %>%
  mutate(GLOBAL_POS = BIN_START + offset)

p_pi_global <- ggplot(theta_cov, aes(x = GLOBAL_POS, y = PI)) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = c(pi_low, pi_high), linetype = "dashed", color = c("blue","red")) +
  labs(title = "Genome-wide π curve (10kb)", x = "Genomic position", y = expression("π")) +
  theme_minimal()

p_td_global <- ggplot(theta_cov, aes(x = GLOBAL_POS, y = tajimaD)) +
  geom_line(color = "tomato") +
  geom_hline(yintercept = c(td_low, td_high), linetype = "dashed", color = c("blue","red")) +
  labs(title = "Genome-wide Tajima's D curve (10kb)", x = "Genomic position", y = "Tajima's D") +
  theme_minimal()

p_cov_global <- ggplot(theta_cov, aes(x = GLOBAL_POS, y = COV)) +
  geom_line(color = "darkgreen") +
  geom_hline(yintercept = cov_low, linetype = "dashed", color = "purple") +
  labs(title = "Genome-wide Coverage curve (10kb)", x = "Genomic position", y = "Coverage") +
  theme_minimal()

grid.arrange(p_pi_global, p_td_global, p_cov_global, ncol = 1)

for (chr in chroms) {
  df <- filter(theta_cov, CHROM == chr)
  p1 <- ggplot(df, aes(x = BIN_START, y = PI)) +
    geom_ribbon(aes(ymin = pi_low, ymax = pi_high), fill = "grey80", alpha = 0.4) +
    geom_line(color = "steelblue") +
    labs(title = paste("pi:", chr), x = "Position", y = expression("pi")) +
    theme_minimal()
  p2 <- ggplot(df, aes(x = BIN_START, y = tajimaD)) +
    geom_ribbon(aes(ymin = td_low, ymax = td_high), fill = "grey80", alpha = 0.4) +
    geom_line(color = "tomato") +
    labs(title = paste("Tajima's D:", chr), x = "Position", y = "Tajima's D") +
    theme_minimal()
  p3 <- ggplot(df, aes(x = BIN_START, y = COV)) +
    geom_line(color = "darkgreen") +
    geom_hline(yintercept = cov_low, linetype = "dashed", color = "purple") +
    labs(title = paste("Coverage:", chr), x = "Position", y = "Coverage") +
    theme_minimal()
  p4 <- ggplot(df, aes(x = PI, y = tajimaD)) +
    geom_point(color = "grey50", size = 1, alpha = 0.5) +
    geom_point(data = filter(df, overlap_flag == "both_high"), color = "red", size = 2, shape = 21, fill = "gold") +
    geom_point(data = filter(df, overlap_flag == "both_low"),  color = "blue", size = 2, shape = 21, fill = "deepskyblue") +
    labs(title = paste("pi vs Tajima's D:", chr), x = expression("π"), y = "Tajima's D") +
    theme_minimal()
  grid.arrange(p1, p2, p3, p4, ncol = 1)
}

dev.off()
cat("Done. Figures saved to:", file.path(outdir, "angsd_pi_tajima_plots_10kb.pdf"), "\n")
