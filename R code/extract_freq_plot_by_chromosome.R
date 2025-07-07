library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(furrr)

plan(multisession, workers = 4) 

slim_dir <- "SLiM_update/output"

calculate_freq <- function(gt_vector) {
  alleles <- unlist(strsplit(gt_vector, "[/:|]"))
  freq <- sum(alleles == "1", na.rm = TRUE) / (length(gt_vector) * 2)
  return(freq)
}

process_vcf_file <- function(vcf_file) {
  filename <- basename(vcf_file)
  chrom <- str_split(vcf_file, "/", simplify = TRUE)[,2]
  generation <- as.numeric(str_match(filename, "_([0-9]+)_generations")[,2])
  replicate <- str_match(filename, "generations\\.(\\d+)\\.vcf$")[,2]
  if (is.na(replicate)) replicate <- "0"

  lines <- readLines(vcf_file)
  header_line <- grep("^#CHROM", lines)
  if (length(header_line) == 0) return(NULL)
  col_names <- strsplit(lines[header_line], "\t")[[1]]

  vcf_data <- tryCatch({
    df <- read.table(vcf_file, skip = header_line, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(df) <- col_names
    df
  }, error = function(e) return(NULL))

  if (is.null(vcf_data) || nrow(vcf_data) == 0) return(NULL)

  gt_matrix <- vcf_data[, 10:ncol(vcf_data)]
  freqs <- apply(gt_matrix, 1, calculate_freq)
  info_field <- vcf_data$INFO
  mt_type <- case_when(
    str_detect(info_field, "MT=3") ~ "MT=3",
    str_detect(info_field, "MT=2") ~ "MT=2",
    str_detect(info_field, "MT=1") ~ "MT=1",
    TRUE ~ "MT=unknown"
  )

  result <- data.frame(
    CHROM = chrom,
    POS = vcf_data$POS,
    Generation = generation,
    Replicate = replicate,
    Frequency = freqs,
    MT = mt_type
  )
  return(result)
}

chrom_dirs <- list.dirs(slim_dir, recursive = FALSE)

for (chrom_dir in chrom_dirs) {
  chrom <- basename(chrom_dir)
  vcf_files <- list.files(chrom_dir, pattern = "_[0-9]+_generations(\\.[0-9]+)?\\.vcf$", full.names = TRUE)
  cat("Processing", chrom, "with", length(vcf_files), "VCF files...\n")

  # 并行处理每个文件
  chrom_freq_data <- future_map_dfr(vcf_files, process_vcf_file, .progress = TRUE)

  if (nrow(chrom_freq_data) == 0) {
    cat("No data found for", chrom, "\n")
    next
  }

  # 保存结果
  write_csv(chrom_freq_data, file.path(chrom_dir, paste0(chrom, "_frequency_over_generations.csv")))

  # 绘图
  p <- ggplot(chrom_freq_data, aes(x = Generation, y = Frequency,
                                   group = interaction(POS, Replicate),
                                   color = MT)) +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    scale_color_manual(values = c("MT=1" = "black", "MT=2" = "red", "MT=3" = "blue", "MT=unknown" = "gray")) +
    labs(title = paste0("Allele Frequency Trajectories - ", chrom),
         subtitle = "Each line = one variant x replicate",
         x = "Generation", y = "Mutant Allele Frequency",
         color = "Mutation Type") +
    theme_minimal(base_size = 12)

  ggsave(file.path(chrom_dir, paste0(chrom, "_frequency_plot.pdf")), p, width = 8, height = 5)
}
