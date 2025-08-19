library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(furrr)

plan(multisession, workers = 4)

slim_dir <- "SLiM_update/outputs"

calculate_freq <- function(gt_vector){
  alleles <- unlist(strsplit(gt_vector, "[/:|]"))
  sum(alleles == "1") / length(alleles)
}

process_vcf_file <- function(vcf_file){
  filename   <- basename(vcf_file)
  chrom      <- basename(dirname(vcf_file))
  generation <- as.numeric(str_match(filename, "_([0-9]+)_generations\\.vcf")[,2])
  replicate  <- str_match(filename, "generations\\.(\\d+)\\.vcf$")[,2]
  if(is.na(replicate)) replicate <- "0"
  lines <- readLines(vcf_file, n = 200)
  header_line <- grep("^#CHROM", lines)
  if(!length(header_line)) return(NULL)
  col_names <- strsplit(lines[header_line], "\t")[[1]]
  vcf_data  <- tryCatch(read.table(vcf_file, skip = header_line, header = FALSE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
  if(is.null(vcf_data) || nrow(vcf_data) == 0) return(NULL)
  colnames(vcf_data) <- col_names
  gt_matrix <- vcf_data[, 10:ncol(vcf_data)]
  freqs     <- apply(gt_matrix, 1, calculate_freq)
  mt_type <- case_when(
    str_detect(vcf_data$INFO, "MT=3") ~ "MT=3",
    str_detect(vcf_data$INFO, "MT=2") ~ "MT=2",
    str_detect(vcf_data$INFO, "MT=1") ~ "MT=1",
    TRUE ~ "MT=unknown"
  )
  tibble(CHROM = chrom, POS = vcf_data$POS, Generation = generation, Replicate = replicate, Frequency = freqs, MT = mt_type)
}

chrom_dirs <- list.dirs(slim_dir, recursive = FALSE)

for(chrom_dir in chrom_dirs){
  chrom     <- basename(chrom_dir)
  vcf_files <- list.files(chrom_dir, pattern = "_[0-9]+_generations(\\.[0-9]+)?\\.vcf$", full.names = TRUE)
  message("Processing ", chrom, " (", length(vcf_files), " VCF)…")
  chrom_freq_data <- future_map_dfr(vcf_files, process_vcf_file, .progress = TRUE)
  if(nrow(chrom_freq_data) == 0){
    message("No data for ", chrom)
    next
  }
  write_csv(chrom_freq_data, file.path(chrom_dir, paste0(chrom, "_frequency_over_generations.csv")))
  p <- ggplot(chrom_freq_data, aes(x = Generation, y = Frequency, group = interaction(POS, Replicate), color = MT)) +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    scale_color_manual(values = c("MT=1"="black","MT=2"="red","MT=3"="blue","MT=unknown"="grey")) +
    labs(title = paste("Allele-frequency trajectories:", chrom), x = "Generation", y = "Mutant-allele frequency") +
    theme_minimal(base_size = 12)
  ggsave(file.path(chrom_dir, paste0(chrom, "_frequency_plot.pdf")), p, width = 8, height = 5)
}
