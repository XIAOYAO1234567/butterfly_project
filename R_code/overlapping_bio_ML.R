library(dplyr)
library(readr)
library(GenomicRanges)

canon_chr <- function(x){x <- as.character(x); x <- sub("^chr","",x,ignore.case=TRUE); sub("\\.\\d+$","",x)}
pick_col <- function(df, ...) {
  nms <- tolower(names(df)); cands <- tolower(c(...))
  for (nm in cands) { i <- match(nm, nms); if (!is.na(i)) return(names(df)[i]) }
  stop("Missing a required column among: ", paste(cands, collapse=", "))
}

std_regions <- function(path){
  x <- read_tsv(path, show_col_types = FALSE)
  CH <- pick_col(x, "chrom","seqid")
  ST <- pick_col(x, "start","bin_start","begin")
  ED <- pick_col(x, "end","bin_end","stop")
  x %>%
    transmute(CHROM = canon_chr(.data[[CH]]),
              START = as.integer(.data[[ST]]),
              END   = as.integer(.data[[ED]])) %>%
    filter(!is.na(CHROM), !is.na(START), !is.na(END), START <= END)
}

ml  <- std_regions("prediction_results_merged.tsv")
bio <- std_regions("Bio_regions_nonblacklist.tsv")

gr_ml  <- GRanges(ml$CHROM,  IRanges(ml$START,  ml$END))
gr_bio <- GRanges(bio$CHROM, IRanges(bio$START, bio$END))

ov <- findOverlaps(gr_ml, gr_bio)
ints <- pintersect(gr_ml[queryHits(ov)], gr_bio[subjectHits(ov)])
res <- tibble(
  CHROM = as.character(seqnames(ints)),
  START = start(ints),
  END   = end(ints)
) %>%
  filter(END >= START) %>%
  distinct() %>%
  arrange(CHROM, START, END)

write_tsv(res, "overlap_overall.tsv")
cat("Saved: overlap_overall.tsv\n")
