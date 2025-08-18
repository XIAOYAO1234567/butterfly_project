library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)

in_dir  <- "angsd_output_10kb"
bio_file <- file.path(in_dir, "Bio_regions_nonblacklist.tsv")
ml_file  <- file.path(in_dir, "prediction_results_merged.tsv")
out_counts <- file.path(in_dir, "venn_region_cluster_counts_colored.tsv")
out_pdf    <- file.path(in_dir, "venn_regions_colored.pdf")
out_png_sw <- file.path(in_dir, "venn_regions_colored_sweep.png")
out_png_bl <- file.path(in_dir, "venn_regions_colored_balancing.png")

canon_chrom <- function(x){x <- as.character(x); x <- sub("^chr","",x,ignore.case=TRUE); sub("\\.\\d+$","",x)}
canon_label <- function(x){y <- tolower(trimws(as.character(x))); case_when(grepl("balanc",y)~"balancing", grepl("sweep",y)~"sweep", grepl("neut",y)~"neutral", TRUE~y)}
circle_df <- function(cx,cy,r,n=360){t <- seq(0,2*pi,length.out=n); data.frame(x=cx+r*cos(t), y=cy+r*sin(t))}

cluster_counts <- function(df_ml, df_bio, label){
  M <- df_ml %>% filter(LABEL==label) %>% select(CHROM,START,END) %>% drop_na()
  B <- df_bio %>% filter(LABEL==label) %>% select(CHROM,START,END) %>% drop_na()
  M$START <- as.integer(M$START); M$END <- as.integer(M$END)
  B$START <- as.integer(B$START); B$END <- as.integer(B$END)
  if(nrow(M)+nrow(B)==0) return(tibble(LABEL=label, ml_total=0, bio_total=0, ml_only=0, bio_only=0, both=0))
  ml_total <- nrow(M); bio_total <- nrow(B)
  chroms <- union(unique(M$CHROM), unique(B$CHROM))
  tally <- c(ml_only=0, bio_only=0, both=0)
  for(chr in chroms){
    A <- M %>% filter(CHROM==chr) %>% select(START,END) %>% as.matrix()
    C <- B %>% filter(CHROM==chr) %>% select(START,END) %>% as.matrix()
    if(nrow(A)==0 && nrow(C)==0) next
    combined <- bind_rows(if(nrow(A)) tibble(s=A[,1],e=A[,2],ml=1L,bio=0L) else tibble(s=integer(),e=integer(),ml=integer(),bio=integer()),
                          if(nrow(C)) tibble(s=C[,1],e=C[,2],ml=0L,bio=1L) else tibble(s=integer(),e=integer(),ml=integer(),bio=integer())) %>% arrange(s,e)
    cur_s <- NA_integer_; cur_e <- NA_integer_; has_ml <- 0L; has_bio <- 0L
    flush_cluster <- function(){ if(is.na(cur_s)) return(invisible(NULL)); if(has_ml==1L && has_bio==1L) tally["both"] <<- tally["both"]+1L else if(has_ml==1L) tally["ml_only"] <<- tally["ml_only"]+1L else if(has_bio==1L) tally["bio_only"] <<- tally["bio_only"]+1L }
    if(nrow(combined)>0){
      cur_s <- combined$s[1]; cur_e <- combined$e[1]; has_ml <- combined$ml[1]; has_bio <- combined$bio[1]
      if(nrow(combined)>1){
        for(i in 2:nrow(combined)){
          s <- combined$s[i]; e <- combined$e[i]; m <- combined$ml[i]; b <- combined$bio[i]
          if(s <= cur_e + 1L){ if(e>cur_e) cur_e <- e; has_ml <- as.integer(has_ml | m); has_bio <- as.integer(has_bio | b) }
          else { flush_cluster(); cur_s <- s; cur_e <- e; has_ml <- m; has_bio <- b }
        }
      }
      flush_cluster()
    }
  }
  tibble(LABEL=label, ml_total=ml_total, bio_total=bio_total, ml_only=as.integer(tally["ml_only"]), bio_only=as.integer(tally["bio_only"]), both=as.integer(tally["both"]))
}

draw_two_set <- function(left_label,right_label,left_total,right_total,left_only,right_only,both,
                         fill_left="#D95F02",fill_right="#1B9E77",edge="#333333",title=NULL){
  sep <- 1.3; r <- 1.2; L <- circle_df(-sep/2,0,r); R <- circle_df(sep/2,0,r)
  ggplot() +
    geom_polygon(data=L,aes(x,y),fill=scales::alpha(fill_left,0.35),color=edge,linewidth=0.6) +
    geom_polygon(data=R,aes(x,y),fill=scales::alpha(fill_right,0.35),color=edge,linewidth=0.6) +
    annotate("text",x=-sep/2-0.05,y=r+0.35,hjust=1,label=sprintf("%s (n=%d)",left_label,left_total),size=4) +
    annotate("text",x= sep/2+0.05,y=r+0.35,hjust=0,label=sprintf("%s (n=%d)",right_label,right_total),size=4) +
    annotate("text",x=-sep/2-0.45,y=0,label=left_only,size=5) +
    annotate("text",x= sep/2+0.45,y=0,label=right_only,size=5) +
    annotate("text",x=0,y=0,label=both,fontface=2,size=5) +
    coord_equal(xlim=c(-2.4,2.4),ylim=c(-1.8,1.8),expand=FALSE) + theme_void() + { if(!is.null(title)) ggtitle(title) else NULL }
}

bio <- read_tsv(bio_file, show_col_types=FALSE) %>% mutate(CHROM=canon_chrom(CHROM), LABEL=canon_label(LABEL), START=as.integer(START), END=as.integer(END))
ml  <- read_tsv(ml_file,  show_col_types=FALSE) %>% mutate(CHROM=canon_chrom(CHROM), LABEL=canon_label(LABEL), START=as.integer(START), END=as.integer(END))

labels <- c("sweep","balancing")
counts <- map_dfr(labels, ~ cluster_counts(ml, bio, .x))
write_tsv(counts, out_counts)

p_sweep <- draw_two_set("ML: sweep","Bio: sweep",
                        counts$ml_total[counts$LABEL=="sweep"], counts$bio_total[counts$LABEL=="sweep"],
                        counts$ml_only[counts$LABEL=="sweep"],  counts$bio_only[counts$LABEL=="sweep"],
                        counts$both[counts$LABEL=="sweep"])
p_bal   <- draw_two_set("ML: balancing","Bio: balancing",
                        counts$ml_total[counts$LABEL=="balancing"], counts$bio_total[counts$LABEL=="balancing"],
                        counts$ml_only[counts$LABEL=="balancing"],  counts$bio_only[counts$LABEL=="balancing"],
                        counts$both[counts$LABEL=="balancing"])

pdf(out_pdf, width=4.2, height=3.8); print(p_sweep); print(p_bal); dev.off()
ggsave(out_png_sw, p_sweep, width=4.2, height=3.8, dpi=400)
ggsave(out_png_bl, p_bal,   width=4.2, height=3.8, dpi=400)

cat("Done.\n","Counts:", out_counts, "\n","PDF:   ", out_pdf, "\n","PNG:   ", out_png_sw, "\n","PNG:   ", out_png_bl, "\n")
