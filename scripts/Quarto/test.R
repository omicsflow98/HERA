#!/bin/env Rscript

library(dplyr)
library(stringr)
library(tidyr)
library(plotly)
library(DT)
library(kableExtra)
library(zoo)

if (vcfstate) {
  side_info <- read.table(vcf_info, sep = "\t", fill = TRUE, header = TRUE)
}

if (genestate) {
  bedfile <- read.table(bedfile, sep = "\t", fill = TRUE, header = FALSE)
  colnames(bedfile)[1:4] <- c("CHR", "START", "END", "gene_name")
  bedsubset <- bedfile %>%
    select(gene_name, CHR, START, END)
}


count_snps_in_regions <- function(regions, snps) {
  snps %>%
    select(CHR, POS) %>%
    inner_join(
      distinct(regions, CHR, POS1, POS2),
      by = join_by(CHR, between(POS, POS1, POS2))
    ) %>%
    count(CHR, POS1, POS2, name = "snp_count")
}


if (hom) {
  ROHfile <- read.table(homfile, sep = "\t", fill = TRUE, header = TRUE) %>%
    arrange(as.integer(sub("ROH_", "", ROH_ID)))
  
  write.table(ROHfile, file = "HOM_combined.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

  indiv_HOM <- read.table(indiv_homfile, sep = "\t", fill = TRUE, header = TRUE) %>%
    rename(Sample = IID) %>%
    rename(nROH = NSEG) %>%
    rename(`Total ROH Length (KB)` = KB) %>%
    rename(`Average ROH Length (KB)` = KBAVG) %>%
    rename(FROH = FRR)

  ROHislands <- read.table(HOMisland, sep = "\t", col.names = c("SNP", "CHR", "POS", "POS2", "count")) %>%
    select(SNP, CHR, POS)

  ROHislands$CHR <- as.character(ROHislands$CHR)

  prepHOMtable <- ROHfile %>%
    select(-c("PHOM", "PHET")) %>%
    mutate(
      KB = round(KB, 2),
      samples = str_count(samples, ",") + 1
    ) %>%
    select(-DENSITY) %>%
    arrange(desc(samples))

  if (genestate) {
    prepHOMtable <- prepHOMtable %>%
      mutate (genes = str_count(genes, " ") + 1)
  }

  prepHOMtable$CHR <- as.character(prepHOMtable$CHR)

  plot_ROHislands <- prepHOMtable %>%
    filter(samples > cutoff) %>%
    count_snps_in_regions(ROHislands) %>%
    filter(snp_count >= as.integer(minsnp)) %>%
    inner_join(prepHOMtable, by = c("CHR", "POS1", "POS2"))

  if (genestate) {
    plot_ROHislands <- select(plot_ROHislands, ROH_ID, POS1, POS2, samples, CHR, KB, NSNP, genes)
  } else {
    plot_ROHislands <- select(plot_ROHislands, ROH_ID, POS1, POS2, samples, CHR, KB, NSNP)
  }

  plot_table_HOM <- prepHOMtable %>%
    slice_head(n = 500)

  if (genestate) {
    num_genes_HOM <- ROHfile %>%
      mutate(gene_count = str_count(genes, " ")) %>%   # count spaces
      select(ROH_ID, gene_count) %>%
      mutate(gene_count = replace_na(gene_count, 0)) %>%  # optional
      mutate(gene_count = as.integer(gene_count)) %>%
      {
        full_range <- min(.$gene_count):max(.$gene_count)
      
        data.frame(
          gene_count = full_range,
          n = as.vector(table(factor(.$gene_count, levels = full_range)))
        )
      }
  
    num_genes_HOM$gene_count <- as.factor(num_genes_HOM$gene_count)
    num_genes_HOM$n <- as.integer(num_genes_HOM$n)

    ROHfile_long <- ROHfile %>%
    separate_rows(genes, sep = " ") %>%
    filter(genes != "") %>%  # Remove any empty strings from extra spaces
    select(genes, ROH_ID, samples)

    gene_df_ROH <- ROHfile_long %>%
    group_by(genes) %>%
    summarize(
      ROHs = paste(ROH_ID, collapse = ";"),
      samples = paste(samples, collapse = ";"),
      .groups = "drop"
    ) %>% 
    mutate(samples = sapply(samples, function(x) {
        clean_samples <- str_split(x, "[,;]")[[1]] %>%
          str_trim() %>%
          unique() %>%
          paste(collapse = ";")
        clean_samples
      })
    ) %>%
    rename(gene = genes) %>%
    select(gene, ROHs, samples)

    final_gene_ROH <- gene_df_ROH %>%
      left_join(bedsubset, by = c("gene" = "gene_name")) %>%
      select(gene, CHR, START, END, ROHs, samples)

    write.table(final_gene_ROH, file = "gene_ROH.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

    final_ROH_count <- final_gene_ROH %>%
    mutate(
      ROHs = str_count(ROHs, ";") + 1,
      samples = str_count(samples, ";") + 1
    ) %>%
    select(gene, ROHs, samples, CHR, START, END) %>%
    slice_head(n = 500)
  }

ROH_info_HOM <- ROHfile %>%
    select(ROH_ID, KB, NSNP, DENSITY)
  
  khom<- max(1, round(0.01 * nrow(ROH_info_HOM)))
  ROH_info_HOM$y1_kb <- rollmean(ROH_info_HOM$KB, k = khom, fill = NA, align = "center")
  ROH_info_HOM$y1_nsnp <- rollmean(ROH_info_HOM$NSNP, k = khom, fill = NA, align = "center")
  ROH_info_HOM$y1_dens <- rollmean(1/ROH_info_HOM$DENSITY, k = khom, fill = NA, align = "center")
  ROH_info_HOM$ROH_num <- as.numeric(sub("ROH_", "", ROH_info_HOM$ROH_ID))
  x_labels <- c(ROH_info_HOM$ROH_num[1], ROH_info_HOM$ROH_num[nrow(ROH_info_HOM)])

}

if (het) {
  HRRfile <- read.table(hetfile, sep = "\t", fill = TRUE, header = TRUE) %>%
    arrange(as.integer(sub("HRR_", "", HRR_ID)))
  
  write.table(HRRfile, file = "HET_combined.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

  indiv_HET <- read.table(indiv_hetfile, sep = "\t", fill = TRUE, header = TRUE) %>%
    rename(Sample = IID) %>%
    rename(nHRR = NSEG) %>%
    rename(`Total HRR Length (KB)` = KB) %>%
    rename(`Average HRR Length (KB)` = KBAVG) %>%
    rename(FHRR = FRR)

  HRRislands <- read.table(HETisland, sep = "\t", col.names = c("SNP", "CHR", "POS", "POS2", "count")) %>%
  select(SNP, CHR, POS)

  HRRislands$CHR <- as.character(HRRislands$CHR)

  prepHETtable <- HRRfile %>%
    select(-c("PHOM", "PHET")) %>%
    mutate(
      DENSITY = round(DENSITY, 2),
      KB = round(KB, 2),
      samples = str_count(samples, ",") + 1
    ) %>%
    select(-DENSITY) %>%
    arrange(desc(samples))
 
  if (genestate) {
    prepHETtable <- prepHETtable %>%
      mutate (genes = str_count(genes, " ") + 1)
  }

  prepHETtable$CHR <- as.character(prepHETtable$CHR)

  plot_HRRislands <- prepHETtable %>%
    filter(samples >= cutoff) %>%
    count_snps_in_regions(HRRislands) %>%
    filter(snp_count >= as.integer(minsnp)) %>%
    inner_join(prepHETtable, by = c("CHR", "POS1", "POS2"))

  if (genestate) {
    plot_HRRislands <- select(plot_HRRislands, HRR_ID, POS1, POS2, samples, CHR, KB, NSNP, genes)
  } else {
    plot_HRRislands <- select(plot_HRRislands, HRR_ID, POS1, POS2, samples, CHR, KB, NSNP)
  }

  plot_table_HET <- prepHETtable %>%
    slice_head(n = 500)
  
  if (genestate) {
    num_genes_HET <- HRRfile %>%
      mutate(gene_count = str_count(genes, " ")) %>%   # count spaces
      select(HRR_ID, gene_count) %>%
      mutate(gene_count = replace_na(gene_count, 0)) %>%  # optional
      mutate(gene_count = as.integer(gene_count)) %>%
      {
        full_range <- min(.$gene_count):max(.$gene_count)
      
        data.frame(
          gene_count = full_range,
          n = as.vector(table(factor(.$gene_count, levels = full_range)))
        )
      }
  
    num_genes_HET$gene_count <- as.factor(num_genes_HET$gene_count)
    num_genes_HET$n <- as.integer(num_genes_HET$n)

    HRRfile_long <- HRRfile %>%
    separate_rows(genes, sep = " ") %>%
    filter(genes != "") %>%  # Remove any empty strings from extra spaces
    select(genes, HRR_ID, samples)

    gene_df_HRR <- HRRfile_long %>%
    group_by(genes) %>%
    summarize(
      ROHs = paste(HRR_ID, collapse = ";"),
      samples = paste(samples, collapse = ";"),
      .groups = "drop"
    )  %>% 
    mutate(samples = sapply(samples, function(x) {
      clean_samples <- str_split(x, "[,;]")[[1]] %>%
        str_trim() %>%
        unique() %>%
        paste(collapse = ";")
      clean_samples
      })
    ) %>%
    rename(gene = genes) %>%
    select(gene, ROHs, samples)

    final_gene_HRR <- gene_df_HRR %>%
      left_join(bedsubset, by = c("gene" = "gene_name")) %>%
      select(gene, CHR, START, END, ROHs, samples)
    write.table(final_gene_HRR, file = "gene_HRR.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

    final_HRR_count <- final_gene_HRR %>%
    mutate(
      ROHs = str_count(ROHs, ";") + 1,
      samples = str_count(samples, ";") + 1
    ) %>%
    select(gene, ROHs, samples, CHR, START, END) %>%
    slice_head(n = 500)
  }
  
  ROH_info_HET <- HRRfile %>%
  select(HRR_ID, KB, NSNP, DENSITY)
  
  khet <- max(1, round(0.01 * nrow(ROH_info_HET)))
  ROH_info_HET$y1_kb <- rollmean(ROH_info_HET$KB, k = khet, fill = NA, align = "center")
  ROH_info_HET$y1_nsnp <- rollmean(ROH_info_HET$NSNP, k = khet, fill = NA, align = "center")
  ROH_info_HET$y1_dens <- rollmean(1/ROH_info_HET$DENSITY, k = khet, fill = NA, align = "center")
  ROH_info_HET$ROH_num <- as.numeric(sub("HRR_", "", ROH_info_HET$HRR_ID))
  x_labels <- c(ROH_info_HET$ROH_num[1], ROH_info_HET$ROH_num[nrow(ROH_info_HET)])
}
