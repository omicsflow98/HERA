#!/bin/env Rscript
# fixROH.R — adapted from Mohamad/HERA-main/scripts/fixROH.R
# Change: samples added to an ROH/HRR because it is a subset of another
# segment are suffixed with MARK (e.g. "S12*"). Samples that originally had
# that exact segment (from splitROH's group_by) stay unmarked.

suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(yaml))

MARK <- "*"   # change here if downstream tools dislike "*"

# Split a comma-separated sample string into a character vector (NA/"" safe)
split_samples <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(str_split(x, ","))
}

option_list = list(
  make_option(c("-T", "--tsv"), default=NULL,
              help = "merged tsv file", metavar = "tsv file"),
  make_option(c("-V", "--version"), default=NULL,
              help = "nextflow process version", metavar = "process version")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tsv)) {
  print_help(opt_parser)
  stop("Please provide the necessary files ", call. = FALSE)
}

one_chrom <- read.table(opt$tsv, sep = "\t", fill = TRUE, header = TRUE,
                        stringsAsFactors = FALSE)
filename  <- stringr::str_remove(opt$tsv, "\\..*")

if (nrow(one_chrom) < 2) {
  # Nothing to combine, so no samples are added and no marks are needed
  ver <- strsplit(R.version.string, " ")[[1]][3]
  settings_list <- list(`Process Version` = list(fixROH = opt$version),
                        `Tool Version`    = list(R = ver))
  write_yaml(settings_list, file = "versions.yml")
  write.table(one_chrom, file = paste0(filename, "_fixed.tsv"),
              quote = FALSE, row.names = FALSE, sep = "\t")
  quit(save = "no", status = 0)
}

# original = exact-match samples per row (never modified)
# current  = running sample lists, kept UNMARKED while the loop runs
original <- lapply(one_chrom$samples, split_samples)
current  <- original

for (row in 1:(nrow(one_chrom) - 1)) {
  POS1 <- one_chrom$POS1[row]
  POS2 <- one_chrom$POS2[row]

  for (scroll in (row + 1):nrow(one_chrom)) {
    POSstart <- one_chrom$POS1[scroll]
    POSend   <- one_chrom$POS2[scroll]

    if (POS2 <= POSstart) {
      break
    }
    if (POS2 >= POSend) {
      # later segment is inside the current one: it inherits the current samples
      current[[scroll]] <- unique(c(current[[scroll]], current[[row]]))
    } else if (POS1 == POSstart && POS2 < POSend) {
      # current segment is inside the later one: it inherits the later samples
      current[[row]] <- unique(c(current[[row]], current[[scroll]]))
    }
  }
}

# Add the marks only at output time: anything not in the original list gets MARK
one_chrom$samples <- vapply(seq_along(current), function(i) {
  cur <- current[[i]]
  paste(ifelse(cur %in% original[[i]], cur, paste0(cur, MARK)), collapse = ",")
}, character(1))

write.table(one_chrom, file = paste0(filename, "_fixed.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")

ver <- strsplit(R.version.string, " ")[[1]][3]
settings_list <- list(`Process version` = list(fixROH = opt$version),
                      `Tool version`    = list(R = ver))
write_yaml(settings_list, file = "versions.yml")