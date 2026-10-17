#!/usr/bin/env Rscript

# ==============================================================================
# SMaHT mitochondrial variant evidence classification
# ==============================================================================
# Purpose:
#   Integrate per-tissue VCF calls from multiple sequencing technologies and
#   variant callers; annotate CrossTech, CrossCaller, and CrossTissue evidence;
#   assign HighConf, LowConf, or LikelyArtifact confidence classes.
#
# Usage:
#   Rscript classify_mtDNA_variants.R <vcf_directory> <output_directory>
#
# Input:
#   One uncompressed .vcf file per tissue. VCF sample columns follow the pattern
#   <tissue>-<technology>-<caller> (see CALLSETS below).
#
# Outputs:
#   variant_classification.tsv    Per-tissue variant evidence and confidence
#   tissue_summary.tsv            Summary counts by tissue
#   evidence_plot_data.tsv        Long-format CrossTech/Caller/Tissue counts
#   confidence_plot_data.tsv      Long-format confidence-class counts
#   sessionInfo.txt                Software/environment information
#
# Requirements: R >= 4.1, data.table (record tested versions in the README).
# ============================================================================== 

suppressPackageStartupMessages(library(data.table))

# ---- 1. Input arguments and configuration -----------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript classify_mtDNA_variants.R ",
    "<vcf_directory> <output_directory>",
    call. = FALSE
  )
}

vcf_dir <- normalizePath(args[[1L]], mustWork = TRUE)
output_dir <- args[[2L]]
if (!dir.exists(vcf_dir)) stop("VCF input is not a directory: ", vcf_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

# A technology/caller is considered available in a tissue if its VCF sample
# column contains at least one non-NA entry. This reproduces the original
# availability criterion (including '.' as a supplied, no-call value).
CALLSETS <- data.table(
  input_suffix = c(
    "ont-mitoscope", "pacbio-mitoscope", "Illumina-mutect2",
    "ont-himito", "pacbio-himito"
  ),
  name = c(
    "ont_mitoscope", "pacbio_mitoscope", "illumina_mutect2",
    "ont_himito", "pacbio_himito"
  ),
  technology = c("ont", "pacbio", "illumina", "ont", "pacbio"),
  caller = c("mitoscope", "mitoscope", "mutect2", "himito", "himito")
)

callset_names <- CALLSETS$name
confidence_levels <- c("HighConf", "LowConf", "LikelyArtifact")

# ---- 2. VCF reading and allele-fraction extraction --------------------------
get_tissue <- function(path) sub("\\.vcf$", "", basename(path))

read_vcf_column_names <- function(path) {
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("No #CHROM header found in: ", path)
    if (startsWith(line, "#CHROM\t")) {
      return(strsplit(line, "\t", fixed = TRUE)[[1L]])
    }
  }
}

# Extract the AF entry indicated by the VCF FORMAT field, rather than assuming
# that the last colon-separated sample value is always AF.
extract_af <- function(sample_values, format_values) {
  if (length(sample_values) != length(format_values)) {
    stop("Sample and FORMAT vectors must have the same length.")
  }

  af <- rep(NA_real_, length(sample_values))
  valid <- !is.na(sample_values) & nzchar(sample_values) &
    !(sample_values %in% c(".", "./.", ".:.")) &
    !is.na(format_values) & nzchar(format_values)

  for (format_string in unique(format_values[valid])) {
    af_position <- match(
      "AF", strsplit(format_string, ":", fixed = TRUE)[[1L]]
    )
    if (is.na(af_position)) next

    idx <- which(valid & format_values == format_string)
    sample_parts <- strsplit(sample_values[idx], ":", fixed = TRUE)
    raw_af <- vapply(
      sample_parts,
      function(parts) {
        if (length(parts) >= af_position) parts[[af_position]] else NA_character_
      },
      character(1)
    )
    raw_af[raw_af %in% c("", ".")] <- NA_character_
    af[idx] <- suppressWarnings(as.numeric(raw_af))
  }

  af
}

read_one_vcf <- function(path) {
  tissue <- get_tissue(path)
  header <- read_vcf_column_names(path)

  dt <- fread(
    path, sep = "\t", header = TRUE, skip = "#CHROM",
    check.names = FALSE, data.table = TRUE,
    na.strings = c("", "NA"), showProgress = FALSE
  )

  if (ncol(dt) != length(header)) {
    stop("Column count does not match #CHROM header in: ", path)
  }
  setnames(dt, header)

  required <- c("#CHROM", "POS", "REF", "ALT", "FORMAT")
  if (!all(required %in% names(dt))) {
    stop("Missing required VCF columns in ", path, ": ",
         paste(setdiff(required, names(dt)), collapse = ", "))
  }
  if (!nrow(dt)) {
    warning("Skipping empty VCF: ", path)
    return(NULL)
  }

  expected <- paste0(tissue, "-", CALLSETS$input_suffix)
  present <- expected %in% names(dt)
  if (!any(present)) {
    warning("Skipping ", tissue, ": no expected sample columns found.")
    return(NULL)
  }

  setnames(dt, old = expected[present], new = callset_names[present])
  for (nm in callset_names[!present]) {
    dt[, (nm) := NA_character_]
  }

  dt[, tissue := tissue]
  dt[, variant_id := paste(`#CHROM`, POS, REF, ALT, sep = "_")]

  # Availability is a tissue-level property, shared by all rows in this VCF.
  for (nm in callset_names) {
    dt[, (paste0("avail_", nm)) := any(!is.na(get(nm)))]
    dt[, (paste0("AF_", nm)) := extract_af(get(nm), FORMAT)]
  }

  dt
}

vcf_files <- sort(list.files(vcf_dir, pattern = "\\.vcf$", full.names = TRUE))
if (!length(vcf_files)) stop("No .vcf files found in: ", vcf_dir)

message("Reading ", length(vcf_files), " VCF file(s) from ", vcf_dir)
vcf_list <- Filter(Negate(is.null), lapply(vcf_files, read_one_vcf))
if (!length(vcf_list)) stop("No usable VCF inputs were found.")
all_vcf <- rbindlist(vcf_list, fill = TRUE, use.names = TRUE)

# ---- 3. Call-level detection and technology support -------------------------
# A caller detects a variant when its numeric AF is greater than zero.
# Missing/non-numeric AF values do not constitute supporting evidence.
for (nm in callset_names) {
  af_name <- paste0("AF_", nm)
  all_vcf[, (paste0("det_", nm)) :=
            !is.na(get(af_name)) & get(af_name) > 0]
}

all_vcf[, `:=`(
  tech_illumina_available = avail_illumina_mutect2,
  tech_ont_available = avail_ont_mitoscope | avail_ont_himito,
  tech_pacbio_available = avail_pacbio_mitoscope | avail_pacbio_himito,
  det_illumina = det_illumina_mutect2,
  det_ont = det_ont_mitoscope | det_ont_himito,
  det_pacbio = det_pacbio_mitoscope | det_pacbio_himito
)]

all_vcf[, `:=`(
  n_tech_available = as.integer(tech_illumina_available) +
    as.integer(tech_ont_available) + as.integer(tech_pacbio_available),
  n_tech_detected = as.integer(det_illumina) +
    as.integer(det_ont) + as.integer(det_pacbio)
)]

all_vcf[, CrossTech := n_tech_available >= 2L & n_tech_detected >= 2L]

# ---- 4. Independent caller support -----------------------------------------
all_vcf[, `:=`(
  caller_mutect2_available = avail_illumina_mutect2,
  caller_mitoscope_available = avail_ont_mitoscope | avail_pacbio_mitoscope,
  caller_himito_available = avail_ont_himito | avail_pacbio_himito,
  det_mutect2 = det_illumina_mutect2,
  det_mitoscope = det_ont_mitoscope | det_pacbio_mitoscope,
  det_himito = det_ont_himito | det_pacbio_himito
)]

all_vcf[, `:=`(
  n_caller_available = as.integer(caller_mutect2_available) +
    as.integer(caller_mitoscope_available) +
    as.integer(caller_himito_available),
  n_caller_detected = as.integer(det_mutect2) +
    as.integer(det_mitoscope) + as.integer(det_himito)
)]

all_vcf[, CrossCaller := n_caller_available >= 2L & n_caller_detected >= 2L]

# ---- 5. Cross-tissue support ------------------------------------------------
all_vcf[, detected_any :=
          det_ont_mitoscope | det_pacbio_mitoscope |
          det_illumina_mutect2 | det_ont_himito | det_pacbio_himito]

# A variant is tissue-shared when its CHROM/POS/REF/ALT combination is detected
# in at least two distinct tissues. This is variant-level evidence, annotated
# on every tissue row containing that variant, matching the original script.
all_vcf[, n_tissues_detected := uniqueN(tissue[detected_any]),
        by = variant_id]
all_vcf[, CrossTissue := n_tissues_detected >= 2L]

# ---- 6. Confidence classification -------------------------------------------
# HighConf      = CrossTech OR (CrossCaller AND CrossTissue)
# LowConf       = CrossCaller OR CrossTissue, excluding HighConf
# LikelyArtifact = all remaining rows
all_vcf[, confidence := fcase(
  CrossTech | (CrossCaller & CrossTissue), "HighConf",
  CrossCaller | CrossTissue, "LowConf",
  default = "LikelyArtifact"
)]
all_vcf[, confidence := factor(confidence, levels = confidence_levels)]

# ---- 7. Tissue-level summaries ---------------------------------------------
summary_counts <- all_vcf[, .(
  CrossTech = sum(CrossTech),
  CrossCaller = sum(CrossCaller),
  CrossTissue = sum(CrossTissue),
  HighConf = sum(confidence == "HighConf"),
  LowConf = sum(confidence == "LowConf"),
  LikelyArtifact = sum(confidence == "LikelyArtifact"),
  n_provided = sum(unlist(.SD[1L]), na.rm = TRUE)
), by = tissue, .SDcols = paste0("avail_", callset_names)]

summary_counts[, total_count := HighConf + LowConf + LikelyArtifact]
setorder(summary_counts, -n_provided, -total_count, tissue)
summary_counts[, tissue := factor(tissue, levels = tissue)]

# Long-form tables are ready for visualization or plotting in ggplot2.
plot_dt1 <- melt(
  summary_counts,
  id.vars = c("tissue", "total_count", "n_provided"),
  measure.vars = c("CrossTech", "CrossCaller", "CrossTissue"),
  variable.name = "category", value.name = "count"
)

plot_dt2 <- melt(
  summary_counts,
  id.vars = c("tissue", "total_count", "n_provided"),
  measure.vars = confidence_levels,
  variable.name = "category", value.name = "count"
)

# ---- 8. Export and run information -----------------------------------------
fwrite(all_vcf, file.path(output_dir, "variant_classification.tsv"), sep = "\t")
fwrite(summary_counts, file.path(output_dir, "tissue_summary.tsv"), sep = "\t")
fwrite(plot_dt1, file.path(output_dir, "evidence_plot_data.tsv"), sep = "\t")
fwrite(plot_dt2, file.path(output_dir, "confidence_plot_data.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Classification completed: ", nrow(all_vcf), " variant-tissue rows; ",
        uniqueN(all_vcf$tissue), " tissues.")
print(summary_counts)
message("Output files written to: ", output_dir)
message("R version: ", getRversion(),
        "; data.table version: ", as.character(packageVersion("data.table")))
