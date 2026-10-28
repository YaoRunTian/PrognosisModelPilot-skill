#!/usr/bin/env Rscript
# Default TCGA acquisition; sample selection and survival derivation are separate.

download_tcga <- function(project, output_dir, files_per_chunk = 20L, backend = NULL) {
  if (length(project) != 1L || is.na(project) ||
      !grepl("^TCGA-[A-Z0-9]+$", project)) stop("Specify one TCGA-* project ID")
  if (length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))
    stop("A nonempty output directory is required")
  if (length(files_per_chunk) != 1L || is.na(files_per_chunk) ||
      !is.finite(files_per_chunk) || files_per_chunk < 1 ||
      files_per_chunk != floor(files_per_chunk)) stop("files_per_chunk must be a positive integer")
  files_per_chunk <- as.integer(files_per_chunk)
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("Install SummarizedExperiment before running this script")
  if (is.null(backend)) {
    if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
      stop("Install TCGAbiolinks before running this script")
    backend <- list(query = TCGAbiolinks::GDCquery,
                    results = TCGAbiolinks::getResults,
                    download = TCGAbiolinks::GDCdownload,
                    prepare = TCGAbiolinks::GDCprepare,
                    clinical = TCGAbiolinks::GDCquery_clinic)
  }
  for (name in c("query", "results", "download", "prepare", "clinical"))
    if (!is.function(backend[[name]])) stop("Missing backend function: ", name)
  outputs <- c("expression_SE.rds", "counts.rds", "tpm.rds", "log2_tpm_plus1.rds",
               "clinical_raw.rds", "clinical_raw.csv", "sample_metadata.rds",
               "gene_metadata.rds", "query.rds", "acquisition_manifest.rds")
  existing <- outputs[file.exists(file.path(output_dir, outputs))]
  if (length(existing)) stop("Existing outputs will not be overwritten: ", paste(existing, collapse = ", "))
  parameters <- list(project = project, data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts", access = "open")
  started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  query <- do.call(backend$query, parameters)
  file_manifest <- backend$results(query)
  if (!is.data.frame(file_manifest) || nrow(file_manifest) == 0L)
    stop("GDC query returned no files")
  if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(output_dir))
    stop("Cannot create output directory")
  cache <- file.path(output_dir, "GDCdata")
  backend$download(query, method = "api", directory = cache, files.per.chunk = files_per_chunk)
  se <- backend$prepare(query, directory = cache, summarizedExperiment = TRUE)
  if (!inherits(se, "SummarizedExperiment")) stop("GDCprepare did not return a SummarizedExperiment")
  assays <- SummarizedExperiment::assayNames(se)
  needed <- c("unstranded", "tpm_unstrand")
  if (!all(needed %in% assays))
    stop("Expected STAR counts and TPM assays; actual assays: ", paste(assays, collapse = ", "))
  counts <- SummarizedExperiment::assay(se, "unstranded")
  tpm <- SummarizedExperiment::assay(se, "tpm_unstrand")
  if (!identical(dimnames(counts), dimnames(tpm)) || any(dim(counts) == 0L))
    stop("Counts and TPM dimensions/identifiers are inconsistent or empty")
  clinical <- backend$clinical(project = project, type = "clinical")
  if (!is.data.frame(clinical) || nrow(clinical) == 0L) stop("Clinical query returned no rows")
  csv_ok <- !any(vapply(clinical, is.list, logical(1)))
  save_object <- function(x, name) saveRDS(x, file.path(output_dir, name))
  save_object(se, "expression_SE.rds")
  save_object(counts, "counts.rds")
  save_object(tpm, "tpm.rds")
  save_object(log2(tpm + 1), "log2_tpm_plus1.rds")
  save_object(clinical, "clinical_raw.rds")
  if (csv_ok) utils::write.csv(clinical, file.path(output_dir, "clinical_raw.csv"), row.names = FALSE)
  else message("Clinical list columns retained in clinical_raw.rds; CSV was not generated")
  save_object(SummarizedExperiment::colData(se), "sample_metadata.rds")
  save_object(SummarizedExperiment::rowData(se), "gene_metadata.rds")
  save_object(query, "query.rds")
  produced <- outputs[outputs != "clinical_raw.csv" | csv_ok]
  manifest <- list(project = project, provider = "TCGAbiolinks_GDC", query_parameters = parameters,
                   download = list(method = "api", files_per_chunk = files_per_chunk, directory = cache),
                   clinical_type = "clinical", clinical_csv_generated = csv_ok,
                   started_utc = started, completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                   gdc_files = file_manifest, assay_names = assays, output_files = produced,
                   samples = ncol(se), genes = nrow(se), clinical_rows = nrow(clinical),
                   sample_selection_performed = FALSE, survival_endpoint_derived = FALSE,
                   session = utils::sessionInfo())
  save_object(manifest, "acquisition_manifest.rds")
  message("Acquisition complete: ", project, "; genes=", nrow(se), "; samples=", ncol(se),
          ". Raw clinical data saved; analysis cohort and survival endpoint are not yet derived.")
  invisible(manifest)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- "Rscript download_tcga.R --project TCGA-LIHC --out data/TCGA-LIHC [--files-per-chunk 20]"
  if (identical(args, "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
  if (length(args) == 0L || length(args) %% 2L != 0L) stop(usage)
  keys <- args[seq.int(1L, length(args), 2L)]
  if (anyDuplicated(keys) || !all(keys %in% c("--project", "--out", "--files-per-chunk"))) stop(usage)
  opts <- setNames(args[seq.int(2L, length(args), 2L)], keys)
  if (!all(c("--project", "--out") %in% names(opts))) stop(usage)
  chunk <- if ("--files-per-chunk" %in% names(opts)) suppressWarnings(as.numeric(opts[["--files-per-chunk"]])) else 20L
  download_tcga(opts[["--project"]], opts[["--out"]], chunk)
}

if (sys.nframe() == 0L) main()
