library(Seurat)
library(ggplot2)

### Configuration & Paths ###
PATHS <- list(
  outputdir = "/scratch/project_2001351/Daniel/MALBAR/RESULTS/",
  intermediatedir = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/"
)

SAMPLES <- list(
  low_conc = list(
    result_file = "low_conc_results.rds",
    label = "0.04 µM",
    output_prefix = "low_conc"
  ),
  high_conc = list(
    result_file = "high_conc_results.rds",
    label = "0.4 µM",
    output_prefix = "high_conc"
  ),
  old_low_conc = list(
    result_file = "old_low_conc_results.rds",
    label = "4 µM",
    output_prefix = "old_low_conc"
  ),
  old_high_conc = list(
    result_file = "old_high_conc_results.rds",
    label = "8 µM",
    output_prefix = "old_high_conc"
  )
)

POSITIVE_QUANTILE <- 0.99
HTO_SPECIES <- c("human", "mouse")
SPECIES_LEVELS <- c("human", "mouse", "mix", "negative")

### Helper functions ###

# Stop with a message if two named objects do not contain exactly the same cell barcodes
assert_same_cells <- function(reference_cells, other_cells, reference_name, other_name) {
  missing_from_other <- setdiff(reference_cells, other_cells)
  extra_in_other <- setdiff(other_cells, reference_cells)

  if (length(missing_from_other) > 0L || length(extra_in_other) > 0L) {
    stop(
      sprintf(
        paste0(
          "Cell-name mismatch between %s and %s. ",
          "Missing from %s: %d; extra in %s: %d."
        ),
        reference_name,
        other_name,
        other_name,
        length(missing_from_other),
        other_name,
        length(extra_in_other)
      )
    )
  }

  invisible(TRUE)
}

# Convert Seurat's HTO classifications into the species labels used here
make_species_hto_labels <- function(object) {
  global_class <- as.character(object$HTO_classification.global)
  max_id <- as.character(object$HTO_maxID)

  unexpected_singlets <- global_class == "Singlet" & !(max_id %in% HTO_SPECIES)
  if (any(unexpected_singlets, na.rm = TRUE)) {
    stop(
      "Some HTO singlets have an HTO_maxID other than 'human' or 'mouse'."
    )
  }

  species_label <- rep(NA_character_, length(global_class))
  species_label[global_class == "Singlet"] <- max_id[global_class == "Singlet"]
  species_label[global_class == "Doublet"] <- "mix"
  species_label[global_class == "Negative"] <- "negative"

  if (anyNA(species_label)) {
    stop("Unrecognized or missing values in HTO_classification.global.")
  }

  names(species_label) <- colnames(object)
  factor(species_label, levels = SPECIES_LEVELS)
}

make_confusion_plot <- function(spec_ann, hto_ann) {
  stopifnot(identical(names(spec_ann), names(hto_ann)))

  cm <- table(
    Species_mRNA = factor(
      spec_ann,
      levels = c("human", "mix", "mouse")
    ),
    Species_tag = factor(
      hto_ann,
      levels = c("human", "mix", "mouse")
    )
  )

  df_cm <- as.data.frame(cm)

  plot <- ggplot(df_cm, aes(x = Species_mRNA, y = Species_tag)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    geom_text(aes(label = Freq), size = 5) +
    scale_fill_gradient(low = "white", high = "orange") +
    theme_bw() +
    theme(legend.position = "none")

  list(table = cm, data = df_cm, plot = plot)
}

save_sample_plot <- function(plot, outputdir, prefix, suffix) {
  ggsave(
    filename = paste0(prefix, suffix),
    plot = plot,
    device = "png",
    width = 5,
    height = 5,
    path = outputdir
  )
}

calculate_multiclass_metrics <- function(cm) {
 
  tp <- diag(cm)
  fp <- colSums(cm) - tp
  fn <- rowSums(cm) - tp
  support <- rowSums(cm)

  precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA_real_)

  recall <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)

  f1 <- ifelse(!is.na(precision) & !is.na(recall) & precision + recall > 0, 2 * precision * recall / (precision + recall), 0)

  per_class <- data.frame(class = c("human", "mix", "mouse"), true_positive = as.numeric(tp), false_positive = as.numeric(fp), false_negative = as.numeric(fn),
                        support = as.numeric(support), precision = as.numeric(precision), recall = as.numeric(recall), F1 = as.numeric(f1), row.names = NULL)

  # Equal weight for every class
  macro_precision <- mean(precision, na.rm = TRUE)
  macro_recall <- mean(recall, na.rm = TRUE)
  macro_F1 <- mean(f1, na.rm = TRUE)

  # Weight each class by the number of true cells
  weighted_F1 <- weighted.mean(f1, w = support, na.rm = TRUE)

  # Pool all class-specific TP, FP and FN
  micro_tp <- sum(tp)
  micro_fp <- sum(fp)
  micro_fn <- sum(fn)

  micro_precision <- micro_tp / (micro_tp + micro_fp)
  micro_recall <- micro_tp / (micro_tp + micro_fn)

  micro_F1 <- 2 * micro_precision * micro_recall / (micro_precision + micro_recall)

  accuracy <- sum(tp) / sum(cm)

  summary <- data.frame(metric = c("macro_precision", "macro_recall", "macro_F1", "weighted_F1", "micro_precision", "micro_recall", "micro_F1", "accuracy"),
    value = c(macro_precision, macro_recall, macro_F1, weighted_F1, micro_precision, micro_recall, micro_F1, accuracy))

  return(list(per_class = per_class, summary = summary))
}

### Demultiplexing function ###

process_demux_sample <- function(config, paths, positive_quantile = 0.90) {
  message("Demultiplexing sample: ", config$label)

  result_path <- file.path(
    paths$intermediatedir,
    "R_ANALYSIS",
    config$result_file
  )

  input <- readRDS(result_path)

  required_components <- c("seurat", "hto", "specAnn")
  missing_components <- setdiff(required_components, names(input))
  if (length(missing_components) > 0L) {
    stop(
      "Missing component(s) in ", config$result_file, ": ",
      paste(missing_components, collapse = ", ")
    )
  }

  smix <- input$seurat
  hto_counts <- input$hto
  spec_ann_all <- input$specAnn
  seurat_cells <- colnames(smix)

  if (is.null(colnames(hto_counts))) {
    stop("The HTO count matrix has no cell names.")
  }
  if (is.null(names(spec_ann_all))) {
    stop("specAnn has no cell names.")
  }

  if (anyDuplicated(seurat_cells)) {
    stop("The Seurat object contains duplicated cell names.")
  }
  if (anyDuplicated(colnames(hto_counts))) {
    stop("The HTO matrix contains duplicated cell names.")
  }
  if (anyDuplicated(names(spec_ann_all))) {
    stop("specAnn contains duplicated cell names.")
  }

  assert_same_cells(
    seurat_cells,
    colnames(hto_counts),
    "Seurat object",
    "HTO matrix"
  )
  assert_same_cells(
    seurat_cells,
    names(spec_ann_all),
    "Seurat object",
    "specAnn"
  )

  if (!setequal(rownames(hto_counts), HTO_SPECIES)) {
    stop(
      "Expected the HTO matrix to contain exactly the rows 'human' and 'mouse'."
    )
  }

  # Establish the Seurat cell order as the canonical order for all objects.
  hto_counts <- hto_counts[HTO_SPECIES, seurat_cells, drop = FALSE]
  spec_ann_all <- spec_ann_all[seurat_cells]

  stopifnot(
    identical(colnames(hto_counts), seurat_cells),
    identical(names(spec_ann_all), seurat_cells)
  )

  # Add the raw HTO counts as an assay and perform CLR normalization.
  smix[["HTO"]] <- CreateAssayObject(counts = hto_counts)
  smix <- NormalizeData(
    smix,
    assay = "HTO",
    normalization.method = "CLR"
  )

  # Demultiplex based on tag enrichment.
  smix <- HTODemux(
    smix,
    assay = "HTO",
    positive.quantile = positive_quantile
  )

  if (!identical(colnames(smix), seurat_cells)) {
    stop("HTODemux changed the cell order unexpectedly.")
  }

  smix$species_HTO <- make_species_hto_labels(smix)

  global_classification <- table(smix$HTO_classification.global)
  species_classification <- table(
    smix$HTO_classification.global,
    smix$species_HTO,
    useNA = "ifany"
  )

  print(global_classification)
  print(species_classification)

  # HTO distributions and scatter plot.
  p_ridge <- RidgePlot(
    smix,
    assay = "HTO",
    features = HTO_SPECIES,
    ncol = 2
  )
  print(p_ridge)

  previous_default_assay <- DefaultAssay(smix)
  DefaultAssay(smix) <- "HTO"
  p_scatter <- FeatureScatter(
    smix,
    feature1 = "mouse",
    feature2 = "human",
    group.by = "species_HTO"
  )
  DefaultAssay(smix) <- previous_default_assay
  print(p_scatter)

  save_sample_plot(
    p_scatter,
    paths$outputdir,
    config$output_prefix,
    "_Demux_Tag_Scatter.png"
  )

  # Compare RNA UMI counts among HTO singlets, doublets and negatives.
  p_rna_counts <- VlnPlot(
    smix,
    features = "nCount_RNA",
    group.by = "HTO_classification.global",
    pt.size = 0.1,
    log = TRUE
  )
  print(p_rna_counts)

  # Remove HTO-negative cells before comparison with transcript species calls.
  keep_cells <- seurat_cells[as.character(smix$HTO_classification.global) != "Negative"]
  smix_nonnegative <- subset(smix, cells = keep_cells)

  nonnegative_meta <- smix_nonnegative[[]]
  hto_ann <- as.character(nonnegative_meta$species_HTO)
  names(hto_ann) <- rownames(nonnegative_meta)
  hto_ann <- hto_ann[keep_cells]

  spec_ann <- spec_ann_all[keep_cells]

  stopifnot(
    identical(colnames(smix_nonnegative), keep_cells),
    identical(names(spec_ann), keep_cells),
    identical(names(hto_ann), keep_cells)
  )

  confusion <- make_confusion_plot(spec_ann, hto_ann)
  print(confusion$table)
  print(confusion$plot)
  cm_metrics <- calculate_multiclass_metrics(confusion$table)
  print(cm_metrics$per_class)
  print(cm_metrics$summary)

  save_sample_plot(
    confusion$plot,
    paths$outputdir,
    config$output_prefix,
    "_Demux_Confusion_Matrix.png"
  )

  list(
    seurat = smix,
    seurat_nonnegative = smix_nonnegative,
    specAnn = spec_ann,
    htoAnn = hto_ann,
    classification_tables = list(
      global = global_classification,
      species = species_classification,
      confusion = confusion$table
    ),
    plots = list(
      ridge = p_ridge,
      tag_scatter = p_scatter,
      rna_counts = p_rna_counts,
      confusion = confusion$plot
    )
  )
}

### Run all samples ###

demux_results <- setNames(
  lapply(
    names(SAMPLES),
    function(sample_id) {
      process_demux_sample(
        config = SAMPLES[[sample_id]],
        paths = PATHS,
        positive_quantile = POSITIVE_QUANTILE
      )
    }
  ),
  names(SAMPLES)
)

### Session information ###
sessionInfo()
