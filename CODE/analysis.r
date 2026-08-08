library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

### Configuration & Paths ###
PATHS <- list(
  outputdir = "/scratch/project_2001351/Daniel/MALBAR/RESULTS/",
  intermediatedir = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/",
  translator = "/scratch/project_2001351/Daniel/MALBAR/DATA/CHROMIUM/translation_3M-3pgex-may-2023.txt",
  mito_homo  = "/scratch/project_2001351/Daniel/MALBAR/DATA/homomitogenes.txt",
  mito_mus   = "/scratch/project_2001351/Daniel/MALBAR/DATA/musmitogenes.txt",
  samples = list(
    low_conc = list(
      ge  = "/scratch/project_2001351/Daniel/MALBAR/DATA/CHROMIUM/P37252_1011/outs/filtered_feature_bc_matrix/",
      hto = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P37252_1015_S2_L005/umi_count/",
      empty = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P37252_1015_S2_L005/EMPTYDROPS/umi_count/",
      label = "0.04µM",
      limits = list(mt = 5, nCt_max = "NA", nCt_min = 5000, nFeat_max = "NA", nFeat_min = 2500)
    ),
    high_conc = list(
      ge  = "/scratch/project_2001351/Daniel/MALBAR/DATA/CHROMIUM/P37252_1021/outs/filtered_feature_bc_matrix/",
      hto = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P37252_1025_S4_L005/umi_count/",
      empty = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P37252_1025_S4_L005/EMPTYDROPS/umi_count/",
      label = "0.4µM",
      limits = list(mt = 5, nCt_max = "NA", nCt_min = 5000, nFeat_max = "NA", nFeat_min = 2500)
    ),
    old_low_conc = list(
    ge = "/scratch/project_2001351/Daniel/MALBAR/DATA/CHROMIUM/X3_24_038_GE/outs/per_sample_outs/X3_24_038_GE/count/sample_filtered_feature_bc_matrix/",
    hto = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P32860_1002_S3_L001/umi_count/",
    empty = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P32860_1002_S3_L001/EMPTYDROPS/umi_count/",
    label = "4µM",
    limits = list(mt = 5, nCt_max = "NA", nCt_min = 1500, nFeat_max = "NA", nFeat_min = 1000)
    ),
    old_high_conc = list(
    ge = "/scratch/project_2001351/Daniel/MALBAR/DATA/CHROMIUM/X3_24_039_GE/outs/per_sample_outs/X3_24_039_GE/count/sample_filtered_feature_bc_matrix/",
    hto = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P32860_1004_S5_L001/umi_count/",
    empty = "/scratch/project_2001351/Daniel/MALBAR/INTERMEDIATE/CITESEQ_COUNT/P32860_1004_S5_L001/EMPTYDROPS/umi_count/",
    label = "8µM",
    limits = list(mt = 5, nCt_max = "NA", nCt_min = 1500, nFeat_max = "NA", nFeat_min = 1000)
    )    
  )
)

### Helper functions ###

# Helper to load standard 10x style matrices
load_10x_matrix <- function(path, is_ge = FALSE) {
  mat <- readMM(paste0(path, "matrix.mtx.gz"))
  feats <- read.delim(paste0(path, "features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(paste0(path, "barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
  
  if (is_ge) {
    feats$V1 <- substring(feats$V1, 8) # Strip GRCm39_/GRCh38_
    barcodes$V1 <- substr(barcodes$V1, 1, 16) # Strip "-1"
  }
  
  colnames(mat) <- barcodes$V1
  rownames(mat) <- feats$V1
  return(mat)
}

# Annotate the species in a cellranger-like way, by first setting the species by highest counts, and then annotating as
# a mix if both species counts are above the nth percentile of their distributions.

annotateSpecies <- function(mat_ge) {

  ge_hc <- colSums(mat_ge[grep("ENSG", rownames(mat_ge)),])
  ge_mc <- colSums(mat_ge[grep("ENSMUSG", rownames(mat_ge)),])

  specAnn <- ifelse(
    ge_hc > ge_mc,
    "human",
    ifelse(ge_mc > ge_hc, "mouse", "ambiguous")
  )

  human_cutoff <- quantile(
    ge_hc[specAnn == "human"],
    probs = 0.05,
    na.rm = TRUE
  )

  mouse_cutoff <- quantile(
    ge_mc[specAnn == "mouse"],
    probs = 0.05,
    na.rm = TRUE
  )

  is_mix <- ge_hc >= human_cutoff &
            ge_mc >= mouse_cutoff
  specAnn[is_mix] <- "mix"
  names(specAnn) <- colnames(mat_ge)
  return(specAnn)
}

# Calculate species classification metrics

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

### Processing and plotting ###

process_sample <- function(config, translator_df, mito_genes, outputdir) {
  message("Processing sample: ", config$label)
  
  # Load data
  mat_hto <- load_10x_matrix(config$hto)
  mat_ge  <- load_10x_matrix(config$ge, is_ge = TRUE)
  mat_empty <- load_10x_matrix(config$empty)
  
  # Annotate the species and filter out ambiguous cases
  
  specAnn <- annotateSpecies(mat_ge)
  mat_ge <- mat_ge[, specAnn != "ambiguous"]
  
  # Create Seurat object
  obj <- CreateSeuratObject(counts = mat_ge, project = config$label)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mito_genes)
  
  # Plot Seurat violin plots
  
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  
  # Plot histogram of total counts per cell and lower filtering cutoff
  
  hist(obj$nCount_RNA, breaks = 100, main = paste(config$label, "- Total Counts"), xlab = "nCount_RNA")
  abline(v = config$limits$nCt_min, col = "red")
  
  # Plot a histogram of detected genes per cell and lower filtering cutoff

  hist(obj$nFeature_RNA, breaks = 100, main = paste(config$label, "- Detected genes"), xlab = "nFeature_RNA")
  abline(v = config$limits$nFeat_min, col = "red")

  # Plot a histogram of mitochondrial reads per cell and filtering cutoff

  hist(obj$percent.mt, breaks = 100, main = paste(config$label, "- Mitochondrial reads"), xlab = "nFeature_RNA")
  abline(v = config$limits$mt, col = "red")
  
  # Filtering 1
  obj <- subset(obj, subset = percent.mt < config$limits$mt & 
                  nCount_RNA > config$limits$nCt_min & 
                  nFeature_RNA > config$limits$nFeat_min)
                  
  # Calculate upper cutoffs based on median + 7 x median absolute deviation
 
   upperc <- median(obj$nCount_RNA) + 7*mad(obj$nCount_RNA)
   upperg <- median(obj$nFeature_RNA) + 7*mad(obj$nFeature_RNA)

  # Plot histogram of total counts per cell and higher filtering cutoff
  
  hist(obj$nCount_RNA, breaks = 100, main = paste(config$label, "- Total Counts"), xlab = "nCount_RNA")
  abline(v = upperc, col = "red")

  # Plot histogram of detected genes per cell and higher filtering cutoff
  
  hist(obj$nFeature_RNA, breaks = 100, main = paste(config$label, "- Total Counts"), xlab = "nCount_RNA")
  abline(v = upperg, col = "red")
  
  # Filtering 2
  obj <- subset(obj, subset = nCount_RNA < upperc & 
                              nFeature_RNA < upperg)

  # Barcode Translation & Intersection
  matches <- match(colnames(mat_hto), translator_df$V2)
  colnames(mat_hto) <- translator_df$V1[na.omit(matches)]
  
  # Check that all cell barcodes are shared between obj and mat_hto
  
  sum(!(colnames(mat_hto) %in% colnames(obj)))
  sum(!(colnames(obj) %in% colnames(mat_hto)))
  
  # Take a common subset for obj, mat_hto, and specAnn
  positive_hto_cells <- colnames(mat_hto)[Matrix::colSums(mat_hto) > 0]
  common_cells <- colnames(obj)[colnames(obj) %in% positive_hto_cells]
  obj <- obj[, common_cells]
  mat_ge <- mat_ge[, common_cells, drop = FALSE]
  mat_hto <- mat_hto[1:2, common_cells, drop = FALSE]
  specAnn <- specAnn[common_cells]
  
  # The mouse and human oligos are swapped between the new and old set of data
  if (config$label %in% c("4µM", "8µM")) {
    rownames(mat_hto) <- c("mouse", "human")
  }
  else {
    rownames(mat_hto) <- c("human", "mouse")
  }
  
  # Separate human and mouse counts in the cells common between GE and tag matrices
  ge_hc <- colSums(mat_ge[grep("ENSG", rownames(mat_ge)), common_cells])
  ge_mc <- colSums(mat_ge[grep("ENSMUSG", rownames(mat_ge)), common_cells])
  
  # Species annotation by tags. Make sure that obj, mat_ge, mat_hto, specAnn, and htoAnn have the same cells in the same order.
  htoAnn <- ifelse(mat_hto["human", ] > 2 * mat_hto["mouse", ], "human", ifelse(mat_hto["mouse", ] > 2 * mat_hto["human", ], "mouse", "mix"))
  names(htoAnn) <- colnames(mat_hto)
  stopifnot(
    identical(colnames(obj), common_cells),
    identical(colnames(mat_ge), common_cells),
    identical(colnames(mat_hto), common_cells),
    identical(names(specAnn), common_cells),
    identical(names(htoAnn), common_cells)
  )
  
  # Species scatter plot
  spec_cols <- c("human" = "black", "mouse" = "green", "mix" = "red")
  spec_factor <- factor(specAnn, levels = names(spec_cols))
  plot(ge_mc, ge_hc, col = spec_cols[as.character(spec_factor)], 
       main = paste(config$label, "- Species Separation"),
       xlab = "Mouse counts", ylab = "Human counts")
  legend("topright", legend=names(spec_cols), text.col=spec_cols)
  
  # Species scatter plot with ggplot
  df_sc <- data.frame(Mouse=ge_mc, Human=ge_hc, Species_mRNA=specAnn, Species_tag=htoAnn)
  p_sc <- ggplot(data = df_sc) + geom_point(mapping = aes(x = Mouse, y = Human, colour = Species_tag)) +
	scale_color_manual(values = c("mix" = "purple","human" = "blue", "mouse" = "red")) + coord_fixed(ratio=1) +
	theme(legend.position = c(.9, .9), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
        theme(legend.text = element_text(size = 10), legend.title = element_text(size = 11))
  print(p_sc)
  
  # Save the plot
  ggsave(filename = paste(config$label, "_Species_Separation.png"), plot = p_sc, device = "png", width = 5, height = 5, path = outputdir)

  # Confusion matrix
  cm <- table(specAnn, htoAnn)
  df1 <- data.frame(cm)
  colnames(df1) <- c("Species_mRNA", "Species_tag", "Freq")
  p_cm <- ggplot(data =  df1, mapping = aes(x = Species_mRNA, y = Species_tag)) +
          geom_tile(aes(fill = Freq), colour = "white") +
            geom_text(size = 13, aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
            theme_bw() + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank()) +
            scale_fill_gradient(low = "white", high = "orange")
  print(p_cm)
    
  # Save the plot
  ggsave(filename = paste(config$label, "_Confusion_Matrix.png"), plot = p_cm, device = "png", width = 5, height = 5, path = outputdir)
  
  # Calculate classification metrics
  cm_metrics <- calculate_multiclass_metrics(cm)
  print(cm_metrics$per_class)
  print(cm_metrics$summary)
  
  # Distribution of total HTO counts

  print(summary(colSums(mat_hto)))
  hist(colSums(mat_hto), breaks = 100, main = paste(config$label, "- Total HTO counts per cell"))

  # Distribution of mouse HTO counts in mouse cells

  hist(mat_hto["mouse",specAnn == "mouse"], breaks = 100, main = paste(config$label, "- Mouse HTO counts in mouse cells"))

  # Distribution of human HTO counts in mouse cells

  hist(mat_hto["human",specAnn == "mouse"], breaks = 100, main = paste(config$label, "- Human HTO counts in mouse cells"))

  # Distribution of human HTO counts in human cells

  hist(mat_hto["human",specAnn == "human"], breaks = 100, main = paste(config$label, "- Human HTO counts in human cells"))

  # Distribution of mouse HTO counts in human cells

  hist(mat_hto["mouse",specAnn == "human"], breaks = 100, main = paste(config$label, "- Mouse HTO counts in human cells"))

  # What is the distribution of mouse vs human tag count differences (proportionally) in non-mixed cells?
  
  nonmixed <- colnames(mat_hto)[specAnn != "mix"]
  hist(log2(mat_hto["human",nonmixed] / mat_hto["mouse",nonmixed]), breaks=seq(-10,10,1), xlim=c(-10, 10), main = paste(config$label, "- Human vs mouse tags per cell"))
  axis(1, at = seq(-10, 10, by = 2))
  
  df_h <- data.frame(props=log2(mat_hto["human",nonmixed] / mat_hto["mouse",nonmixed]))
  p_h <- ggplot(df_h, aes(x=props)) +  geom_histogram(aes(y=..density..), colour="black", fill="white") + geom_density(alpha=.2, fill="#FF6666") +
         scale_x_continuous(breaks=seq(-10,10,1)) + theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) 
  plot(p_h)
  ggsave(filename = paste(config$label, "_mouse_vs_human_tags_hist.png"), plot = p_h, device = "png", width = 5, height = 5, path = outputdir)

  # Scatter plot of human vs mouse tag counts

  plot(log10(mat_hto["human", nonmixed]), log10(mat_hto["mouse", nonmixed]), main = paste(config$label, "- Human vs mouse tags per cell (log10)"))
  plot(mat_hto["human", nonmixed], mat_hto["mouse", nonmixed], main = paste(config$label, "- Human vs mouse tags per cell"))

  # Tag scatter plot with ggplot
  df_sc2 <- data.frame(Mouse=mat_hto["mouse", nonmixed], Human=mat_hto["human", nonmixed], mRNA_counts=obj$nCount_RNA[nonmixed])
  p_sc2 <- ggplot(data = df_sc2) + geom_point(mapping = aes(x = Mouse, y = Human, colour = mRNA_counts)) + scale_colour_viridis_c(option="magma") +
           theme(legend.position = c(.8, .8), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) 
  print(p_sc2)
  ggsave(filename = paste(config$label, "_mouse_vs_human_tags_scat.png"), plot = p_sc2, device = "png", width = 5, height = 5, path = outputdir)

  # Scatter plot of total mRNA counts vs total tag counts

  plot(ge_mc + ge_hc, colSums(mat_hto))

  # Distribution of total HTO counts in empty droplets

  print(summary(colSums(mat_empty)))
  hist(colSums(mat_empty), breaks = 100, main = paste(config$label, "- Total tag counts per empty droplet"))

  return(list(seurat = obj, hto = mat_hto, specAnn = specAnn, htoAnn = htoAnn))
}

### Run code ###

# Load Global Requirements
translator <- read.delim(PATHS$translator, header = FALSE, stringsAsFactors = FALSE)
mitogenes <- c(read.delim(PATHS$mito_mus, header = FALSE)$V1, 
               read.delim(PATHS$mito_homo, header = FALSE)$V1)

# Process all samples into a list
results <- lapply(PATHS$samples, function(s) process_sample(s, translator, mitogenes, PATHS$outputdir))

# Create a directory for saved R objects
results_dir <- file.path(PATHS$intermediatedir, "R_ANALYSIS")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Save one complete result list per sample
for (sample_id in names(results)) {
  saveRDS(
    results[[sample_id]],
    file = file.path(
      results_dir,
      paste0(sample_id, "_results.rds")
    )
  )
}

# Get session info

sessionInfo()
