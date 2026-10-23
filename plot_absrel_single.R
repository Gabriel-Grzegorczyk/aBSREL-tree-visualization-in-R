# aBSREL tree visualization: one tree and one JSON file per run.
# Edit USER SETTINGS below, then run the entire script. See README.md.
# Adapted from https://github.com/weikevinhc/kakstree (hyphy.tree).
# Original reference: Zakerzade et al. (2025), doi:10.1371/journal.pgen.1011549.
# Adaptation: Gabriel Valaski Grzegorczyk.
# Please cite the original article and acknowledge this adaptation.
# Associated study: Grzegorczyk and colleagues, unpublished manuscript.
# See README.md for citation details and updates.

# ======================= USER SETTINGS =======================
# Use forward slashes in Windows paths, e.g. "C:/Users/name/project".
data_dir <- "."
family <- "CYP9R"
treefile <- file.path(data_dir, paste0(family, "_geneTree.newick"))
jsonfile <- file.path(data_dir, paste0(family, "_CDS_ALINHADO_ABSREL.JSON"))
output_pdf <- file.path(data_dir, paste0(family, "_absrel.pdf"))
output_csv <- file.path(data_dir, paste0(family, "_absrel_branches.csv"))
export_table <- TRUE
p_cutoff <- 0.05
strong_cutoff <- 0.01
branch_width <- 10
label_size <- 1
label_offset <- 0.03
label_adj_x <- 0.05
pdf_width <- 12

# Set to NULL for black labels, or customize the prefix-to-color mapping.
species_colors <- c(Am = "#FFD700", Bt = "#FFA500", Fv = "#4169E1",
                    Mb = "#8B4513", Mr = "#2E8B57", Mq = "#A52A2A", Ta = "#FF4500")

# ======================= FUNCTION DEFINITIONS =======================
# Match internal branches by their sorted descendant tip labels.
.absrel_keys <- function(tree) {
  vapply(tree$edge[, 2], function(node) {
    if (node <= length(tree$tip.label)) return(tree$tip.label[node])
    tips <- phangorn::Descendants(tree, node, type = "tips")[[1]]
    paste(sort(tree$tip.label[tips]), collapse = "|||")
  }, character(1))
}

# Read a scalar safely, including missing/null JSON fields.
.absrel_number <- function(record, field) {
  value <- record[[field]]
  if (is.null(value) || length(value) != 1L || is.list(value)) return(NA_real_)
  suppressWarnings(as.numeric(value))
}

# Plot one external tree using one aBSREL JSON file (partition 0 only).
# Colors: baseline MG94xREV omega, clipped to [0.1, 10], then log10.
# Missing omega: use omega = 1 ONLY for plotting (not a biological estimate).
# Stars: superp <= corrected p < p (*); corrected p < superp (**).
# The function draws without a legend, title, scale bar or rerooting.
# File export is handled in RUN ONE TREE below.
# tip_colors: optional named vector; names must be complete normalized tip labels.
# Returns invisibly: list(tree = phylo, branches = data.frame).
plot_absrel <- function(treefile, jsonfile, p = 0.05, superp = 0.01,
                        branch_lwd = 10, tip_cex = 1,
                        tip_offset = 0.03, tip_adj_x = 0.05,
                        tip_colors = NULL) {
  packages <- c("ape", "phytools", "phangorn", "rjson", "RColorBrewer")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  if (length(p) != 1L || length(superp) != 1L ||
      !is.finite(p) || !is.finite(superp) || !(0 < superp && superp <= p && p <= 1)) {
    stop("Thresholds must satisfy 0 < superp <= p <= 1.")
  }
  if (!file.exists(treefile) || !file.exists(jsonfile)) stop("Input file not found.")
  normalize <- function(x) gsub("-", "_", x, fixed = TRUE)
  check_tree <- function(tree) {
    if (!inherits(tree, "phylo")) stop("Provide exactly one tree.")
    tree$tip.label <- normalize(tree$tip.label)
    if (anyNA(tree$tip.label) || any(!nzchar(tree$tip.label)) ||
        anyDuplicated(tree$tip.label) || any(grepl("|||", tree$tip.label, fixed = TRUE))) {
      stop("Tip labels must be unique and nonempty after normalization, without |||.")
    }
    ape::reorder.phylo(tree, "cladewise")
  }
  tree <- check_tree(ape::read.tree(treefile))
  json <- rjson::fromJSON(file = jsonfile)
  trees <- json$input$trees
  attributes <- json$`branch attributes`
  if (!setequal(setdiff(names(trees), "attributes"), "0") ||
      !setequal(setdiff(names(attributes), "attributes"), "0") ||
      is.null(trees[["0"]]) || is.null(attributes[["0"]])) {
    stop("Expected a single data partition named '0'; metadata 'attributes' is ignored.")
  }
  newick <- sub(";+[[:space:]]*$", "", trimws(trees[["0"]]))
  json_tree <- check_tree(ape::read.tree(text = paste0(newick, ";")))
  if (!setequal(tree$tip.label, json_tree$tip.label)) stop("Tree and JSON tip sets differ.")
  keys <- .absrel_keys(tree)
  json_keys <- .absrel_keys(json_tree)
  if (anyDuplicated(keys) || anyDuplicated(json_keys)) stop("Ambiguous descendant signatures.")
  if (!setequal(keys, json_keys)) {
    warning("Descendant sets differ: check topology/root representation. Unmatched branches retain NA.", call. = FALSE)
  }
  n <- length(json_tree$tip.label)
  edge_names <- vapply(json_tree$edge[, 2], function(node) {
    if (node <= n) return(json_tree$tip.label[node])
    i <- node - n
    if (is.null(json_tree$node.label) || i > length(json_tree$node.label)) return(NA_character_)
    normalize(json_tree$node.label[i])
  }, character(1))
  valid_names <- edge_names[!is.na(edge_names) & nzchar(edge_names)]
  if (anyDuplicated(valid_names)) stop("Duplicate branch identifiers in the JSON tree.")
  records <- attributes[["0"]]
  record_names <- normalize(names(records))
  if (!length(record_names) || anyNA(record_names) ||
      any(!nzchar(record_names)) || anyDuplicated(record_names)) stop("Ambiguous JSON record identifiers.")
  branch_names <- edge_names[match(keys, json_keys)]
  index <- match(branch_names, record_names)
  index[is.na(branch_names) | !nzchar(branch_names)] <- NA_integer_
  if (all(is.na(index))) stop("No branches could be matched to JSON records.")
  omega <- unname(vapply(records, .absrel_number, numeric(1),
                         field = "Baseline MG94xREV omega ratio")[index])
  corrected_p <- unname(vapply(records, .absrel_number, numeric(1),
                               field = "Corrected P-value")[index])
  if (any(omega < 0, na.rm = TRUE)) stop("Negative omega values in JSON.")
  if (any(!is.na(corrected_p) & (!is.finite(corrected_p) | corrected_p < 0 | corrected_p > 1))) {
    stop("Corrected p-values must lie in [0, 1].")
  }
  if (anyNA(omega) || anyNA(corrected_p)) {
    warning(sprintf(paste0("%d branch(es) lack omega: plotted with the omega=1 color. ",
                           "%d lack corrected p: no stars. Original NA values are retained."),
                    sum(is.na(omega)), sum(is.na(corrected_p))), call. = FALSE)
  }
  omega_plot <- omega
  omega_plot[is.na(omega_plot)] <- 1
  omega_plot <- log10(pmin(10, pmax(0.1, omega_plot)))
  stars <- rep("", length(corrected_p))
  stars[which(corrected_p < p)] <- "*"
  stars[which(corrected_p < superp)] <- "**"
  label_colors <- rep("black", length(tree$tip.label))
  if (!is.null(tip_colors)) {
    if (is.null(names(tip_colors)) || anyDuplicated(names(tip_colors))) {
      stop("tip_colors must have unique names matching complete tip labels.")
    }
    selected <- unname(tip_colors[tree$tip.label])
    label_colors[!is.na(selected)] <- selected[!is.na(selected)]
  }
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0) || !any(tree$edge.length > 0)) {
    stop("External tree requires finite, nonnegative branch lengths, with at least one positive length.")
  }
  depth <- max(ape::node.depth.edgelength(tree))
  palette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))
  old_par <- graphics::par(mar = c(1, 1, 1, 8))
  on.exit(graphics::par(old_par), add = TRUE)
  phytools::plotBranchbyTrait(
    tree, omega_plot, mode = "edges", palette = palette,
    show.tip.label = FALSE, xlims = c(-1, 1), legend = FALSE,
    edge.width = branch_lwd, lend = 1, x.lim = c(0, depth * 1.3)
  )
  ape::tiplabels(text = tree$tip.label, col = label_colors, frame = "none",
                 offset = tip_offset, adj = c(tip_adj_x, 0.2), cex = tip_cex)
  significant <- which(nzchar(stars))
  if (length(significant)) ape::edgelabels(
    text = stars[significant], edge = significant, frame = "none",
    adj = c(0.5, 0.75), cex = 2
  )
  branches <- data.frame(
    edge = seq_len(nrow(tree$edge)), parent = tree$edge[, 1], child = tree$edge[, 2],
    json_branch = branch_names, matched = !is.na(index), descendants = keys,
    baseline_omega = omega, corrected_p = corrected_p,
    omega_missing = is.na(omega), p_missing = is.na(corrected_p),
    plot_log10_omega = omega_plot, stars = stars
  )
  invisible(list(tree = tree, branches = branches))
}

# ======================= RUN ONE TREE =======================
# Check dependencies before reading inputs or opening the PDF.
required <- c("ape", "phytools", "phangorn", "rjson", "RColorBrewer")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))
if (!file.exists(treefile) || !file.exists(jsonfile)) stop("Input file not found. Check USER SETTINGS.")
input_tree <- ape::read.tree(treefile)
if (!inherits(input_tree, "phylo")) stop("Provide exactly one tree.")
labels <- gsub("-", "_", input_tree$tip.label, fixed = TRUE)
tip_colors <- NULL
if (!is.null(species_colors)) {
  tip_colors <- setNames(unname(species_colors[substr(labels, 1, 2)]), labels)
}

# Closing the device is guaranteed even if plotting fails.
result <- local({
  grDevices::pdf(output_pdf, width = pdf_width, height = max(6, length(labels) * 0.25))
  tryCatch(
    plot_absrel(treefile, jsonfile, p = p_cutoff, superp = strong_cutoff,
                branch_lwd = branch_width, tip_cex = label_size,
                tip_offset = label_offset, tip_adj_x = label_adj_x,
                tip_colors = tip_colors),
    finally = grDevices::dev.off()
  )
})
if (export_table) {
  write.csv(result$branches, output_csv, row.names = FALSE, na = "NA")
}
message("PDF saved: ", output_pdf)
if (export_table) message("Branch table saved: ", output_csv)

# Optional inspection:
# result$branches
# subset(result$branches, !matched | omega_missing | p_missing)
