# aBSREL tree visualization in R

Plot one phylogenetic tree using existing HyPhy aBSREL results. This single-file R script generates a clean PDF with branch colors based on baseline MG94×REV omega and asterisks based on corrected aBSREL p-values.

The figure contains tip labels and significance markers, with no legend, title or scale bar. An optional CSV preserves the values used for each branch.

The script **does not run aBSREL, estimate omega or recalculate p-value corrections**. It reads results already present in a HyPhy JSON file. It does not use `jq`.

## Requirements

Install these R packages once:

```r
install.packages(c("ape", "phytools", "phangorn", "rjson", "RColorBrewer"))
```

RStudio is optional. Input data are not included.

## Quick start

1. Download `plot_absrel_single.R`.
2. Place your Newick tree and aBSREL JSON in a directory accessible to R.
3. Edit the `USER SETTINGS` section at the top of the script:

```r
data_dir <- "C:/Users/your_name/project"
family <- "CYP9R"
```

By default, these settings expect:

```text
CYP9R_geneTree.newick
CYP9R_CDS_ALINHADO_ABSREL.JSON
```

For different filenames, edit `treefile` and `jsonfile` directly. File extensions and capitalization must match your files. `data_dir <- "."` refers to the current R working directory (`getwd()`), not automatically the script directory.

4. Run the entire script in R/RStudio. Alternatively, after editing its settings:

```bash
Rscript plot_absrel_single.R
```

Each run processes one tree and one JSON. Running the script again with the same output paths replaces the previous output files.

## Outputs

Files are saved in `data_dir` by default:

| Output | Content |
| --- | --- |
| `CYP9R_absrel.pdf` | Vector tree figure |
| `CYP9R_absrel_branches.csv` | Branch-level data, when `export_table = TRUE` |

After execution, `result$tree` holds the plotted tree and `result$branches` holds the branch table. Original missing values remain `NA` in the table.

```r
# Inspect unmatched branches or missing results.
subset(result$branches, !matched | omega_missing | p_missing)
```

Set `export_table <- FALSE` to skip CSV export. The table remains available in `result$branches`.

## Settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `family` | `"CYP9R"` | Prefix for default input/output names |
| `export_table` | `TRUE` | Write the branch table to CSV |
| `p_cutoff` | `0.05` | Corrected p threshold for significance |
| `strong_cutoff` | `0.01` | Corrected p threshold for two asterisks |
| `branch_width` | `10` | Plotted branch width |
| `label_size` | `1` | Tip-label text size |
| `label_offset` | `0.03` | Tip-label offset |
| `label_adj_x` | `0.05` | Horizontal tip-label alignment |
| `pdf_width` | `12` | PDF width in inches |
| `species_colors` | Named vector | Tip-label colors based on two-character prefixes |

PDF height is set to `max(6, number_of_tips * 0.25)` inches. Thresholds must satisfy `0 < strong_cutoff <= p_cutoff <= 1`.

The example label-color map uses `Am`, `Bt`, `Fv`, `Mb`, `Mr`, `Mq` and `Ta`. Edit it for your data, or set `species_colors <- NULL` for black labels. Unrecognized prefixes are black. This mapping changes tip-label colors, not branch colors.

## What the colors and asterisks mean

### Branch colors

Colors use the branch attribute `Baseline MG94xREV omega ratio`. This is a branch-specific estimate under the baseline model, not the full distribution of site classes inferred by aBSREL and not a weighted mean calculated from those classes.

For visualization, omega is clipped to the interval [0.1, 10] and then log10-transformed:

```r
log10(pmin(10, pmax(0.1, omega)))
```

The reversed `RdYlBu` palette maps low values to blue, omega = 1 to the light central color, and high values to red. Values at or below 0.1 share the lower endpoint; values at or above 10 share the upper endpoint. Original values are retained in the table.

### Significance markers

Markers use `Corrected P-value` directly from the JSON. With the defaults:

| Marker | Condition |
| --- | --- |
| `**` | Corrected p < 0.01 |
| `*` | 0.01 <= corrected p < 0.05 |
| None | Corrected p >= 0.05, or missing p |

Branch color is not a significance test. A low baseline omega can coexist with a significant branch-site test. No uncorrected p-values are substituted when corrected p is unavailable.

### Missing values

To retain the chosen plotting convention, **missing omega is displayed using the omega = 1 color**. This is only a graphical placeholder: it is not an estimate and does not establish neutrality. The figure alone cannot distinguish such a branch from one with an estimated omega of 1.

Missing corrected p-values receive no asterisk. If corrected p is available and significant but omega is missing, an asterisk can still be drawn. No omega or p-value is copied from a sister branch.

Warnings report missing values. Check the branch table before interpreting the figure, and explain this convention in a figure caption if any displayed omega values are missing. No gray category is added.

## Matching branches

The external Newick file supplies the plotted topology, root representation and branch lengths. The tree stored in `input$trees[["0"]]` provides the identifiers used by HyPhy.

Branches are matched by their sorted sets of descendant tip labels. This allows internal labels in the external tree to differ from HyPhy node identifiers. Hyphens are normalized to underscores for matching; normalization must not create duplicate names.

The script supports one data partition named `0`. The `attributes` metadata entry is not counted as a partition.

Matching requires compatible descendant sets. Different root representations, extra root-adjacent segments or topology differences can leave branches unmatched even when both trees contain the same tips. The script warns about these differences; it does not reroot trees or transfer results via complementary bipartitions.

## Branch table

| Column | Meaning |
| --- | --- |
| `edge` | Row index in the returned tree's edge matrix |
| `parent`, `child` | Node numbers defining that edge |
| `json_branch` | Corresponding HyPhy branch identifier |
| `matched` | Whether a JSON result record was found |
| `descendants` | Sorted descendant labels, joined by `|||` |
| `baseline_omega` | Original numeric baseline omega, including NA |
| `corrected_p` | Original numeric corrected p-value, including NA |
| `omega_missing`, `p_missing` | Missing-value indicators |
| `plot_log10_omega` | Transformed value used for branch coloring |
| `stars` | Significance marker |

Edge/node numbers belong to the returned tree representation and should not be treated as stable identifiers across different trees.

## Troubleshooting

- **Input file not found:** check `data_dir`, filenames and capitalization.
- **Tree and JSON tip sets differ:** check that the tree corresponds to the alignment analyzed by aBSREL.
- **Descendant sets differ:** inspect root representation and topology; identical tip sets do not guarantee identical branch mapping.
- **Missing omega or corrected p:** inspect `result$branches` and the corresponding JSON fields. A matched record may still lack a requested field.
- **Ambiguous identifiers:** resolve duplicate tip names, duplicate branch identifiers or duplicate descendant signatures.
- **Invalid branch lengths:** the external tree must contain finite, nonnegative lengths, with at least one positive length.
- **Clipped labels:** increase `pdf_width`, reduce `label_size`, or adjust label spacing.

## Provenance and citation

Adapted by **Gabriel Valaski Grzegorczyk** from the `hyphy.tree` workflow in [weikevinhc/kakstree](https://github.com/weikevinhc/kakstree).

The upstream repository requests citation of:

Zakerzade R, Chang CH, Chatla K, Krishnapura A, Appiah SP, et al. (2025). *Diversification and recurrent adaptation of the synaptonemal complex in Drosophila*. PLOS Genetics 21(1): e1011549. https://doi.org/10.1371/journal.pgen.1011549

Changes include descendant-based matching, corrected p-values for significance markers, robust extraction of missing fields, preserved original values, a branch-level table, configurable label colors, English documentation and a single-tree workflow. Optional rerooting and the multi-family loop were removed. The omega = 1 plotting convention for missing omega was retained.

When reporting use of this adaptation, acknowledge the upstream code and its requested reference, and identify the exact version or commit of this adaptation. Cite the aBSREL method separately for the analysis itself. Repository and release citation metadata should be added once publication details are finalized.

## Validation and licensing status

This consolidated script has not been executed in R in the preparation environment. Validate it with your tree/JSON pair before treating it as a tested release, especially branch matching, missing-value reports and significant-branch counts against the HyPhy results. No tested R/package version matrix is currently provided.

Upstream redistribution and licensing terms remain to be clarified. No new license is asserted over the original code. Before public redistribution, resolve the applicable terms and include the appropriate license and notices. Attribution does not replace permission.
