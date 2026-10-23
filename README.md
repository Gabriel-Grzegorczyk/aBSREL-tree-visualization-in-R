# aBSREL tree visualization

An R script adapted from [kakstree](https://github.com/weikevinhc/kakstree), the plotting code associated with:

**Zakerzade et al. (2025).** *Diversification and recurrent adaptation of the synaptonemal complex in Drosophila*. PLOS Genetics 21(1): e1011549. [https://doi.org/10.1371/journal.pgen.1011549](https://doi.org/10.1371/journal.pgen.1011549)

**If you use this adaptation, please cite the original article and acknowledge this repository and the associated study by Gabriel Valaski Grzegorczyk and colleagues on detoxification-related gene superfamilies in Meliponini bees (unpublished manuscript). Please cite the associated paper once it becomes publicly available.** Until then, identify the repository URL and the version or commit used; publication details will be added here when available.

Adaptation by **Gabriel Valaski Grzegorczyk**. Changes include branch matching by descendant tips, corrected p-values for significance markers, and a single-tree workflow with optional CSV export.

## What it does

The script plots one external Newick tree using existing HyPhy aBSREL JSON results. It generates a clean PDF without a legend and optionally exports a branch-level CSV. It does not run aBSREL or recalculate statistical tests.

- **Colors:** baseline MG94×REV omega, clipped to [0.1, 10] and log10-transformed. Colors do not represent the full aBSREL omega-class distribution.
- **Asterisks:** `*` for 0.01 ≤ corrected p < 0.05 and `**` for corrected p < 0.01 by default.
- **Missing values:** missing omega uses the omega = 1 color only as a graphical placeholder, not an estimate of neutrality. Missing p receives no asterisk. Original values, including `NA`, are preserved in the table.

## Files

| File | Description |
| --- | --- |
| `plot_absrel_single.R` | Single-file plotting script |
| `CYP9R_geneTree.newick` | Example CYP9R tree |
| `CYP9R_CDS_ALINHADO_ABSREL.JSON` | Example aBSREL results |

The tree and JSON reproduce the visualization. Rerunning aBSREL also requires the codon alignment and analysis settings.

## Run the example

Install the required R packages once:

```r
install.packages(c("ape", "phytools", "phangorn", "rjson", "RColorBrewer"))
```

Download the repository, set its folder as your R working directory, and run the entire script. The default settings are:

```r
data_dir <- "."
family <- "CYP9R"
```

Alternatively, set `data_dir` to the full path of your input folder. For other datasets, change `family` or edit `treefile` and `jsonfile` directly.

Outputs are `CYP9R_absrel.pdf` and `CYP9R_absrel_branches.csv`. Set `export_table <- FALSE` to disable CSV export. Reusing output paths overwrites existing files.

## Customize and check

Edit the `USER SETTINGS` section to change significance thresholds, branch width, label size, PDF width or species-label colors. Set `species_colors <- NULL` for black labels.

The external tree and JSON must have matching tip sets and compatible descendant sets. Differences in topology or root representation can leave branches unmatched. Only partition `0` is supported; no rerooting is performed. Inspect warnings and missing values with:

```r
subset(result$branches, !matched | omega_missing | p_missing)
```

## License status

Upstream redistribution terms remain to be clarified. No new license is asserted over the original code; citation does not replace permission to redistribute it.
