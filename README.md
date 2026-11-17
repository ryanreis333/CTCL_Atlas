# CTCL Atlas

Code and intermediate results for:

> **Large-scale integrated single-cell analysis reveals molecular subtypes and microenvironmental remodeling in cutaneous T cell lymphoma**
> Reis RJ, Mou E, Liu V, Brunner PM, Borcherding N, Jabbari A. *(manuscript in preparation)*

We integrated public and contributed single-cell RNA-seq (scRNA-seq) data from cutaneous T cell lymphoma (CTCL), benign inflammatory dermatoses (BID) and healthy control (HC) skin and blood into an atlas of ~1.57 million cells. The analysis then:

- identifies malignant T cells using expression and TCR data;
- tests for recurrent malignant gene expression changes with sample-aware pseudobulk mixed models;
- defines two molecular subtypes of skin CTCL;
- characterizes tumor-microenvironment (TME) remodeling;
- evaluates the subtype signatures in an independent bulk RNA-seq survival cohort (GSE168508).

---

## Repository layout

```
CTCL_Atlas/
├── final/                  # one folder per main figure
│   ├── fig1/               # atlas overview, cell-type annotation, composition
│   ├── fig2/               # T/ILC compartment, clonality, example samples
│   ├── fig3/               # cell-level (per-sample) malignant vs. reactive DEG
│   ├── fig4/               # pseudobulk CTCL vs. CD4 DEG, consensus-clustered subtypes
│   │   └── surv.Rmd        #   preparation of the GSE168508 validation cohort
│   ├── fig5/               # skin TME: CTCL vs. BID vs. HC
│   ├── fig6/               # TME by lesion type and molecular subtype
│   └── fig7/               # clinical trajectories, signature scoring, survival
│       ├── figN.Rmd        #   analysis notebook
│       ├── results/        #   intermediate results (.rds / .csv) used by the figures
│       └── figs/           #   rendered figure panels (.pdf / .png)
├── int/                    # small shared inputs (GSE168508 matrix + metadata, subtype lists)
├── docs/
│   └── software_refs.bib   # citations for the software used
├── renv/ , renv.lock       # pinned R environment (R 4.5.3)
├── .Rprofile               # activates renv
└── CTCL.Rproj
```

### Figure → notebook → key results

| Figure | Notebook | Key outputs in `results/` |
|---|---|---|
| 1 | `final/fig1/fig1.Rmd` | `markers.celltype.all.rds` (lineage markers) |
| 2 | `final/fig2/fig2.Rmd` | figure panels only |
| 3 | `final/fig3/fig3.Rmd` | `deg_all.*` (pooled Wilcoxon), `deg.per.sample.all.rds`, `shrunk_mat_sig.rds` (ashr-shrunken log2FC), `gene_summary.csv`, `gene_stats.csv` |
| 4 | `final/fig4/fig4.Rmd` | `bulk.se.rds` (pseudobulk), `res_df_cd4.*` (CTCL vs. CD4 dream), `varPart.csv`, `gsea_res_cd4.*`, `subtype_df.rds` (subtype labels), `res_skin_subtypes.rds` (subtype DEG), `functional_enrichment*.rds` |
| 4 / 7 | `final/fig4/surv.Rmd` | builds `int/surv.mat.rds`, `int/surv.meta*.rds` from GSE168508 |
| 5 | `final/fig5/fig5.Rmd` | `da_prop.*` (propeller DA), `tme_bulk_deg.*`, `deg_celltype_*` (cell-type pseudobulk DEG), `gsea_*` |
| 6 | `final/fig6/fig6.Rmd` | `da_les.rds`, `da_sub.*`, `res_lesion_tbl.*`, `les_bulk_deg.*`, `subtype_bulk_deg.*`, `subtype_celltype_deg.*`, `niche_*.rds` |
| 7 | `final/fig7/fig7.Rmd` | `deg_traj.rds`, `deg_sub_traj.rds`, `deg_sub_rev.rds`, `score_df.rds`, `cox_continuous_subtype_scores.*`, `surv_*.csv` |

---

## Reproducing the analysis

1. Clone the repository and open `CTCL.Rproj` in RStudio. `.Rprofile` activates renv.
2. Restore the R environment:
   ```r
   renv::restore()
   ```
   fig4 also uses two custom helper packages, [`R2easyviz`](https://github.com/ryanreis333/R2easyviz) (downsampling) and `R2Pseudosc` (pseudobulk aggregation):
   ```r
   remotes::install_github("ryanreis333/R2easyviz")
   remotes::install_github("ryanreis333/R2Pseudosc")
   ```
3. Place the large input objects (see below) at the paths the notebooks expect.
4. Run the notebooks in order: **fig1 → fig2, fig3, fig4 (+ surv.Rmd) → fig5 → fig6 → fig7.**
   - Every later notebook reads `final/fig1/results/obj.fig1.rds`.
   - fig6 reads `final/fig5/unified.rds` and the fig4 subtype labels.
   - fig7 reads fig4 results and the files in `int/`.

All `results/` files needed to regenerate the figure panels are included. Most panels can be re-plotted directly from them without the large objects.

### Large objects (not included)

These Seurat objects exceed GitHub's file-size limits and are not in this repository.

| Object | Size | Expected path (relative to repo root) | Produced by |
|---|---|---|---|
| Full atlas counts (BPCells) + metadata + UMAP | — | `../all_cells_final_data/` | upstream preprocessing |
| `obj.fig1.rds` | 234 MB | `final/fig1/results/` | fig1 |
| `fig3_t.cells.rds` | 315 MB | `final/fig3/results/` | fig3 |
| `immune.rds`, `stromal.rds`, `unified.rds` | 1.8–4.1 GB | `final/fig5/` | fig5 |

---

## Data sources

The atlas includes 18 cohorts:

| Cohort / accession | Libraries | Biological samples | Patients | Tissue | Diagnoses |
|---|---:|---:|---:|---|---|
| E-MTAB-12303 | 40 | 8 | 8 | skin | MF |
| E-MTAB-8142 | 12 | 12 | 12 | skin | BID, HC |
| GSE122703 | 2 | 1 | 1 | blood | SS |
| GSE124898 | 2 | 1 | 1 | blood | SS |
| GSE126310 | 2 | 2 | 2 | blood | HC, SS |
| GSE128531 | 9 | 9 | 9 | skin | HC, MF, SS |
| GSE146586 | 5 | 5 | 5 | blood | SS |
| GSE165623 | 3 | 3 | 1 | blood, lymph node, skin | MF |
| GSE171811 | 14 | 14 | 8 | blood, skin | HC, MF, SS |
| GSE173205 | 27 | 27 | 14 | skin | FMF, HC, MF |
| GSE182861 | 8 | 8 | 7 | blood, skin | MF, SS |
| GSE207679 | 8 | 8 | 8 | blood, skin | MF, SS |
| GSE284075 | 40 | 26 | 17 | blood, skin | MF, SS |
| HRA000166 | 16 | 16 | 12 | skin | MF, PCALCL |
| PRJNA842142 | 15 | 15 | 6 | blood | SS |
| SRP332550 | 15 | 15 | 8 | skin | MF |
| Shandong | 4 | 4 | 4 | blood | HC, SS |
| WashU | 16 | 16 | 6 | blood | SS |

Libraries from the same specimen (e.g., sorted fractions) were merged into one biological sample (`orig.merge`). A few patients appear in more than one cohort.

Independent validation cohort: bulk RNA-seq of advanced-stage MF skin biopsies, [GSE168508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168508). The expression file and processed metadata are in `int/`.

MF, mycosis fungoides; FMF, folliculotropic MF; SS, Sézary syndrome; PCALCL, primary cutaneous anaplastic large cell lymphoma; BID, benign inflammatory dermatoses; HC, healthy control.

---

## Software

R 4.5.3 with Seurat 5.4.0, BPCells 0.3.1, harmony 1.2.4, edgeR 4.8.2, limma 3.66.0, variancePartition 1.41.4 (dream), ashr 2.2-63, fgsea 1.36.2, msigdbr 26.1.0, decoupleR 2.16.0, ConsensusClusterPlus, speckle 1.10.0, survival 3.8-6, survminer 0.5.2 and GEOquery 2.78.0. Exact versions are in `renv.lock`, and citations are in `docs/software_refs.bib`.

## Citation

If you use this code or these results, please cite the manuscript above. The citation will be updated on publication.

## License

See [LICENSE](LICENSE).
