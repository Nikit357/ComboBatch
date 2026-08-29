# Examples

## The example dataset

`example_dataset/` holds five public GEO cohorts — 144 samples across 4 platforms, spanning
FL, DLBCL and normal B cells, both fresh-frozen and FFPE. It is deliberately **ragged**:
gene coverage differs sharply between platforms, which is what makes the imputation axis
meaningful rather than decorative. See its `PROVENANCE.md` for accessions and attribution.

## Configurations

| Config | Shows |
|---|---|
| `configs/minimal.yaml` | The smallest useful run — one imputer, one method, no R required. |
| `configs/full_crossproduct.yaml` | The 4 × 12 cross-product with metrics, for the dispatcher. |
| `configs/hyperparameters.yaml` | Tuning method and imputer hyperparameters, including a sweep of one method at two values of `k`. |
| `configs/shambhala.yaml` | `20_shambhala`, with its raw-scale and parallelism constraints spelled out. |
| `configs/best_approaches/` | The source benchmark's recommended configurations, one per scenario, with the eight caveats that qualify them. |

```bash
combobatch run --config examples/configs/minimal.yaml                    # ~12 s
combobatch dispatch --config examples/configs/full_crossproduct.yaml --dry-run
```

Paths inside these files are relative, so run them from the repository root.

## A worked result

Four methods on `strict` imputation, from the example dataset (`r2` is the share of
principal-component variance explained by a factor — lower batch is better mixing, higher
biology is better preservation):

| Method | `r2_RNA_BATCH` | `r2_Diagnosis_cell_type_unified` |
|---|---|---|
| `01_raw` (baseline) | 0.297 | 0.313 |
| `02_median_scaling` | **0.038** | **0.133** |
| `15_fsqn_py` | 0.296 | 0.394 |
| `17_quantile` | 0.298 | **0.411** |

Median scaling wins decisively on batch mixing *and* destroys more than half the
biological signal in the process. That is the whole argument for benchmarking rather than
picking a method: a single batch-effect number would have ranked it first.

`marker_panel_fl.csv` is an **example** marker panel for metric groups L and M — a sample
of what `--marker-panel` expects (a `gene` column, optionally `is_housekeeping`), not a
dependency of the tool. 617 of its 633 genes are present in the example matrix.
