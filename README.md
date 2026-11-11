# Nanopore_Downstream_Analysis
These scripts are for downstream analysis for Vector Nanopore Sequencing
# Vector amplicon QC dashboard

Generates a QC report for a multiplexed amplicon panel (VGSC I–IV, ACE1 I/II, GSTe2, RDL1/2, CO1, ITS2) from per-sample coverage values.

**Outputs** (written to `--outdir`):

| File | Contents |
|---|---|
| `QC_Dashboard_<run>.html` | Interactive dashboard: per-gene and per-reaction coverage, pass rates, sample QC, negative controls, plate heatmaps |
| `QC_Summary_<run>.xlsx` | Gene / reaction / site stats, per-sample QC, failure registry, negative controls, plate coverage |
| `plots/*.png` | Static copies of every figure (skipped with `--no-images`) |

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Python 3.9+. `kaleido` is only needed for PNG export. The HTML loads plotly.js from a CDN, so viewing it needs an internet connection.

## Quick start (synthetic example)

```bash
python examples/make_example_data.py        # writes synthetic inputs (not real data)

python qc_dashboard.py \
    --coverage   examples/example_coverage.csv \
    --plate-file examples/example_plate_layout.xlsx \
    --plate-sheet plate_layout \
    --run-name   example \
    --exclude-pattern "^unassigned$" \
    --outdir     qc_output
```

## Usage

```bash
python qc_dashboard.py --coverage coverage.csv \
    --plate-file sample_sheet.xlsx --plate-sheet 1 \
    --run-name run4 --outdir qc_run4
```

| Option | Default | Meaning |
|---|---|---|
| `--coverage` | *(required)* | Long-form CSV/XLSX with columns `Sample, Gene, Coverage` |
| `--plate-file` | none | Excel workbook with the plate layout; omit to skip plate heatmaps |
| `--plate-sheet` | `1` | Sheet index (0-based) or sheet name containing the layout |
| `--outdir` | `qc_output` | Output directory |
| `--run-name` | `run` | Label used in titles and file names |
| `--threshold` | `10` | Biological samples: a measurement passes if coverage is **greater than** this |
| `--nc-threshold` | `50` | Same rule, applied to negative controls only |
| `--min-failed-targets` | `6` | A sample fails QC if at least this many of its 11 targets fail |
| `--strip-pattern` | none | Regex; text it matches is deleted from every sample name |
| `--control-prefix` | `NC` | Samples whose name starts with this are negative controls |
| `--exclude-pattern` | none | Regex; matching samples are left out of QC and listed on the `Excluded_Rows` sheet |
| `--site-pattern` / `--region-pattern` | none | Regex whose match (or first group) in the sample name defines site / region; otherwise "Unknown" |
| `--no-images` | off | Skip PNG export |

## Input formats

**Coverage table** – one row per sample × target:

```
Sample,Gene,Coverage
ALP-10-R1-H1-1-L1,VGSCI,33.7
ALP-10-R1-H1-1-L1,cox1,120.2
```

- **Sample names are used exactly as written** – nothing is assumed to come before or after them. If your table carries extra text (for example pipeline file paths such as `flagstat/<sample>_coverage_mean`), remove it with `--strip-pattern`, e.g. `--strip-pattern "^.*flagstat/|_(resistance|cox|its)_coverage_mean$"`.
- The names in the coverage table must match the names in the plate layout.
- Gene names are normalised (`cox1` → `CO1`, `its_region` → `ITS2`, `rdl1` → `RDL1`, …); see `target_map` in the script. Edit it if your panel uses other names.

**Plate layout sheet** – blocks headed `Plate_N`, followed by rows `A`–`H` and 12 columns of sample IDs:

```
Plate_1   1        2     3 ...  12
A         sampleA  ...
B         ...
```

**Site / region** grouping is off unless you give `--site-pattern` / `--region-pattern`, e.g. `--site-pattern "^[A-Za-z]+-\d+"`.

## Notes on the QC rules

- Pass/Fail uses `coverage > threshold`. For negative controls the same rule is applied at `--nc-threshold`, so a control with low coverage is labelled "Fail" – read that as "clean". Flip the logic in `add_status` if you prefer it the other way round.
- Multiplex reactions are defined in `reaction_map` / `reaction_labels` at the top of the script.

## Reproducibility

- Inputs are passed explicitly on the command line; there are no hard-coded paths.
- `examples/make_example_data.py` uses a fixed random seed, so the example files are identical every time.
- Pin exact package versions with `pip freeze > requirements.lock` when you publish results.

## Data

Do not commit real sample IDs, metadata or coverage tables unless you are cleared to share them – `.gitignore` excludes `data/` and generated reports by default.
