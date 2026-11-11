#!/usr/bin/env python3
"""
Amplicon-panel sequencing QC dashboard.

Reads a long-form coverage table (Sample, Gene, Coverage) and an optional plate-layout
workbook, and produces:
  * an interactive HTML dashboard
  * an Excel summary workbook
  * (optionally) static PNG copies of every figure

Run `python qc_dashboard.py --help` for all options.
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

# --------------------------------------------------------
# TARGET MAP
# --------------------------------------------------------
target_map = {
    'vgsci': 'VGSCI', 'vgscii': 'VGSCII', 'vgsciii': 'VGSCIII', 'vgsciv': 'VGSCIV',
    'ace1_i': 'ACE1_I',
    'vgsc1': 'VGSCI', 'vgsc2': 'VGSCII', 'vgsc3': 'VGSCIII', 'vgsc4': 'VGSCIV',
    'vgsc_1': 'VGSCI', 'vgsc_2': 'VGSCII', 'vgsc_3': 'VGSCIII', 'vgsc_4': 'VGSCIV',
    'ace1_1': 'ACE1_I', 'ace1_2': 'ACE1_II', 'ace1_i': 'ACE1_I', 'ace1_ii': 'ACE1_II',
    'gste2': 'GSTe2', 'rdl1': 'RDL1', 'rdl2': 'RDL2',
    'its': 'ITS2', 'its2': 'ITS2', 'its_region': 'ITS2', 'cox1': 'CO1'
}

# --------------------------------------------------------
# REACTION MAPPING
# --------------------------------------------------------
reaction_map = {
    'VGSCII': 1, 'VGSCIII': 1, 'VGSCIV': 1,
    'VGSCI': 2, 'ITS2': 2,
    'CO1': 3, 'ACE1_II': 3, 'RDL2': 3,
    'GSTe2': 4, 'RDL1': 4,
    'ACE1_I': 5
}

reaction_labels = {
    1: "Reaction 1 (VGSCII, VGSCIII, VGSCIV)",
    2: "Reaction 2 (VGSCI, ITS2)",
    3: "Reaction 3 (CO1, ACE1_II, RDL2)",
    4: "Reaction 4 (GSTe2, RDL1)",
    5: "Reaction 5 (ACE1_I)"
}

# --------------------------------------------------------
# GENE PALETTE
# --------------------------------------------------------
gene_palette = {
    "VGSCI": "#1f77b4", "VGSCII": "#4c78a8", "VGSCIII": "#6baed6", "VGSCIV": "#9ecae1",
    "ACE1_I": "#2ca02c", "ACE1_II": "#98df8a", "GSTe2": "#ff7f0e", "RDL1": "#d62728",
    "RDL2": "#9467bd", "ITS2": "#e377c2", "CO1": "#8c564b"
}

PLATE_ROWS = list("ABCDEFGH")
PLATE_COLS = list(range(1, 13))

# --------------------------------------------------------
# HELPERS
# --------------------------------------------------------
def prepare_dataframe(df):
    if df.empty:
        return df
    df = df.copy()
    df["Gene"] = df["Gene"].astype(str).str.lower().map(target_map).fillna(df["Gene"])
    df["Sample"] = df["Sample"].astype(str).str.strip()
    return df

def clean_sample_id(series, strip_pattern=None):
    """Sample names are used exactly as given (only surrounding whitespace is removed).
    If the names carry extra text, pass a regex via strip_pattern; whatever it matches is deleted,
    e.g. r"^.*flagstat/|_(resistance|cox|its)_coverage_mean$"."""
    out = series.astype(str).str.strip()
    if strip_pattern:
        out = out.str.replace(strip_pattern, "", regex=True).str.strip()
    return out

def load_coverage(filepath, strip_pattern=None):
    # Supports the long-form CSV (or an Excel workbook, as in earlier runs)
    if str(filepath).lower().endswith((".xlsx", ".xls")):
        xl = pd.ExcelFile(filepath)
        frames = [xl.parse(sh) for sh in xl.sheet_names]
        frames = [f for f in frames if all(c in f.columns for c in ["Sample", "Gene", "Coverage"])]
        if not frames:
            raise ValueError("No valid data sheets.")
        combined = pd.concat(frames, ignore_index=True)
    else:
        combined = pd.read_csv(filepath)
        missing = [c for c in ["Sample", "Gene", "Coverage"] if c not in combined.columns]
        if missing:
            raise ValueError(f"Missing columns in {filepath}: {missing}")

    combined["Sample"] = clean_sample_id(combined["Sample"], strip_pattern)
    combined = prepare_dataframe(combined)
    combined["Coverage"] = pd.to_numeric(combined["Coverage"], errors="coerce")
    combined = combined.dropna(subset=["Coverage"])
    return combined

def compute_stats(df, group_cols, threshold=10):
    stats = df.groupby(group_cols).agg(
        Total=("Coverage", "size"),
        Mean=("Coverage", "mean"),
        Median=("Coverage", "median"),
        Passed=("Coverage", lambda x: (x > threshold).sum()),
        Failed=("Coverage", lambda x: (x <= threshold).sum())
    ).reset_index()
    stats["SuccessRate"] = (stats["Passed"] / stats["Total"] * 100).round(1)
    return stats

def add_status(df, threshold=10, nc_threshold=50):
    """Bio samples use `threshold`; negative controls (Is_Control) use `nc_threshold`."""
    cutoff = np.where(df["Is_Control"], nc_threshold, threshold)
    df["Status"] = np.where(df["Coverage"] > cutoff, "Pass", "Fail")
    return df

def load_plate_layout(filepath, sheet=1):
    """Parse the plate-layout sheet: blocks headed 'Plate_N' followed by rows A-H x columns 1-12."""
    raw = pd.read_excel(filepath, sheet_name=sheet, header=None)
    records, plate = [], None
    for i in range(len(raw)):
        first = raw.iloc[i, 0]
        if isinstance(first, str) and re.match(r"^Plate[_ ]?\d+", first.strip()):
            plate = int(re.findall(r"\d+", first)[0])
            continue
        if plate is not None and isinstance(first, str) and first.strip().upper() in PLATE_ROWS:
            for col in range(1, 13):
                val = raw.iloc[i, col] if col < raw.shape[1] else np.nan
                if pd.notna(val) and str(val).strip():
                    records.append({"Plate": plate, "Row": first.strip().upper(),
                                    "Column": col, "Sample": str(val).strip()})
    return pd.DataFrame(records)

# --------------------------------------------------------
# MAIN
# --------------------------------------------------------
def main(args):
    # Settings (from the command line; see parse_args)
    THRESHOLD = args.threshold
    NC_THRESHOLD = args.nc_threshold
    MIN_FAILED_TARGETS = args.min_failed_targets
    COVERAGE_FILE = args.coverage
    PLATE_FILE = args.plate_file
    PLATE_SHEET = args.plate_sheet
    SAVE_PLOTS = not args.no_images
    RUN_NAME = args.run_name
    OUTDIR = args.outdir
    PLOTS_DIR = os.path.join(OUTDIR, "plots")
    HTML_OUTPUT = os.path.join(OUTDIR, f"QC_Dashboard_{RUN_NAME}.html")
    EXCEL_OUTPUT = os.path.join(OUTDIR, f"QC_Summary_{RUN_NAME}.xlsx")
    os.makedirs(OUTDIR, exist_ok=True)

    if SAVE_PLOTS:
        os.makedirs(PLOTS_DIR, exist_ok=True)

    print("Loading coverage data...")
    coverage_df = load_coverage(COVERAGE_FILE, args.strip_pattern)
    print(f"Loaded {len(coverage_df)} records from {coverage_df['Sample'].nunique()} samples.")

    # Optional grouping: derive Study_Site / Region from the sample name only if a regex is given
    # (first capture group, or the whole match if there is no group). Otherwise everything is "Unknown".
    def extract_group(pattern):
        if not pattern:
            return pd.Series("Unknown", index=coverage_df.index)
        rx = pattern if re.compile(pattern).groups else f"({pattern})"
        return coverage_df["Sample"].str.extract(rx)[0].fillna("Unknown")
    coverage_df["Study_Site"] = extract_group(args.site_pattern)
    coverage_df["Region"] = extract_group(args.region_pattern)

    # Negative controls: sample names starting with --control-prefix (case-insensitive)
    coverage_df["Is_Control"] = coverage_df["Sample"].str.upper().str.startswith(args.control_prefix.upper())
    coverage_df.loc[coverage_df["Is_Control"], ["Study_Site", "Region"]] = "Negative Control"

    # Rows to leave out of QC entirely (e.g. reads that did not demultiplex); reported on their own sheet
    if args.exclude_pattern:
        coverage_df["Is_Unassigned"] = coverage_df["Sample"].str.contains(args.exclude_pattern, regex=True, case=False)
    else:
        coverage_df["Is_Unassigned"] = False
    unassigned_df = coverage_df[coverage_df["Is_Unassigned"]]
    coverage_df = coverage_df[~coverage_df["Is_Unassigned"]].copy()

    # Add reaction info and status
    coverage_df["Reaction"] = coverage_df["Gene"].map(reaction_map)
    coverage_df["Reaction_Label"] = coverage_df["Reaction"].map(reaction_labels)
    coverage_df = add_status(coverage_df, THRESHOLD, NC_THRESHOLD)

    # Enforce strict chronological order for reactions
    ordered_reaction_labels = [
        reaction_labels[1],
        reaction_labels[2],
        reaction_labels[3],
        reaction_labels[4],
        reaction_labels[5]
    ]
    coverage_df["Reaction_Label"] = pd.Categorical(
        coverage_df["Reaction_Label"], 
        categories=ordered_reaction_labels, 
        ordered=True
    )

    # Separate biological samples from negative controls for general stats
    bio_df = coverage_df[~coverage_df["Is_Control"]]
    nc_df = coverage_df[coverage_df["Is_Control"]]

    # --------------------------------------------------------
    # STATISTICS
    # --------------------------------------------------------
    gene_stats = compute_stats(bio_df, ["Gene"], THRESHOLD)
    reaction_stats = compute_stats(bio_df, ["Reaction_Label"], THRESHOLD)
    site_stats = compute_stats(bio_df, ["Study_Site"], THRESHOLD)
    region_stats = compute_stats(bio_df, ["Region"], THRESHOLD)

    # Sample-level stats
    sample_stats = bio_df.groupby("Sample").agg(
        Total_Targets=("Gene", "nunique"),
        Total_Measurements=("Coverage", "size"),
        Passed=("Status", lambda x: (x == "Pass").sum()),
        Failed=("Status", lambda x: (x == "Fail").sum())
    ).reset_index()
    sample_stats["SuccessRate"] = (sample_stats["Passed"] / sample_stats["Total_Measurements"] * 100).round(1)
    sample_stats["QC_Status"] = np.where(sample_stats["Failed"] >= MIN_FAILED_TARGETS, "Fail", "Pass")
    sample_stats = sample_stats.sort_values("Failed", ascending=False)

    # Negative Control Summary
    nc_summary = pd.DataFrame()
    nc_bar_fig = None
    nc_box_fig = None
    if not nc_df.empty:
        nc_summary = nc_df[["Sample", "Gene", "Coverage", "Status"]].sort_values(["Sample", "Gene"])
        
        nc_bar_fig = px.bar(nc_df, x='Sample', y='Coverage', color='Gene',
                            barmode='group', title="Negative Control Coverage per Target Gene",
                            labels={'Coverage': 'Coverage (x)', 'Sample': 'Negative Control ID'},
                            color_discrete_map=gene_palette)
        nc_bar_fig.add_hline(y=NC_THRESHOLD, line_dash="dash", line_color="red", annotation_text=f"NC threshold ({NC_THRESHOLD}x)")
        
        nc_box_fig = px.box(nc_df, x='Gene', y='Coverage', color='Gene',
                            title="Negative Control Coverage Distribution across Genes",
                            color_discrete_map=gene_palette)
        nc_box_fig.add_hline(y=NC_THRESHOLD, line_dash="dash", line_color="red", annotation_text=f"NC threshold ({NC_THRESHOLD}x)")
        nc_box_fig.update_layout(showlegend=False)

    failure_registry = bio_df[bio_df["Status"] == "Fail"][
        ["Sample", "Gene", "Coverage", "Study_Site", "Region"]
    ].sort_values(["Sample", "Gene"])

    # Completely failed biological samples
    all_targets = ['VGSCI', 'VGSCII', 'VGSCIII', 'VGSCIV', 'ACE1_I', 'ACE1_II',
                   'GSTe2', 'RDL1', 'RDL2', 'CO1', 'ITS2']
    sample_gene_sets = bio_df.groupby("Sample")["Gene"].apply(set).reset_index(name="gene_set")
    sample_gene_sets["has_all"] = sample_gene_sets["gene_set"].apply(lambda s: all(t in s for t in all_targets))
    complete_samples = sample_gene_sets[sample_gene_sets["has_all"]]["Sample"].tolist()
    completely_failed = bio_df[
        (bio_df["Sample"].isin(complete_samples)) &
        (bio_df["Status"] == "Fail")
    ].groupby("Sample")["Gene"].count().reset_index(name="Failed_Count")
    completely_failed = completely_failed[completely_failed["Failed_Count"] == len(all_targets)]
    completely_failed_samples = completely_failed["Sample"].tolist()

    # --------------------------------------------------------
    # PLOTLY FIGURES
    # --------------------------------------------------------
    print("Generating figures...")

    mean_bar_fig = px.bar(gene_stats, x='Gene', y='Mean', color='Mean',
                          title="Mean Coverage per Gene (Biological Samples)", labels={'Mean': 'Mean Coverage (x)'},
                          color_continuous_scale='Blues')
    mean_bar_fig.update_layout(showlegend=False)

    median_bar_fig = px.bar(gene_stats, x='Gene', y='Median', color='Median',
                            title="Median Coverage per Gene (Biological Samples)", labels={'Median': 'Median Coverage (x)'},
                            color_continuous_scale='Greens')
    median_bar_fig.update_layout(showlegend=False)

    success_bar_fig = px.bar(gene_stats, x='Gene', y='SuccessRate', text='SuccessRate',
                             title=f"Success Rate per Gene (coverage > {THRESHOLD}x)",
                             labels={'SuccessRate': 'Success Rate (%)'},
                             color='SuccessRate', color_continuous_scale='RdYlGn')
    success_bar_fig.update_traces(texttemplate='%{text:.1f}%', textposition='outside')
    success_bar_fig.update_layout(showlegend=False)

    fail_count_fig = px.bar(gene_stats, x='Gene', y='Failed', title="Failed Measurements per Gene",
                            labels={'Failed': 'Failed Measurements'},
                            color='Failed', color_continuous_scale='Reds')
    fail_count_fig.update_layout(showlegend=False)

    box_fig = px.box(bio_df, x='Gene', y='Coverage', color='Gene',
                     log_y=True, points='outliers',
                     color_discrete_map=gene_palette,
                     title="Coverage Distribution per Gene (Biological Samples)")
    box_fig.add_hline(y=THRESHOLD, line_dash="dash", line_color="red", annotation_text=f"Threshold ({THRESHOLD}x)")
    box_fig.update_yaxes(type="log", range=[0, 5.3], dtick=1)
    box_fig.update_layout(showlegend=False)

    violin_fig = px.violin(bio_df, x='Gene', y='Coverage', color='Gene',
                           log_y=True, box=True, points='all',
                           color_discrete_map=gene_palette,
                           title="Violin + Box plots per Gene")
    violin_fig.add_hline(y=THRESHOLD, line_dash="dash", line_color="red")
    violin_fig.update_yaxes(type="log", range=[0, 5.3], dtick=1)
    violin_fig.update_layout(showlegend=False)

    reaction_box_fig = px.box(bio_df, x='Reaction_Label', y='Coverage', color='Reaction_Label',
                              log_y=True, points='all',
                              title="Coverage per Multiplex Reaction",
                              labels={'Coverage': 'Coverage (x)', 'Reaction_Label': 'Reaction'},
                              category_orders={'Reaction_Label': ordered_reaction_labels})
    reaction_box_fig.add_hline(y=THRESHOLD, line_dash="dash", line_color="red", annotation_text=f"Threshold ({THRESHOLD}x)")
    reaction_box_fig.update_yaxes(type="log", range=[0, 5.3], dtick=1)
    reaction_box_fig.update_layout(showlegend=False)

    reaction_mean_fig = px.bar(reaction_stats, x='Reaction_Label', y='Mean', text='Mean',
                               title="Mean Coverage per Reaction",
                               labels={'Mean': 'Mean Coverage (x)'},
                               color='Mean', color_continuous_scale='Purples',
                               category_orders={'Reaction_Label': ordered_reaction_labels})
    reaction_mean_fig.update_traces(texttemplate='%{text:.1f}x', textposition='outside')
    reaction_mean_fig.update_layout(showlegend=False)

    sample_stack = sample_stats.melt(
        id_vars=['Sample'],
        value_vars=['Passed', 'Failed'],
        var_name='Status',
        value_name='Count'
    )
    sample_qc_fig = px.bar(
        sample_stack,
        x='Sample',
        y='Count',
        color='Status',
        title="Target Status per Sample (Passed vs Failed) – sorted by failed count descending",
        labels={'Count': 'Number of Targets', 'Sample': 'Sample ID'},
        color_discrete_map={'Passed': '#2ca02c', 'Failed': '#d62728'},
        barmode='stack',
        category_orders={'Sample': sample_stats['Sample'].tolist()}
    )
    sample_qc_fig.update_layout(
        height=450,
        width=2500,
        xaxis={
            'categoryorder': 'array',
            'categoryarray': sample_stats['Sample'].tolist(),
            'tickangle': -45,
            'automargin': True,
            'title': 'Sample ID'
        },
        yaxis_title='Number of Targets',
        legend_title='Status',
        margin=dict(b=150)
    )

    cf_df = pd.DataFrame({'Sample': completely_failed_samples, 'Count': [len(all_targets)] * len(completely_failed_samples)})
    cf_fig = None
    if not cf_df.empty:
        cf_fig = px.bar(cf_df, x='Sample', y='Count',
                        title=f"Samples with ALL {len(all_targets)} targets failing (≤{THRESHOLD}x)",
                        labels={'Count': 'Number of Failed Targets'})
        cf_fig.update_traces(marker_color='red')

    # Dummy Placeholders for regional/site analysis since metadata is excluded
    site_perf_fig = px.bar(site_stats, x='Study_Site', y='SuccessRate', text='SuccessRate',
                           title="Success Rate by Study Site (from sample ID prefix)",
                           labels={'SuccessRate': 'Success Rate (%)'})
    site_box_fig = px.box(bio_df, x='Study_Site', y='Coverage', color='Study_Site',
                          log_y=True, title="Coverage Distribution by Site Group")
    region_perf_fig = px.bar(region_stats, x='Region', y='SuccessRate', text='SuccessRate',
                             title="Success Rate by Site Code (from sample ID prefix)",
                             labels={'SuccessRate': 'Success Rate (%)'})
    region_box_fig = px.box(bio_df, x='Region', y='Coverage', color='Region',
                            log_y=True, title="Coverage Distribution by Region")

    # --------------------------------------------------------
    # PLATE HEATMAPS (per plate, all wells incl. NCs)
    # --------------------------------------------------------
    plate_df = pd.DataFrame()
    plate_figs = []          # (fig, name) for saving
    plate_summary = pd.DataFrame()
    if PLATE_FILE and os.path.exists(PLATE_FILE):
        plate_df = load_plate_layout(PLATE_FILE, PLATE_SHEET)
        if plate_df.empty:
            print("  Warning: no plate blocks (rows headed 'Plate_N', then rows A-H) found in the plate sheet.")
        print(f"Plate layout parsed: {len(plate_df)} wells across {plate_df['Plate'].nunique()} plates.")
        missing = sorted(set(coverage_df["Sample"]) - set(plate_df["Sample"]))
        if missing:
            print(f"  Note: {len(missing)} samples in coverage not found on any plate: {missing[:10]}")
    else:
        print("No plate file supplied/found; skipping plate heatmaps.")

    if not plate_df.empty:
        per_sample = coverage_df.groupby("Sample").agg(
            Avg_Coverage=("Coverage", "mean"),
            Failed_Targets=("Status", lambda x: (x == "Fail").sum())
        ).reset_index()
        plate_merged = plate_df.merge(per_sample, on="Sample", how="left")
        plate_summary = plate_merged.sort_values(["Plate", "Row", "Column"])

        for plate_num in sorted(plate_merged["Plate"].unique()):
            pdata = plate_merged[plate_merged["Plate"] == plate_num].set_index(["Row", "Column"])
            z_log, z_fail, label, hover = [], [], [], []
            for r in PLATE_ROWS:
                zl, zf, lb, hv = [], [], [], []
                for c in PLATE_COLS:
                    if (r, c) in pdata.index:
                        rec = pdata.loc[(r, c)]
                        avg, nf, name = rec["Avg_Coverage"], rec["Failed_Targets"], rec["Sample"]
                        is_nc = str(name).upper().startswith(args.control_prefix.upper())
                        zl.append(np.log10(avg + 1) if pd.notna(avg) else np.nan)
                        zf.append(nf if pd.notna(nf) else np.nan)
                        lb.append("NC" if is_nc else "")
                        cutoff = NC_THRESHOLD if is_nc else THRESHOLD
                        hv.append(f"{name}<br>Well {r}{c}<br>Avg coverage: {avg:,.1f}x<br>"
                                  f"Targets failed (≤{cutoff}x): {int(nf) if pd.notna(nf) else 'NA'}/11")
                    else:
                        zl.append(np.nan); zf.append(np.nan); lb.append(""); hv.append("Empty")
                z_log.append(zl); z_fail.append(zf); label.append(lb); hover.append(hv)

            tick_vals = [0, 1, 2, 3, 4, 5]
            cov_fig = go.Figure(go.Heatmap(
                z=z_log, x=PLATE_COLS, y=PLATE_ROWS, colorscale="Viridis",
                text=label, texttemplate="%{text}", textfont=dict(color="red", size=14, family="Arial Black"),
                customdata=hover, hovertemplate="%{customdata}<extra></extra>",
                colorbar=dict(title="Avg coverage (x)", tickvals=tick_vals,
                              ticktext=[f"{10**v - 1:,.0f}" for v in tick_vals]),
                zmin=0, zmax=5, xgap=2, ygap=2))
            cov_fig.update_yaxes(autorange="reversed", title="Row")
            cov_fig.update_xaxes(dtick=1, title="Column", side="top")
            cov_fig.update_layout(title=f"Plate {plate_num} – Average Coverage per Well (log colour scale; NC wells labelled)",
                                  height=500)
            plate_figs.append((cov_fig, f"plate_{plate_num}_avg_coverage"))

            fail_fig = go.Figure(go.Heatmap(
                z=z_fail, x=PLATE_COLS, y=PLATE_ROWS, colorscale="Reds", zmin=0, zmax=11,
                text=label, texttemplate="%{text}", textfont=dict(color="blue", size=14, family="Arial Black"),
                customdata=hover, hovertemplate="%{customdata}<extra></extra>",
                colorbar=dict(title="Failed targets"), xgap=2, ygap=2))
            fail_fig.update_yaxes(autorange="reversed", title="Row")
            fail_fig.update_xaxes(dtick=1, title="Column", side="top")
            fail_fig.update_layout(title=f"Plate {plate_num} – Failed Targets per Well (of 11; bio ≤{THRESHOLD}x, NC ≤{NC_THRESHOLD}x)",
                                   height=500)
            plate_figs.append((fail_fig, f"plate_{plate_num}_failed_targets"))

    # --------------------------------------------------------
    # SAVE PLOTS
    # --------------------------------------------------------
    if SAVE_PLOTS:
        all_figs = [
            (mean_bar_fig, "mean_coverage_per_gene"),
            (median_bar_fig, "median_coverage_per_gene"),
            (success_bar_fig, "success_rate_per_gene"),
            (fail_count_fig, "failed_count_per_gene"),
            (box_fig, "coverage_boxplots"),
            (violin_fig, "coverage_violin"),
            (reaction_box_fig, "reaction_boxplot"),
            (reaction_mean_fig, "reaction_mean_bar"),
            (sample_qc_fig, "sample_qc_stacked"),
            (cf_fig, "completely_failed_samples"),
            (site_perf_fig, "site_performance"),
            (site_box_fig, "site_coverage_box"),
            (region_perf_fig, "region_performance"),
            (region_box_fig, "region_coverage_box"),
            (nc_bar_fig, "negative_control_bar"),
            (nc_box_fig, "negative_control_box"),
        ]
        all_figs.extend(plate_figs)
        for fig, name in all_figs:
            if fig is not None:
                try:
                    pio.write_image(fig, os.path.join(PLOTS_DIR, f"{name}.png"), width=2000, height=800, scale=2)
                    print(f"  Saved: {name}.png")
                except Exception as e:
                    print(f"  Could not save {name}: {e}")

    # --------------------------------------------------------
    # EXPORT EXCEL
    # --------------------------------------------------------
    with pd.ExcelWriter(EXCEL_OUTPUT, engine='openpyxl') as writer:
        gene_stats.to_excel(writer, sheet_name='Gene_Stats', index=False)
        reaction_stats.to_excel(writer, sheet_name='Reaction_Stats', index=False)
        site_stats.to_excel(writer, sheet_name='Site_Stats', index=False)
        region_stats.to_excel(writer, sheet_name='Region_Stats', index=False)
        sample_stats.to_excel(writer, sheet_name='Sample_QC', index=False)
        if not failure_registry.empty:
            failure_registry.to_excel(writer, sheet_name='Failure_Registry', index=False)
        if not nc_summary.empty:
            nc_summary.to_excel(writer, sheet_name='Negative_Controls', index=False)
        if not plate_summary.empty:
            plate_summary.to_excel(writer, sheet_name='Plate_Coverage', index=False)
        if not unassigned_df.empty:
            unassigned_df[['Sample', 'Gene', 'Coverage']].to_excel(writer, sheet_name='Excluded_Rows', index=False)
        if completely_failed_samples:
            pd.DataFrame({'Completely_Failed_Samples': completely_failed_samples}).to_excel(
                writer, sheet_name='Completely_Failed', index=False)

    print(f"Excel saved: {EXCEL_OUTPUT}")

    # --------------------------------------------------------
    # HTML DASHBOARD
    # --------------------------------------------------------
    def fig_to_html(fig, div_id):
        if fig is None:
            return "<p class='text-muted'>Plot not available.</p>"
        return fig.to_html(full_html=False, include_plotlyjs='cdn', div_id=div_id)

    plate_html = "".join(
        f'<div class="row"><div class="col-lg-6">{fig_to_html(f, name)}</div>'
        f'<div class="col-lg-6">{fig_to_html(g, name2)}</div></div>'
        for (f, name), (g, name2) in zip(plate_figs[0::2], plate_figs[1::2])
    ) or "<p class='text-muted'>Plate information not available.</p>"

    html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{RUN_NAME} – QC Dashboard</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{ background: #f8f9fa; padding: 20px; }}
        .card {{ margin-bottom: 30px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .card-header {{ background-color: #e9ecef; font-weight: bold; }}
        .table-container {{ overflow-x: auto; }}
        .plot-container {{ overflow-x: auto; }}
    </style>
</head>
<body>
<div class="container-fluid">
    <h1 class="text-center mb-4">{RUN_NAME} – QC Dashboard</h1>

    <!-- Overall stats -->
    <div class="row">
        <div class="col-12">
            <div class="card">
                <div class="card-header">Overall Statistics (Biological Samples)</div>
                <div class="card-body">
                    <p><strong>Total Biological Samples:</strong> {bio_df['Sample'].nunique()}</p>
                    <p><strong>Total Negative Control Samples:</strong> {nc_df['Sample'].nunique()}</p>
                    <p><strong>Total Targets:</strong> {bio_df['Gene'].nunique()}</p>
                    <p><strong>Pass Rate (Biological):</strong> {(bio_df['Status'] == 'Pass').mean() * 100:.1f}%</p>
                </div>
            </div>
        </div>
    </div>

    <!-- Gene-level -->
    <div class="row">
        <div class="col-md-4"><div class="card"><div class="card-header">Mean Coverage</div><div class="card-body">{fig_to_html(mean_bar_fig, 'mean_bar')}</div></div></div>
        <div class="col-md-4"><div class="card"><div class="card-header">Median Coverage</div><div class="card-body">{fig_to_html(median_bar_fig, 'median_bar')}</div></div></div>
        <div class="col-md-4"><div class="card"><div class="card-header">Success Rate</div><div class="card-body">{fig_to_html(success_bar_fig, 'success_bar')}</div></div></div>
    </div>
    <div class="row">
        <div class="col-md-6"><div class="card"><div class="card-header">Failed Measurements</div><div class="card-body">{fig_to_html(fail_count_fig, 'fail_count')}</div></div></div>
        <div class="col-md-6"><div class="card"><div class="card-header">Sample QC (Stacked)</div><div class="card-body"><div class="plot-container">{fig_to_html(sample_qc_fig, 'sample_qc')}</div></div></div></div>
    </div>

    <!-- Distributions -->
    <div class="row">
        <div class="col-md-6"><div class="card"><div class="card-header">Boxplots</div><div class="card-body">{fig_to_html(box_fig, 'box')}</div></div></div>
        <div class="col-md-6"><div class="card"><div class="card-header">Violin Plots</div><div class="card-body">{fig_to_html(violin_fig, 'violin')}</div></div></div>
    </div>

    <!-- Reaction -->
    <div class="row">
        <div class="col-12"><div class="card"><div class="card-header">Coverage by Multiplex Reaction</div><div class="card-body">
            <div class="row">
                <div class="col-md-6">{fig_to_html(reaction_box_fig, 'reaction_box')}</div>
                <div class="col-md-6">{fig_to_html(reaction_mean_fig, 'reaction_mean')}</div>
            </div>
        </div></div></div>
    </div>

    <!-- Completely failed -->
    <div class="row">
        <div class="col-12"><div class="card border-danger">
            <div class="card-header bg-danger text-white">Completely Failed Samples (ALL {len(all_targets)} targets ≤ {THRESHOLD}x)</div>
            <div class="card-body">
                {fig_to_html(cf_fig, 'cf')}
                {f"<p><strong>Count:</strong> {len(completely_failed_samples)}</p>" if completely_failed_samples else "<p>None found.</p>"}
                {f"<p>Samples: {', '.join(completely_failed_samples)}</p>" if completely_failed_samples else ""}
            </div>
        </div></div>
    </div>

    <!-- Negative Controls Quality Section -->
    <div class="row">
        <div class="col-12"><div class="card border-warning">
            <div class="card-header bg-warning text-dark">Negative Controls Quality Assessment (NC Samples, cutoff {NC_THRESHOLD}x)</div>
            <div class="card-body">
                <div class="row">
                    <div class="col-md-6">{fig_to_html(nc_bar_fig, 'nc_bar')}</div>
                    <div class="col-md-6">{fig_to_html(nc_box_fig, 'nc_box')}</div>
                </div>
            </div>
        </div></div>
    </div>

    <!-- Plate heatmaps -->
    <div class="row">
        <div class="col-12"><div class="card border-primary">
            <div class="card-header bg-primary text-white">Plate Layout – Average Coverage and Failed Targets per Well (NC wells labelled)</div>
            <div class="card-body">
                {plate_html}
            </div>
        </div></div>
    </div>

    <!-- Tables -->
    <div class="row">
        <div class="col-md-6"><div class="card"><div class="card-header">Gene Stats</div><div class="card-body table-container">{gene_stats.to_html(index=False, classes="table table-sm table-striped", float_format="%.1f")}</div></div></div>
        <div class="col-md-6"><div class="card"><div class="card-header">Reaction Stats</div><div class="card-body table-container">{reaction_stats.to_html(index=False, classes="table table-sm table-striped", float_format="%.1f")}</div></div></div>
    </div>
    
    <!-- Negative Controls Table & Failure Registry -->
    <div class="row">
        <div class="col-md-6"><div class="card border-warning">
            <div class="card-header bg-warning text-dark">Negative Controls Inspection Data</div>
            <div class="card-body table-container">
                {nc_summary.to_html(index=False, classes="table table-sm table-striped") if not nc_summary.empty else "<p>No negative controls found.</p>"}
            </div>
        </div></div>
        <div class="col-md-6"><div class="card"><div class="card-header">Failure Registry (Biological Samples Only)</div><div class="card-body table-container">{failure_registry.to_html(index=False, classes="table table-sm table-danger table-striped") if not failure_registry.empty else "<p>No failures.</p>"}</div></div></div>
    </div>

</div>
</body>
</html>
"""

    with open(HTML_OUTPUT, "w") as f:
        f.write(html)

    print(f"HTML dashboard saved: {HTML_OUTPUT}")
    print("Done.")

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate a QC dashboard (HTML + Excel + PNGs) from a long-form coverage table.")
    p.add_argument("--coverage", required=True,
                   help="CSV (or XLSX) with columns Sample, Gene, Coverage (long form)")
    p.add_argument("--plate-file", default=None,
                   help="Optional Excel workbook containing the plate layout")
    p.add_argument("--plate-sheet", default="1",
                   help="Sheet holding the plate layout: 0-based index or sheet name (default: 1)")
    p.add_argument("--outdir", default="qc_output", help="Output directory (default: qc_output)")
    p.add_argument("--run-name", default="run", help="Label used in titles and output filenames")
    p.add_argument("--threshold", type=float, default=10,
                   help="Coverage cutoff for biological samples; Pass means > cutoff (default: 10)")
    p.add_argument("--nc-threshold", type=float, default=50,
                   help="Coverage cutoff for negative controls (default: 50)")
    p.add_argument("--min-failed-targets", type=int, default=6,
                   help="A sample fails QC when at least this many targets fail (default: 6)")
    p.add_argument("--strip-pattern", default=None,
                   help="Optional regex; text it matches is removed from every sample name "
                        "(default: names are used exactly as given)")
    p.add_argument("--control-prefix", default="NC",
                   help="Samples whose name starts with this are negative controls (default: NC)")
    p.add_argument("--exclude-pattern", default=None,
                   help="Optional regex; samples whose name matches are excluded from QC "
                        "and listed on the Excluded_Rows sheet (e.g. '^unassigned$')")
    p.add_argument("--site-pattern", default=None,
                   help="Optional regex whose match (or first group) in the sample name defines the study site")
    p.add_argument("--region-pattern", default=None,
                   help="Optional regex whose match (or first group) in the sample name defines the region")
    p.add_argument("--no-images", action="store_true",
                   help="Skip writing PNG copies of the figures (needs the kaleido package)")
    args = p.parse_args(argv)
    args.plate_sheet = int(args.plate_sheet) if str(args.plate_sheet).isdigit() else args.plate_sheet
    return args

if __name__ == "__main__":
    main(parse_args())
