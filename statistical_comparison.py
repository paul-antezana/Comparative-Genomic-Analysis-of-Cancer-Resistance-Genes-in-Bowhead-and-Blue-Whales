"""
statistical_comparison.py ---------------------------------------------------

Author: Paul Antezana

This script reads the comparison CSV and:
    1. Prints a formatted table (cetaceans vs suids, cancer genes vs control genes)
    2. Runs a t-test to check if whale cancer genes are significantly more divergent
       than control (housekeeping) genes
    3. Generates a visual table figure saved as a PNG

### BEFORE RUNNING:#############################
  run fetch_and_compare.py
  
  install necessary modules
################################################

"""

# import necessary modules
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# build file paths relative to this script's location
script_dir = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(script_dir, "../results/comparison.csv")
table_figure = os.path.join(script_dir, "../results/figures/comparison_table.png")

# create the figures folder if it doesn't exist yet
os.makedirs(os.path.join(script_dir, "../results/figures"), exist_ok=True)

# these lists will hold the identity values for the 4 cells of the 2x2 table
cancer_cetacean = []
cancer_suid = []
control_cetacean = []
control_suid = []

# these lists keep track of gene names for display
cancer_genes = []
control_genes = []

# for the interaction t-test: per-gene differences (suid identity - cetacean identity)
cancer_diffs = []
control_diffs = []
cancer_diff_genes = []
control_diff_genes = []

# store all rows for building the visual table
all_rows = []

# if true, it adds a title on top of the table graph
TITLE_ALLOWED = False

# read the CSV
with open(csv_file, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        all_rows.append(row)
        gene_type = row["Gene_Type"]
        cetacean = row["Cetacean_Identity"]
        suid = row["Suid_Identity"]

        # convert to data into a float, or None if data is missing
        if cetacean != "N/A":
            cetacean_val = float(cetacean)
        else:
            cetacean_val = None
        if suid != "N/A":
            suid_val = float(suid)
        else:
            suid_val = None

        # sort each value into the correct group
        if gene_type == "cancer":
            cancer_genes.append(row["Gene"])
            if cetacean_val is not None:
                cancer_cetacean.append(cetacean_val)
            if suid_val is not None:
                cancer_suid.append(suid_val)
        else:
            control_genes.append(row["Gene"])
            if cetacean_val is not None:
                control_cetacean.append(cetacean_val)
            if suid_val is not None:
                control_suid.append(suid_val)

        # if both species pairs have data for this gene, compute the difference (suid_identity - cetacean_identity), and add to list
        # positive results means whales are more diverged, negative means suids are more diverged
        if cetacean_val is not None and suid_val is not None:
            diff = suid_val - cetacean_val
            if gene_type == "cancer":
                cancer_diffs.append(diff)
                cancer_diff_genes.append(row["Gene"])
            else:
                control_diffs.append(diff)
                control_diff_genes.append(row["Gene"])

### Console result table creation ###
print()
print("#" * 70)
print()
print("=" * 70)
print("        COMPARISON TABLE: Pairwise Sequence Identity (%)")
print("=" * 70)
print(f"{'':15} {'Cetaceans':>20}   {'Suids':>20}")
print(f"{'':15} {'(Blue vs Bowhead)':>20}   {'(Boar vs Warthog)':>20}")
print("-" * 70)

# print cancer genes section
print("\nCANCER GENES:")
for row in all_rows:
    if row["Gene_Type"] == "cancer":
        cetacean_str = row["Cetacean_Identity"]
        suids_str = row["Suid_Identity"]
        if cetacean_str != "N/A":
            cetacean_str = cetacean_str + "%"
        if suids_str != "N/A":
            suids_str = suids_str + "%"
        print(f"  {row['Gene']:<13} {cetacean_str:>20}   {suids_str:>20}")

# print control genes section
print("\nCONTROL GENES:")
for row in all_rows:
    if row["Gene_Type"] == "control":
        cetacean_str = row["Cetacean_Identity"]
        suids_str = row["Suid_Identity"]
        if cetacean_str != "N/A":
            cetacean_str = cetacean_str + "%"
        if suids_str != "N/A":
            suids_str = suids_str + "%"
        print(f"  {row['Gene']:<13} {cetacean_str}   {suids_str:>20}")

print("-" * 70)

### Compute and print means for each cell ###

# cancer means
if cancer_cetacean:
    cancer_cet_mean = sum(cancer_cetacean) / len(cancer_cetacean)
    print(f"  {'Mean (cancer)':<13} {cancer_cet_mean:>19.2f}%", end="")
else:
    print(f"  {'Mean (cancer)':<13} {'N/A':>20}", end="")

if cancer_suid:
    cancer_suid_mean = sum(cancer_suid) / len(cancer_suid)
    print(f"   {cancer_suid_mean:>19.2f}%")
else:
    print(f"   {'N/A':>20}")

# control means
if control_cetacean:
    control_cet_mean = sum(control_cetacean) / len(control_cetacean)
    print(f"  {'Mean (control)':<13} {control_cet_mean:>19.2f}%", end="")
else:
    print(f"  {'Mean (control)':<13} {'N/A':>20}", end="")

if control_suid:
    control_suid_mean = sum(control_suid) / len(control_suid)
    print(f"   {control_suid_mean:>19.2f}%")
else:
    print(f"   {'N/A':>20}")

### T-test ###

print("\n" + "=" * 70)
print(" " * 30 + "T-TEST")
print("=" * 70)

# show per-gene differences. positive differences means whales are more diverged for that specific gene
print("\n  Gene difference (suid identity - cetacean identity):")
print("  CANCER GENES:")
for i in range(len(cancer_diff_genes)):
    print(f"    {cancer_diff_genes[i]:<10} {cancer_diffs[i]:>+.2f}%")

print("\n  CONTROL GENES:")
for i in range(len(control_diff_genes)):
    print(f"    {control_diff_genes[i]:<10} {control_diffs[i]:>+.2f}%")

# needs at least 2 values in each group to run t-test
if len(cancer_diffs) >= 2 and len(control_diffs) >= 2:
    cancer_diff_mean = sum(cancer_diffs) / len(cancer_diffs)
    control_diff_mean = sum(control_diffs) / len(control_diffs)

    # Welch's t-test
    t_stat, p_value = stats.ttest_ind(cancer_diffs, control_diffs, equal_var=False)

    # print out mean differences for genes, t-statistic, and p-value
    print(f"\n  Mean difference (cancer genes):  {cancer_diff_mean:>+.2f}%  (n={len(cancer_diffs)})")
    print(f"  Mean difference (control genes): {control_diff_mean:>+.2f}%  (n={len(control_diffs)})")
    print(f"\n  t-statistic: {t_stat:.4f}")
    print(f"  p-value:     {p_value:.4f}")

    # quickly interpret the result using 0.05 signifcant level
    if p_value <= 0.05:
        print("\n  Result: SIGNIFICANT (p <= 0.05)")
        print("  Whale cancer genes are significantly more diverged than expected compared to suid baseline")
    else:
        print("\n  Result: NOT SIGNIFICANT (p > 0.05)")
        print("  No significant evidence that whale cancer gene divergence differs from what is expected based on suid baseline")
# if not enough values are provided for t-test, print out an error message
else:
    print("\n  ERROR: Not enough paired data to run t-test.")
    print(f"  Cancer gene pairs: {len(cancer_diffs)}, Control gene pairs: {len(control_diffs)}")
    t_stat = None
    p_value = None
print("=" * 70)

### Visual 2x2 Table Figure Creation ###

fig, ax = plt.subplots(figsize=(10, 5.5))
ax.axis("off")

# build the table row by row
col_labels = ["Gene", "Cetaceans\n(Blue vs Bowhead)", "Suids\n(Boar vs Warthog)"]
cell_data = []
cell_colors = []

# add a header row for cancer genes
cell_data.append(["CANCER GENES", "", ""])
cell_colors.append(["#87CEEB", "#87CEEB", "#87CEEB"])

# add each cancer gene as a row
for row in all_rows:
    if row["Gene_Type"] == "cancer":
        cet = row["Cetacean_Identity"]
        sui = row["Suid_Identity"]
        if cet != "N/A":
            cet = cet + "%"
        if sui != "N/A":
            sui = sui + "%"
        cell_data.append([row["Gene"], cet, sui])
        cell_colors.append(["lightcyan", "lightcyan", "lightcyan"])

# add cancer mean row
if cancer_cetacean:
    cet_mean_str = f"{sum(cancer_cetacean) / len(cancer_cetacean):.2f}%"
else:
    cet_mean_str = "N/A"
if cancer_suid:
    sui_mean_str = f"{sum(cancer_suid) / len(cancer_suid):.2f}%"
else:
    sui_mean_str = "N/A"
cell_data.append(["Mean", cet_mean_str, sui_mean_str])
cell_colors.append(["#B0E0E6", "#B0E0E6", "#B0E0E6"])

# add a header row for control genes
cell_data.append(["CONTROL GENES", "", ""])
cell_colors.append(["#FFD700", "#FFD700", "#FFD700"])

# add each control gene as a row
for row in all_rows:
    if row["Gene_Type"] == "control":
        cet = row["Cetacean_Identity"]
        sui = row["Suid_Identity"]
        if cet != "N/A":
            cet = cet + "%"
        if sui != "N/A":
            sui = sui + "%"
        cell_data.append([row["Gene"], cet, sui])
        cell_colors.append(["lightyellow", "lightyellow", "lightyellow"])

# add control mean row
if control_cetacean:
    cet_mean_str = f"{sum(control_cetacean) / len(control_cetacean):.2f}%"
else:
    cet_mean_str = "N/A"
if control_suid:
    sui_mean_str = f"{sum(control_suid) / len(control_suid):.2f}%"
else:
    sui_mean_str = "N/A"
cell_data.append(["Mean", cet_mean_str, sui_mean_str])
cell_colors.append(["#FFEEBB", "#FFEEBB", "#FFEEBB"])

# create the matplotlib table
table = ax.table(
    cellText=cell_data,
    colLabels=col_labels,
    cellColours=cell_colors,
    loc="center",
    cellLoc="center"
)

table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.2, 1.8)

# style the column header row (row 0 in the table object)
for col in range(3):
    table[0, col].set_facecolor("#4472C4")
    table[0, col].set_text_props(color="white", fontweight="bold")
    table[0, col].set_height(0.11)

# build a title that includes the p-value if available
if TITLE_ALLOWED == True:
    if p_value is not None:
        title_text = (
            "Table 1: Sequence Identity (%) Comparison of\n"
            "Cancer and Control Genes between Cetaceans and Suids\n"
            f"(p-value: {p_value:.4f})"
        )
    else:
        title_text = (
            "Table 1: Sequence Identity (%) Comparison of\n"
            "Cancer and Control Genes between Cetaceans and Suids\n"
        )
    ax.set_title(title_text, fontsize=14, fontweight="bold", y=1.1)

# adjust title and figure 

plt.tight_layout()
plt.savefig(table_figure, dpi=150, bbox_inches="tight")
# print where figure is saved to
print(f"\nTable figure saved to \nPATH: {table_figure}\n")
print("#" * 70)
print()
