"""
analyze_results.py ---------------------------------------------------------------

Author: Paul Antezana

This program serves to read BLAST results, classify each gene as either
shared (CONSERVED) or divergent (SPECIES-SPECIFIC) by utilizing a set threshold.

It then uses these calculations to create visual charts representing the results.

BEFORE RUNNING:#########
  run blast_compare.py
########################
  
"""

# import necessary modules
import csv
import os
import matplotlib

# this is so the chart can save into a file later
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# build file paths relative to this script's location so it works from any directory
script_dir = os.path.dirname(os.path.abspath(__file__))
results_file = os.path.join(script_dir, "../results/blast_results.csv")
bar_chart = os.path.join(script_dir, "../results/figures/gene_identity_chart.png")
pie_chart = os.path.join(script_dir, "../results/figures/data_availability_chart.png")

# create the figures folder if it doesn't exist yet
os.makedirs(os.path.join(script_dir, "../results/figures"), exist_ok=True)

# genes with identity equal to or above this threshold are considered "shared" (conserved)
SHARED_THRESHOLD = 95.0

# list to hold all the genes of interest
GENES_TO_COMPARE = ["CIRBP", "RPA2", "PCNA", "ERCC1", "TP53"]

# these lists will be filled as we read the .csv
genes_with_data = []
identities = []
genes_no_data = []
all_genes_ordered = []

# read the csv and sort genes into two groups: those WITH data and those WITHOUT
with open(results_file, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene = row["Gene"]
        if gene not in GENES_TO_COMPARE:
            continue
        identity_str = row["Percent_Identity"]
        all_genes_ordered.append(gene)

        if identity_str == "N/A":
            genes_no_data.append(gene)
        else:
            genes_with_data.append(gene)
            identities.append(float(identity_str))

# this will be filled with genes either classified as shared or diverged based on threshold
shared_genes = []
diverged_genes = []

# formatting a table header for easier reading purposes
print("=" * 60)
print(f"{'Gene':<12} {'% Identity':>12}  {'Classification'}")
print("=" * 60)

# this runs for every gene in our ordered gene list one at a time
for gene in all_genes_ordered:
    if gene in genes_no_data:
        print(f"{gene:<12} {'—':>12}   No bowhead sequence in NCBI")
    else:
        idx = genes_with_data.index(gene)
        identity = identities[idx]
        # if the given % identity is EQUAL to or ABOVE the threshold, we consider it SHARED (i.e. conserved)
        if identity >= SHARED_THRESHOLD:
            label = "SHARED (conserved)"
            shared_genes.append(gene)
        # if the given % identity is BELOW the threshold, we consider it DIVERGED (i.e. species-specific)
        else:
            label = "DIVERGED (species-specific)"
            diverged_genes.append(gene)
        # print out the genes and their % identities
        print(f"{gene:<12} {identity:>11.1f}%  {label}")

# print a summary of the counts
print("=" * 60)
print(f"\nReal data:     {len(genes_with_data)} genes")
print(f"Shared:        {len(shared_genes)} genes: {shared_genes}")
print(f"Diverged:      {len(diverged_genes)} genes: {diverged_genes}")
print(f"No NCBI data:  {len(genes_no_data)} genes")

### GRAPH CREATION ###

# assign a color to each bar: blue = shared, red = diverged
colors = []
for i in identities:
    if i >= SHARED_THRESHOLD:
        colors.append("skyblue")
    else:
        colors.append("salmon")

# create and set up the bar chart
fig, ax = plt.subplots(figsize=(7, 5))
bars = ax.bar(genes_with_data, identities, color=colors, width=0.5)

# add a percent label on top of each bar
for bar, val in zip(bars, identities):
    ax.text(bar.get_x() + bar.get_width() / 2, val - 0.3,
            f"{val}%", ha="center", va="bottom", fontsize=11, fontweight="bold")

# draw a dashed line at the shared/diverged threshold. make it semi-transparent too so it can be easily read
ax.axhline(y=SHARED_THRESHOLD, color="green", linestyle="--", linewidth=1.2, alpha=0.75,
           label=f"Shared threshold ({SHARED_THRESHOLD}%)")

# label axes and title
ax.set_xlabel("Gene", fontsize=12, fontweight="bold")
ax.set_ylabel("Percent Identity (%)", fontsize=12, fontweight="bold")
ax.set_title("Bowhead vs. Blue Whale: Sequence Identity\n"
             "(Pairwise alignment of mRNA sequences)", fontsize=12, fontweight="bold")
# set a limit on the y-axis
ax.set_ylim(0, 110)

# enable legends and adjust position (so it's not in front of the bars)
ax.legend(bbox_to_anchor=(1.00, -0.15), loc='lower right', borderaxespad=0)

# save the figure
plt.tight_layout()
plt.savefig(bar_chart)
print(f"\nBar chart saved to {bar_chart}")
