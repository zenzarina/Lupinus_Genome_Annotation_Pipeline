## script per plottare distribuzione degli introni ed identificare quanti scartare 
## made by Chiara Zanoli 2026-06-09
## Functional Genomics Lab Univr 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

# --- CONFIGURATION ---
input_file = "junctions.bed"  # Your input file name
output_lengths = "intron_lengths.txt"
output_table = "intron_stats_summary.csv"
output_plot = "intron_distribution.png"

def calculate_intron_length(row):
    """Calculates the exact gap between the two exons in a HISAT2 BED12 row."""
    block_sizes = [int(x) for x in str(row[10]).strip(',').split(',')]
    block_starts = [int(x) for x in str(row[11]).strip(',').split(',')]
    
    # The intron starts at the end of the first exon and ends at the start of the second
    # Intron Length = (Start of Block 2) - (End of Block 1)
    exon1_end = block_starts[0] + block_sizes[0]
    exon2_start = block_starts[1]
    return exon2_start - exon1_end

# 1. Load data
print("Reading BED file...")
df = pd.read_csv(input_file, sep='\t', header=None)

# 2. Calculate lengths
# We apply the calculation to every row
df['intron_len'] = df.apply(calculate_intron_length, axis=1)

# Save Task 1: Lengths file
df['intron_len'].to_csv(output_lengths, index=False, header=False)
print(f"File 1 created: {output_lengths}")

# 3. Task 2: Generate Summary Table
stats_dict = {
    "Metric": ["Min", "Median", "Mean", "95th Percentile", "99th Percentile", "Total Count"],
    "Value_bp": [
        df['intron_len'].min(),
        df['intron_len'].median(),
        df['intron_len'].mean(),
        df['intron_len'].quantile(0.95),
        df['intron_len'].quantile(0.99),
        len(df)
    ]
}

stats_df = pd.DataFrame(stats_dict)
# Calculate cumulative counts for the specific breakpoints
stats_df['Cumulative_Count'] = stats_df['Value_bp'].apply(lambda x: (df['intron_len'] <= x).sum())
stats_df.to_csv(output_table, index=False)
print(f"File 2 created: {output_table}")

# 4. Task 3: Professional Visualization
print("Generating conference-ready plot...")
sns.set_theme(style="whitegrid")
fig, ax1 = plt.subplots(figsize=(10, 6))

# Histogram (Frequency) - Log Scale on X
sns.histplot(df['intron_len'], bins=100, log_scale=True, color="skyblue", alpha=0.7, ax=ax1, label="Frequency")
ax1.set_ylabel("Frequency (N Count)", fontsize=12, fontweight='bold')
ax1.set_xlabel("Intron Length (bp) [Log10 Scale]", fontsize=12, fontweight='bold')
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Twin axis for Cumulative Line
ax2 = ax1.twinx()
sorted_lens = np.sort(df['intron_len'])
cumulative = np.arange(len(sorted_lens)) + 1
ax2.plot(sorted_lens, cumulative, color="firebrick", lw=2, label="Cumulative Count")
ax2.set_ylabel("Cumulative N Count", color="firebrick", fontsize=12, fontweight='bold')
ax2.set_yscale('linear') # Keep cumulative count linear for clarity

# Annotations (Mean and 99th)
mean_val = df['intron_len'].mean()
p99_val = df['intron_len'].quantile(0.99)
mean_n = (df['intron_len'] <= mean_val).sum()
p99_n = (df['intron_len'] <= p99_val).sum()

for val, n, label, col in [(mean_val, mean_n, 'Mean', 'black'), (p99_val, p99_n, '99th%', 'darkgreen')]:
    ax1.axvline(val, color=col, linestyle='--', alpha=0.8)
    ax1.text(val, ax1.get_ylim()[1]*0.9, f' {label}\n {int(val)}bp', color=col, fontweight='bold', ha='left')

plt.title("Intron Length Distribution & Cumulative Frequency", fontsize=14, pad=20)
fig.tight_layout()
plt.savefig(output_plot, dpi=300)
print(f"File 3 created: {output_plot}")
print("Analysis complete.")
