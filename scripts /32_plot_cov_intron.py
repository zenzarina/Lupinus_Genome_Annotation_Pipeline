# script for plotting intron coverage distribution 
# made by Chiara Zanoli 2026-06-09 
# Functional Genomics lab

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load your full dataset
# In a standard BED from regtools, score is column 5 (index 4)
df = pd.read_csv('junctions.bed', sep='\t', header=None, usecols=[4], names=['coverage'])

# 1. Calculate frequency distribution (How many introns for each coverage level)
freq_dist = df['coverage'].value_counts().sort_index()

# 2. Set up a professional plotting environment
plt.style.use('seaborn-v0_8-muted') 
fig, ax = plt.subplots(figsize=(12, 7))

# 3. Create a Bar Plot
# On a log scale, we adjust the width of the bars so they remain visible
# Width is set to a fraction of the x-value to keep them looking consistent
widths = freq_dist.index * 0.01 
ax.bar(freq_dist.index, freq_dist.values, 
       width=widths, color="#138ce9", alpha=0.7, 
       edgecolor="#0b538e", linewidth=0.5, label='Intron Frequency')

# 4. Use Log Scales for both axes
ax.set_xscale('log')
ax.set_yscale('log')

# 5. Aesthetics
ax.set_title('Global Intron Coverage Distribution (Barplot)', fontsize=16, fontweight='bold', pad=20)
ax.set_xlabel('Read Support (Coverage) [Log Scale]', fontsize=13)
ax.set_ylabel('Number of Unique Introns [Log Scale]', fontsize=13)

# Add gridlines for the log scale
ax.grid(True, which="major", ls="-", alpha=0.3, color='grey')
ax.grid(True, which="minor", ls=":", alpha=0.1, color='grey')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 6. Add the vertical filtering guide at coverage = 10
ax.axvline(x=10, color='#d7191c', linestyle='--', linewidth=2, alpha=0.8, 
           label='Threshold (10 Reads)')

# Annotate the threshold for clarity
ax.text(11, ax.get_ylim()[1]*0.5, 'Filter Zone', color='#d7191c', fontweight='bold', fontsize=10)

ax.legend(frameon=True, facecolor='white', framealpha=1)
plt.tight_layout()
plt.savefig('intron_global_coverage_barplot.png') 
