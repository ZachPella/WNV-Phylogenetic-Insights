import baltic as bt
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from datetime import datetime

# Load your tree using baltic's loadNewick function
ll = bt.loadNewick('/content/tree_2025.nwk')

# Load your metadata
metadata = pd.read_csv('/content/updated_metadata.tsv', sep='\t')
print("Metadata shape:", metadata.shape)
print("Metadata columns:", metadata.columns.tolist())

# Extract year from date column
def extract_year(date_str):
    """Extract year from date string, handling various formats"""
    if pd.isna(date_str) or date_str == '':
        return None
    
    date_str = str(date_str)
    
    # Handle YYYY-XX-XX format
    if 'XX' in date_str:
        year = date_str.split('-')[0]
        try:
            return int(year)
        except:
            return None
    
    # Handle other date formats
    for fmt in ['%Y-%m-%d', '%Y/%m/%d', '%Y', '%m/%d/%Y']:
        try:
            return datetime.strptime(date_str, fmt).year
        except:
            continue
    
    # Try to extract just the year if it's a 4-digit number
    if len(date_str) >= 4 and date_str[:4].isdigit():
        return int(date_str[:4])
    
    return None

# Apply year extraction
metadata['year'] = metadata['date'].apply(extract_year)

# Remove rows without valid years
metadata_with_years = metadata.dropna(subset=['year'])
print(f"\nRows with valid years: {len(metadata_with_years)} out of {len(metadata)}")

# Get year counts
year_counts = metadata_with_years['year'].value_counts().sort_index()
print(f"\nYear range: {year_counts.index.min()} to {year_counts.index.max()}")
print(f"Years with data: {len(year_counts)}")

# Create the visualization
fig = plt.figure(figsize=(15, 10))
gs = gridspec.GridSpec(2, 2, height_ratios=[2, 1], width_ratios=[3, 1])

# Main bar plot - all years
ax1 = fig.add_subplot(gs[0, :])

# Create color array - highlight 2023 in a different color
colors = ['orange' if year == 2023 else 'steelblue' for year in year_counts.index]
bars = ax1.bar(year_counts.index, year_counts.values, alpha=0.7, color=colors, edgecolor='navy', linewidth=0.5)

# Add vertical line at 2020 cutoff
ax1.axvline(x=2019, color='red', linestyle='--', linewidth=2, alpha=0.8, label='2019 Cutoff')

# Create custom legend
from matplotlib.patches import Patch
legend_elements = [
    plt.Line2D([0], [0], color='red', linestyle='--', linewidth=2, alpha=0.8, label='2019 Cutoff'),
    Patch(facecolor='steelblue', alpha=0.7, label='Published genomes'),
    Patch(facecolor='orange', alpha=0.7, label='Genomes generated as part of this study')
]

# Customize main plot
ax1.set_xlabel('Year', fontsize=12, fontweight='bold')
ax1.set_ylabel('Number of Genomes', fontsize=12, fontweight='bold')
ax1.set_title('West Nile Virus Genome Counts by Year\nDemonstrating 2019 Cutoff Rationale', 
              fontsize=14, fontweight='bold', pad=20)
ax1.grid(True, alpha=0.3)
ax1.legend(handles=legend_elements, fontsize=10, loc='upper left')

# Set x-axis ticks for every year
all_years = list(range(int(year_counts.index.min()), int(year_counts.index.max()) + 1))
ax1.set_xticks(all_years)
ax1.set_xticklabels(all_years, rotation=45, ha='right')

# Add value labels on ALL bars
for bar in bars:
    height = bar.get_height()
    year = bar.get_x() + bar.get_width()/2
    ax1.text(year, height + max(year_counts.values) * 0.005, 
            f'{int(height)}', ha='center', va='bottom', fontsize=7, fontweight='bold')

# Bottom left: Pre-2020 vs Post-2020 comparison
ax2 = fig.add_subplot(gs[1, 0])
pre_2019 = metadata_with_years[metadata_with_years['year'] <= 2019]
post_2019 = metadata_with_years[metadata_with_years['year'] > 2019]

comparison_data = [len(pre_2019), len(post_2019)]
comparison_labels = [f'≤2019\n({len(pre_2019)} genomes)', f'>2019\n({len(post_2019)} genomes)']
colors = ['lightblue', 'lightcoral']

bars2 = ax2.bar(comparison_labels, comparison_data, color=colors, alpha=0.7, edgecolor='black')
ax2.set_title('Genome Distribution: Pre vs Post 2019', fontweight='bold')
ax2.set_ylabel('Number of Genomes')

# Add percentage labels
total = sum(comparison_data)
for bar, count in zip(bars2, comparison_data):
    percentage = (count/total) * 100
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(comparison_data) * 0.02,
             f'{percentage:.1f}%', ha='center', va='bottom', fontweight='bold')

# Bottom right: Recent years detail (2015-2025)
ax3 = fig.add_subplot(gs[1, 1])
recent_years = year_counts[year_counts.index >= 2015]
bars3 = ax3.bar(recent_years.index, recent_years.values, 
                color=['lightblue' if year <= 2019 else 'lightcoral' for year in recent_years.index],
                alpha=0.7, edgecolor='black')

ax3.axvline(x=2019, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax3.set_title('Recent Years Detail (2015+)', fontweight='bold')
ax3.set_xlabel('Year')
ax3.set_ylabel('Count')
ax3.tick_params(axis='x', rotation=45)

# Add value labels on all bars in the detail plot
for bar in bars3:
    height = bar.get_height()
    if height > 0:  # Only label non-zero bars
        ax3.text(bar.get_x() + bar.get_width()/2, height + max(recent_years.values) * 0.02,
                f'{int(height)}', ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.show()

# Print summary statistics
print("\n" + "="*50)
print("SUMMARY STATISTICS")
print("="*50)
print(f"Total genomes with year data: {len(metadata_with_years)}")
print(f"Genomes ≤2019: {len(pre_2019)} ({len(pre_2019)/len(metadata_with_years)*100:.1f}%)")
print(f"Genomes >2019: {len(post_2019)} ({len(post_2019)/len(metadata_with_years)*100:.1f}%)")

print(f"\nPost-2019 breakdown:")
post_2019_counts = post_2019['year'].value_counts().sort_index()
for year, count in post_2019_counts.items():
    print(f"  {int(year)}: {count} genomes")

print(f"\nAverage genomes per year:")
print(f"  Pre-2019 (1999-2019): {len(pre_2019)/22:.1f} genomes/year")
if len(post_2019) > 0:
    post_2019_years = len(post_2019_counts)
    print(f"  Post-2019: {len(post_2019)/post_2019_years:.1f} genomes/year")

print(f"\nPeak years:")
top_5_years = year_counts.nlargest(5)
for year, count in top_5_years.items():
    print(f"  {int(year)}: {count} genomes")
