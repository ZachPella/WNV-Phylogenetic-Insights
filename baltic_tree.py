import baltic as bt
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import matplotlib as mpl

# Load your tree using baltic's loadNewick function
ll = bt.loadNewick('/content/tree_2025.nwk')

# Load your metadata
metadata = pd.read_csv('/content/updated_metadata.tsv', sep='\t')
print("Metadata shape:", metadata.shape)
print("Metadata columns:", metadata.columns.tolist())
print("\nRegion value counts:")
print(metadata['Region'].value_counts())

# Filter for Nebraska regions and 2023 samples
ne_regions = ['NE_Central', 'NE_West', 'NE_East']
ne_metadata = metadata[metadata['Region'].isin(ne_regions)]

# Extract year from date column and filter for 2023
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    # For YYYY-XX-XX format, just take the first 4 characters
    import re
    year_match = re.search(r'^(\d{4})', str(date_str))
    if year_match:
        return int(year_match.group(1))
    return None

metadata['year'] = metadata['date'].apply(extract_year)
ne_2023 = metadata[(metadata['Region'].isin(ne_regions)) & (metadata['year'] == 2023)]

print(f"\nFiltered to {len(ne_2023)} Nebraska 2023 samples")
if len(ne_2023) > 0:
    print("NE 2023 Region breakdown:")
    print(ne_2023['Region'].value_counts())

# Create strain-to-metadata mapping
strain_to_region = dict(zip(metadata['strain'], metadata['Region']))
strain_to_year = dict(zip(metadata['strain'], metadata['year']))
strain_to_date = dict(zip(metadata['strain'], metadata['date']))

# Create a mapping for Nebraska 2023 samples
ne_2023_strains = set(ne_2023['strain'])

# Convert dates to decimal years for time axis
def date_to_decimal_year(date_str):
    if pd.isna(date_str):
        return None
    try:
        # Parse YYYY-XX-XX format
        parts = str(date_str).split('-')
        if len(parts) >= 3:
            year = int(parts[0])
            month = int(parts[1]) if parts[1] != 'XX' else 6  # Default to mid-year if XX
            day = int(parts[2]) if parts[2] != 'XX' else 15    # Default to mid-month if XX
            
            # Convert to decimal year
            from datetime import datetime
            date_obj = datetime(year, month, day)
            year_start = datetime(year, 1, 1)
            year_end = datetime(year + 1, 1, 1)
            year_fraction = (date_obj - year_start).total_seconds() / (year_end - year_start).total_seconds()
            return year + year_fraction
        else:
            return int(parts[0])  # Just return year if only year provided
    except:
        return extract_year(date_str)  # Fallback to year extraction

# Add decimal years to metadata
metadata['decimal_year'] = metadata['date'].apply(date_to_decimal_year)
strain_to_decimal_year = dict(zip(metadata['strain'], metadata['decimal_year']))

# First, set absoluteTime for all nodes based on tree structure
# We'll use a recursive approach to set times from tips to root
def set_node_times(node):
    if node.branchType == 'leaf':
        # For tips, use the date from metadata
        if hasattr(node, 'name') and node.name:
            decimal_year = strain_to_decimal_year.get(node.name, None)
            if decimal_year is not None:
                node.absoluteTime = decimal_year
            else:
                node.absoluteTime = 2020  # Default
    else:
        # For internal nodes, set time based on children
        if hasattr(node, 'children') and node.children:
            # First, ensure all children have times set
            for child in node.children:
                set_node_times(child)
            # Set node time as the average of children times minus branch length
            child_times = [child.absoluteTime for child in node.children if hasattr(child, 'absoluteTime')]
            if child_times:
                node.absoluteTime = min(child_times) - (node.length if hasattr(node, 'length') else 0.1)
            else:
                node.absoluteTime = 2010  # Default for internal nodes

# Apply time setting to all nodes
for node in ll.Objects:
    set_node_times(node)

# Add traits to tree nodes for coloring
for node in ll.Objects:
    if hasattr(node, 'name') and node.name:
        strain = node.name
        region = strain_to_region.get(strain, 'other')
        year = strain_to_year.get(strain, None)
        
        # Create traits dictionary
        if not hasattr(node, 'traits'):
            node.traits = {}
        
        # Assign special highlighting for Nebraska 2023
        if strain in ne_2023_strains:
            if region == 'NE_Central':
                node.traits['highlight_color'] = '#FCB614'  # Yellow/Gold
                node.traits['is_ne_2023'] = True
            elif region == 'NE_West':
                node.traits['highlight_color'] = '#005E63'  # Teal
                node.traits['is_ne_2023'] = True
            elif region == 'NE_East':
                node.traits['highlight_color'] = '#AD122A'  # Red
                node.traits['is_ne_2023'] = True
            else:
                node.traits['highlight_color'] = '#666666'  # Dark grey
                node.traits['is_ne_2023'] = False
        else:
            # Muted colors for non-2023 samples
            node.traits['highlight_color'] = '#CCCCCC'  # Light grey
            node.traits['is_ne_2023'] = False
        
        node.traits['region'] = region
        node.traits['year'] = year

# Set up the figure - MODIFIED: Reduced height to compress the tree
fig, ax = plt.subplots(figsize=(20, 10), facecolor='w')  # Changed from 15 to 10

# Get tree dimensions
L = len(list(filter(lambda k: k.branchType == 'leaf', ll.Objects)))

# Define attribute functions
x_attr = lambda k: k.absoluteTime  # Use absolute time (years)

# Set up tree layout
ll.drawTree()  # This sets up the y coordinates

# MODIFIED: Compress the y-coordinates to make the tree more compact
y_compression_factor = 0.6  # Adjust this value to compress more (smaller) or less (larger)
for node in ll.Objects:
    if hasattr(node, 'y'):
        node.y = node.y * y_compression_factor

# Update ySpan after compression
ll.ySpan = ll.ySpan * y_compression_factor

# Function to determine branch color based on descendants
def get_branch_color(node):
    if hasattr(node, 'branchType') and node.branchType == 'leaf':
        # For tips, use their own color
        if hasattr(node, 'traits') and node.traits.get('is_ne_2023', False):
            return node.traits['highlight_color']
        else:
            return '#CCCCCC'
    else:
        # For internal nodes, check if they lead to Nebraska 2023 samples
        has_ne_2023 = False
        if hasattr(node, 'leaves'):
            for leaf in node.leaves:
                if hasattr(leaf, 'traits') and leaf.traits.get('is_ne_2023', False):
                    has_ne_2023 = True
                    break
        return '#666666' if has_ne_2023 else '#CCCCCC'

# Draw branches manually with proper coloring
for node in ll.Objects:
    if hasattr(node, 'parent') and node.parent:
        x = node.absoluteTime
        y = node.y
        x_parent = node.parent.absoluteTime
        y_parent = node.parent.y
        
        # Get color for this branch
        branch_color = get_branch_color(node)
        
        # Draw horizontal branch (from parent's x to node's x, at node's y)
        ax.plot([x_parent, x], [y, y], color=branch_color, linewidth=2, zorder=10, alpha=0.8)
        
        # Draw vertical connector at parent's x position
        if hasattr(node.parent, 'children') and len(node.parent.children) > 1:
            # Get y-coordinates of all siblings
            sibling_ys = [child.y for child in node.parent.children]
            # Draw vertical line from min to max y at parent's x position
            ax.plot([x_parent, x_parent], [min(sibling_ys), max(sibling_ys)], 
                   color=get_branch_color(node.parent), linewidth=2, zorder=9, alpha=0.8)

# Plot all tip points
for node in ll.Objects:
    if hasattr(node, 'branchType') and node.branchType == 'leaf':
        if hasattr(node, 'traits') and node.traits.get('is_ne_2023', False):
            # Nebraska 2023 samples - larger and colored
            ax.scatter(node.absoluteTime, node.y, s=100, c=node.traits['highlight_color'], 
                      zorder=20003, edgecolors='black', linewidth=1, alpha=1.0)
        else:
            # Other samples - smaller and grey
            ax.scatter(node.absoluteTime, node.y, s=30, c='#BBBBBB', 
                      zorder=20001, alpha=0.6)

# Customize the plot
ax.set_ylim(-10, ll.ySpan + 10)
ax.set_xlim(1993, 2025)  # Extend x-axis to start from 1993
[ax.spines[loc].set_visible(False) for loc in ['left', 'right', 'top']]

ax.grid(axis='x', ls='-', color='grey', alpha=0.3)
ax.tick_params(axis='y', size=0)
ax.tick_params(axis='x', labelsize=16)
ax.set_yticklabels([])
ax.set_xlabel('Time (Years)', fontsize=18)
ax.set_title('Phylogenetic Tree - Nebraska 2023 Samples Highlighted', fontsize=20, fontweight='bold', pad=20)

# MODIFIED: Create custom legend and move to top left
from matplotlib.patches import Patch
legend_elements = [
    plt.scatter([], [], c='#FCB614', s=150, edgecolors='black', linewidth=1, label='NE_Central (2023)'),
    plt.scatter([], [], c='#005E63', s=150, edgecolors='black', linewidth=1, label='NE_West (2023)'),
    plt.scatter([], [], c='#AD122A', s=150, edgecolors='black', linewidth=1, label='NE_East (2023)'),
    plt.scatter([], [], c='#BBBBBB', s=75, edgecolors='black', linewidth=0.5, label='Other samples'),
]

# MODIFIED: Changed location from 'upper right' to 'upper left' and adjusted bbox_to_anchor
ax.legend(handles=legend_elements, loc='upper left', fontsize=16, 
          title='Sample Types', title_fontsize=18, frameon=True, fancybox=True, shadow=True,
          bbox_to_anchor=(0.02, 0.98))

plt.tight_layout()
plt.show()

# Print statistics
print(f"\nTree Statistics:")
print(f"Total tips in tree: {len(list(filter(lambda k: k.branchType == 'leaf', ll.Objects)))}")
print(f"Total internal nodes: {len(list(filter(lambda k: k.branchType == 'node', ll.Objects)))}")

# Check which Nebraska 2023 samples are in the tree
tree_tips = set([tip.name for tip in ll.Objects if hasattr(tip, 'branchType') and tip.branchType == 'leaf'])
ne_2023_in_tree = ne_2023_strains.intersection(tree_tips)
print(f"\nNebraska 2023 samples in tree: {len(ne_2023_in_tree)}")
if len(ne_2023_in_tree) > 0:
    print("Sample names:", list(ne_2023_in_tree)[:10], "..." if len(ne_2023_in_tree) > 10 else "")

# Count by region in tree
ne_2023_regions_in_tree = {}
for strain in ne_2023_in_tree:
    region = strain_to_region.get(strain, 'unknown')
    ne_2023_regions_in_tree[region] = ne_2023_regions_in_tree.get(region, 0) + 1

print(f"Nebraska 2023 samples by region in tree:")
for region, count in ne_2023_regions_in_tree.items():
    print(f"  {region}: {count}")

# Optional: Save the figure
# plt.savefig('phylogenetic_tree_nebraska_2023_highlighted.png', dpi=300, bbox_inches='tight')
