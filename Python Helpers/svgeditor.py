import re
import os
import argparse
from pathlib import Path

def extract_map_attributes(file_path):
    """Extracts cluster and clustercolor from map.dot."""
    clusters = {}
    with open(file_path, 'r') as file:
        node_data = {}
        for line in file:
            if 'cluster=' in line:
                node_data['cluster'] = int(line.split('=')[1].split(',')[0].strip())
            if 'clustercolor="' in line:
                node_data['clustercolor'] = line.split('=')[1].replace('"', '').strip()
                clusters[node_data['cluster']] = node_data['clustercolor']
    return clusters

def extract_drawing_attributes(file_path):
    """Extracts department and clustercolor from drawing.dot in the correct order."""
    departments = {}
    with open(file_path, 'r') as file:
        node_data = {}
        for line in file:
            if 'clustercolor="' in line:
                node_data['clustercolor'] = line.split('=')[1].replace('"', '').strip()
            if 'department=' in line:
                node_data['department'] = int(line.split('=')[1].split(',')[0].strip())
                if 'clustercolor' in node_data:
                    departments[node_data['department']] = node_data['clustercolor']
                    node_data = {}
    return departments

def replace_clustercolor(drawing_file, map_data, drawing_data):
    """Replaces clustercolor in drawing_file based on matching department and cluster."""
    with open(drawing_file, 'r') as file:
        drawing_content = file.read()
    
    replacements_made = 0
    # Replace colors based on department-cluster matching
    for department, drawing_color in drawing_data.items():
        if department in map_data:
            replacement_color = map_data[department]
            # Replace all occurrences of the drawing color with the replacement color
            drawing_content = drawing_content.replace(drawing_color[:-1], replacement_color[:-1])
            replacements_made += 1

    # Write the updated content back to drawing.dot
    with open(drawing_file, 'w') as file:
        file.write(drawing_content)
    
    print(f"Total replacements made: {replacements_made}")
# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('drawing_file', type=str, help='drawing.dot from base3 or....')
parser.add_argument('map_file', type=str, help='map.dot from base2')
parser.add_argument('final_svg',  type=str, help='final svg')
args = parser.parse_args()

# Set up input and output directories

# Paths to the .dot files
drawing_file = Path(args.drawing_file.rstrip())
map_file = Path(args.map_file.rstrip())
final_svg=Path(args.final_svg.rstrip())
final_svg2=os.path.join(os.path.dirname(args.final_svg.rstrip()), "drawing.svg")

# Extract clustercolor values from map.dot using cluster attribute
map_data = extract_map_attributes(map_file)
# Extract clustercolor values from drawing.dot using department attribute
drawing_data = extract_drawing_attributes(drawing_file)

# Perform the color replacement in drawing.dot
replace_clustercolor(drawing_file, map_data, drawing_data)

print(f"Cluster colors have been updated in drawing.dot based on map.dot. {drawing_file}-{final_svg}")
os.system(f'neato -Gforcelabels=false  -Ecolor=grey28  -Nshape=ellipse -n2 -Tsvg  {drawing_file} > {final_svg}')

# Additional Replacements for Final SVG
try:
    # Regex for specific replacements
    rx_re = re.compile(r'rx="([^"]+)"')
    ry_re = re.compile(r'ry="([^"]+)"')
    font_size_re = re.compile(r'font-size="([^"]+)"')

    # Modify the intermediate SVG file
    with open(final_svg, 'r', encoding='utf-8') as replacement_file, open(final_svg2, 'w', encoding='utf-8') as final_file:
        for line in replacement_file:
            # Apply attribute replacements
            line = rx_re.sub('rx="20"', line)
            line = ry_re.sub('ry="35"', line)
            line = font_size_re.sub('font-size="14"', line)
            final_file.write(line)

    print(f"Final SVG saved to: {final_svg}")

except FileNotFoundError as e:
    print(f"Error: {e}")
    exit(1)

print("Replacements completed successfully.")