import re
import os
import argparse
from pathlib import Path
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('drawing', type=str, help='drawing from gpu from our2 or....')
parser.add_argument('drawing_mapdot', type=str, help='drawing  from base')
parser.add_argument('base_drawing', type=str, help='drawing  from base2')
parser.add_argument('node_shape', type=str, help='drawing  from base2')
parser.add_argument('replacement_drawing',  type=str, help='final svg')
parser.add_argument('source_path_G',  type=str, help='source_path_G G.doot')
parser.add_argument('destination_path_G',  type=str, help='destination_path_G G.dot')
parser.add_argument('source_path_H',  type=str, help='source_path_H H.doot')
parser.add_argument('destination_path_H',  type=str, help='destination_path_H H.dot')
args = parser.parse_args()

source_path_G= Path(args.source_path_G.rstrip())
destination_path_G= Path(args.destination_path_G.rstrip())
source_path_H= Path(args.source_path_H.rstrip())
destination_path_H= Path(args.destination_path_H.rstrip())
# Set up input and output directories

# Paths to the .dot files

# Define file paths
#file1_path = Path(args.drawing)
# File paths for the two SVG files
file_path_drawing = Path(args.drawing.rstrip())
dot_fil=Path(args.drawing_mapdot.rstrip())
file_path_bdrawing = Path(args.base_drawing.rstrip())

output_path = Path(args.node_shape.rstrip())

final_output_path =  Path(args.replacement_drawing.rstrip())


#os.system(f'neato -Gforcelabels=false -Ecolor=grey  -Nshape=ellipse -n2 -Tsvg  {dot_fil} > {file_path_bdrawing}')
# Dictionary to store stroke color from drawing.svg based on titles (both orientations)
drawing_edge_colors = {}

# Regular expressions to capture edge group, title, and path stroke
edge_group_re = re.compile(r'<g[^>]*class="edge"[^>]*>')
title_re = re.compile(r'<title>(.*?)</title>')
stroke_re = re.compile(r'stroke="([^"]+)"')

# Step 1: Read drawing.svg line by line and store edge colors by title (in both orientations)
with open(file_path_drawing, 'r', encoding='utf-8') as drawing_file:
    current_title = None
    inside_edge_group = False
    for line in drawing_file:
        if edge_group_re.search(line):
            inside_edge_group = True
            current_title = None  # Reset title at each edge group
        elif inside_edge_group and title_re.search(line):
            current_title = title_re.search(line).group(1)
            # Store in both orientations (e.g., "142--257" and "257--142")
            src_dst = current_title
            dst_src = "&#45;&#45;".join(current_title.split("&#45;&#45;")[::-1])
        elif inside_edge_group and 'path' in line and current_title:
            stroke_match = stroke_re.search(line)
            if stroke_match:
                stroke = stroke_match.group(1)
                drawing_edge_colors[src_dst] = stroke
                drawing_edge_colors[dst_src] = stroke
            inside_edge_group = False  # End of this edge group

# Step 2: Read Bdrawing.svg, update matching edges, and collect unmatched titles
unmatched_titles = []
with open(file_path_bdrawing, 'r', encoding='utf-8') as bdrawing_file, open(output_path, 'w', encoding='utf-8') as output_file:
    current_title = None
    inside_edge_group = False
    for line in bdrawing_file:
        # Detect start of an edge group
        if edge_group_re.search(line):
            inside_edge_group = True
            current_title = None

        # Detect the title line in edge group
        elif inside_edge_group and title_re.search(line):
            current_title = title_re.search(line).group(1)

        # If inside edge group and a path line is found, apply the stroke if available
        elif inside_edge_group and 'path' in line:
            if current_title in drawing_edge_colors:
                stroke = drawing_edge_colors[current_title]
                line = stroke_re.sub(f'stroke="{stroke}"', line)
            else:
                # Collect unmatched titles
                unmatched_titles.append(current_title)
            inside_edge_group = False  # End of this edge group

        # Detect end of an edge group
        if '</g>' in line:
            inside_edge_group = False

        # Write the modified or unmodified line to the output file
        output_file.write(line)

# Print unmatched titles
print("Unmatched titles in Bdrawing.svg:")
for title in unmatched_titles:
    print(title)

import re

# File path for the modified SVG file



# Regular expressions for replacements
rx_re = re.compile(r'rx="27"')
ry_re = re.compile(r'ry="18"')
font_size_re = re.compile(r'font-size="14\.00"')

# Step: Read the SVG file, apply replacements, and save to a new file
with open(output_path, 'r', encoding='utf-8') as svg_file, open(final_output_path, 'w', encoding='utf-8') as final_output_file:
    for line in svg_file:
        # Apply replacements
        line = rx_re.sub('rx="5"', line)
        line = ry_re.sub('ry="8"', line)
        line = font_size_re.sub('font-size="8"', line)
        
        # Write modified line to the new file
        final_output_file.write(line)

print("Replacements completed. The final file is saved as Bdrawing_final_adjustment.svg.")
shutil.copy(source_path_G, destination_path_G)
shutil.copy(source_path_H, destination_path_H)