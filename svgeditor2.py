import re
import os
import argparse
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument('drawing', type=str, help='drawing from gpu from base3 or....')
parser.add_argument('base_drawing', type=str, help='drawing  from base2')
parser.add_argument('replacement_drawing',  type=str, help='final_temp svg')
parser.add_argument('replacement_drawing_complete',  type=str, help='final svg')
parser.add_argument('source_path_G',  type=str, help='source_path_G G.dot')
parser.add_argument('destination_path_G',  type=str, help='destination_path_G G.dot')
parser.add_argument('source_path_H',  type=str, help='source_path_H H.doot')
parser.add_argument('destination_path_H',  type=str, help='destination_path_H H.dot')
args = parser.parse_args()
source_path_G= Path(args.source_path_G.rstrip())
destination_path_G= Path(args.destination_path_G.rstrip())
source_path_H= Path(args.source_path_H.rstrip())
destination_path_H= Path(args.destination_path_H.rstrip())
final_output_path=Path(args.replacement_drawing_complete.rstrip())
# Set up input and output directories

# Paths to the .dot files

# Define file paths
file1_path = Path(args.drawing.rstrip())
file2_path = Path(args.base_drawing.rstrip())
output_path =Path(args.replacement_drawing.rstrip())

# Read file contents
with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
    file1_content = file1.read()
    file2_content = file2.read()

# Find all title and polygon pairs in file2 with their fill and stroke values
file2_polygons = re.findall(r'<title>(\d+)</title>\s*<polygon fill="([^"]+)" stroke="([^"]+)"', file2_content)

# Replace matching polygons in file1
for title, fill, stroke in file2_polygons:
    # Construct a regex pattern to find the corresponding polygon in file1 with the specified title
    pattern = rf'(<title>{title}</title>\s*<polygon fill=")([^"]+)(" stroke=")([^"]+)'
    # Replacement string with new fill and stroke values from file2
    replacement = rf'\1{fill}\3{stroke}'
    # Substitute in file1 content
    file1_content = re.sub(pattern, replacement, file1_content)

# Write the modified content to a new file
with open(output_path, 'w') as output_file:
    output_file.write(file1_content)



import argparse
from pathlib import Path
import shutil
# Regular expressions for replacements
rx_re = re.compile(r'rx="([^"]+)"')
ry_re = re.compile(r'ry="([^"]+)"')
font_size_re = re.compile(r'font-size="([^"]+)"')

# Step: Read the SVG file, apply replacements, and save to a new file
with open(output_path, 'r', encoding='utf-8') as svg_file, open(final_output_path, 'w', encoding='utf-8') as final_output_file:
    for line in svg_file:
        # Apply replacements
        line = rx_re.sub('rx="20"', line)
        line = ry_re.sub('ry="35"', line)
        line = font_size_re.sub('font-size="14"', line)
        
        # Write modified line to the new file
        final_output_file.write(line)

print("Replacements completed. The final file is saved as Bdrawing_final_adjustment.svg.")

method_graph_name = os.path.basename(os.path.dirname(args.replacement_drawing_complete.rstrip()))

import cairosvg

cairosvg.svg2pdf(url=args.replacement_drawing_complete.rstrip(), write_to=os.path.dirname(args.replacement_drawing_complete.rstrip())+".pdf")