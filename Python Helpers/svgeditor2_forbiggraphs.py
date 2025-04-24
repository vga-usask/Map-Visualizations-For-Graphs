import re
import argparse
from pathlib import Path

# Argument Parsing
parser = argparse.ArgumentParser()
parser.add_argument('drawing', type=str, help='SVG drawing file from GPU (e.g., base3 or others).')
parser.add_argument('base_drawing', type=str, help='Base SVG file with attributes (e.g., base2).')
parser.add_argument('replacement_drawing', type=str, help='Path to save intermediate SVG output.')
parser.add_argument('replacement_drawing_complete', type=str, help='Path to save the final SVG output.')
parser.add_argument('source_path_G', type=str, help='Source path for G.dot.')
parser.add_argument('destination_path_G', type=str, help='Destination path for G.dot.')
parser.add_argument('source_path_H', type=str, help='Source path for H.dot.')
parser.add_argument('destination_path_H', type=str, help='Destination path for H.dot.')
args = parser.parse_args()

# Resolve paths
source_path_G = Path(args.source_path_G.strip())
destination_path_G = Path(args.destination_path_G.strip())
source_path_H = Path(args.source_path_H.strip())
destination_path_H = Path(args.destination_path_H.strip())
drawing_path = Path(args.drawing.strip())
base_drawing_path = Path(args.base_drawing.strip())
replacement_path = Path(args.replacement_drawing.strip())
final_output_path = Path(args.replacement_drawing_complete.strip())

# Debug: Verify input paths
print(f"Drawing Path: {drawing_path}")
print(f"Base Drawing Path: {base_drawing_path}")
print(f"Intermediate Output Path: {replacement_path}")
print(f"Final Output Path: {final_output_path}")

# Read and Modify Polygons
try:
    with open(drawing_path, 'r') as drawing_file, open(base_drawing_path, 'r') as base_file:
        drawing_content = drawing_file.read()
        base_content = base_file.read()

    # Extract cluster and clustercolor from base_drawing
    polygons = re.findall(r'\[cluster\s*=\s*(\d+),\s*clustercolor\s*=\s*"([^"]+)",', base_content)
    
    # Debug: Log extracted polygons
#    print("Extracted polygons from base_drawing:")
 #   for title, fill in polygons:
  #      print(f"Cluster: {title}, Fill: {fill}")

    # Replace matching polygons in drawing
    for title, fill in polygons:
        # Regex pattern to find matching <title> and <polygon>
        pattern = rf'(<title>{title}</title>\s*<polygon[^>]*fill=")([^"]+)(" stroke=")([^"]+)'
        replacement = rf'\1{fill}\3\4'  # Replace the fill value, retain stroke
        
        # Debug: Log before and after replacement
        #if re.search(pattern, drawing_content):
       #     print(f"Replacing title {title} with fill {fill}")
       # else:
        #    print(f"Title {title} not found in drawing.")

        # Apply replacement
        drawing_content = re.sub(pattern, replacement, drawing_content)

    # Write intermediate SVG
    with open(replacement_path, 'w') as replacement_file:
        replacement_file.write(drawing_content)

    print(f"Intermediate SVG saved to: {replacement_path}")

except FileNotFoundError as e:
    print(f"Error: {e}")
    exit(1)

# Additional Replacements for Final SVG
try:
    # Regex for specific replacements
    rx_re = re.compile(r'rx="([^"]+)"')
    ry_re = re.compile(r'ry="([^"]+)"')
    font_size_re = re.compile(r'font-size="([^"]+)"')

    # Modify the intermediate SVG file
    with open(replacement_path, 'r', encoding='utf-8') as replacement_file, open(final_output_path, 'w', encoding='utf-8') as final_file:
        for line in replacement_file:
            # Apply attribute replacements
            line = rx_re.sub('rx="20"', line)
            line = ry_re.sub('ry="35"', line)
            line = font_size_re.sub('font-size="14"', line)
            final_file.write(line)

    print(f"Final SVG saved to: {final_output_path}")

except FileNotFoundError as e:
    print(f"Error: {e}")
    exit(1)

print("Replacements completed successfully.")
