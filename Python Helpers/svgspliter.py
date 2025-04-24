import os

# Input SVG file
input_svg = "drawing.svg"

# Define the output directory
output_dir = "splited"
os.makedirs(output_dir, exist_ok=True)

# Original viewBox values (from your SVG)
x_start = -31576.7
y_start = -1.80436e+06
total_width = 1.98453e+06
total_height = 1.83235e+06

# Define grid dimensions
rows = 2  # Number of horizontal rows
cols = 5  # Number of vertical columns

# Calculate piece dimensions
piece_width = total_width / cols
piece_height = total_height / rows

# Read the original SVG content
with open(input_svg, "r") as f:
    svg_content = f.read()

# Generate SVG pieces
for row in range(rows):
    for col in range(cols):
        x_offset = x_start + col * piece_width
        y_offset = y_start + row * piece_height
        viewBox = f"{x_offset} {y_offset} {piece_width} {piece_height}"
        
        # Replace the viewBox in the SVG content
        piece_svg = svg_content.replace(
            f'viewBox="{x_start} {y_start} {total_width} {total_height}"',
            f'viewBox="{viewBox}"'
        )
        
        # Save each piece as a separate SVG file
        output_file = os.path.join(output_dir, f"piece_{row}_{col}.svg")
        with open(output_file, "w") as output:
            output.write(piece_svg)
        print(f"Generated: {output_file}")
