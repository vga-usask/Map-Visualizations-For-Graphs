import os
import cairosvg
import argparse

def convert_svg_to_pdf(directory_path):
    # Traverse the directory tree
    for root, dirs, files in os.walk(directory_path):
        # Check if 'drawing.svg' exists in the current directory
        if 'drawing.svg' in files:
            # Get the relative path from the main directory (subdirectory name)
            relative_path = os.path.relpath(root, directory_path)
            # Use only the last part of the relative path (subdirectory name) for the PDF name
            subdirectory_name = os.path.basename(root)
            # Define the output PDF filename directly in the evaluation2 directory
            pdf_filename = os.path.join(directory_path, f"{subdirectory_name}.pdf")
            # Ensure the directory exists for the PDF (should be in the main directory)
            os.makedirs(os.path.dirname(pdf_filename), exist_ok=True)
            # Path to the SVG file
            svg_file = os.path.join(root, 'drawing.svg')
            try:
                # Convert SVG to PDF
                cairosvg.svg2pdf(url=svg_file, write_to=pdf_filename)
                print(f"Converted: {svg_file} -> {pdf_filename}")
            except Exception as e:
                print(f"Error converting {svg_file}: {e}")

# Set up argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Convert SVG files to PDFs in subdirectories.")
    parser.add_argument('main_directory', type=str, help="Directory to search for drawing.svg files")
    return parser.parse_args()

if __name__ == "__main__":
    # Parse the command-line arguments
    args = parse_args()
    
    # Get the main directory path from the command-line argument
    main_directory = args.main_directory
    
    # Convert all drawing.svg files to PDFs in subdirectories
    convert_svg_to_pdf(main_directory)
