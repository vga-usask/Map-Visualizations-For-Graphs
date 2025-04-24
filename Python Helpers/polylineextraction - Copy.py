import xml.etree.ElementTree as ET
from shapely.geometry import Polygon, Point
import sys

def extract_polylines_with_ids(svg_file, svg_file2, output_file):
    # Parse both SVG files
    tree1 = ET.parse(svg_file)
    root1 = tree1.getroot()
    tree2 = ET.parse(svg_file2)
    root2 = tree2.getroot()

    polylines = []
    texts = []

    # Extract polylines from the first SVG file
    for elem in root1.iter():
        if elem.tag.endswith('polyline'):
            points = elem.attrib.get('points', '').strip()
            if points:
                # Parse points and ensure all are captured
                points_list = []
                for point in points.split():
                    try:
                        x, y = point.split(',')
                        points_list.append((float(x), float(y)))
                    except ValueError:
                        print(f"Skipping invalid point format: {point}")
                if points_list and points_list[0] != points_list[-1]:
                    points_list.append(points_list[0])
                if points_list:
                    print(f"Extracted polyline points: {points_list}")  # Debug: Show extracted points
                    polylines.append(points_list)

    # Extract text elements and <title> elements as potential IDs from the second SVG file
    polygon_ids = {}
    current_id = None

    for elem in root2.iter():
        if elem.tag.endswith('title'):
            # Read <title> content as the polygon ID
            current_id = elem.text.strip() if elem.text else None

        if elem.tag.endswith('polygon') and current_id:
            # Associate the <polygon> points with the ID from the <title>
            points = elem.attrib.get('points', '').strip()
            if points:
                points_list = []
                for point in points.split():
                    try:
                        x, y = point.split(',')
                        points_list.append((float(x), float(y)))
                    except ValueError:
                        print(f"Skipping invalid point format in polygon: {point}")
                
                if points_list:
                    polygon_ids[current_id] = points_list
            current_id = None  # Reset after associating

        if elem.tag.endswith('text'):
            # Extract text position and content for additional checking
            try:
                x = float(elem.attrib.get('x', 0))
                y = float(elem.attrib.get('y', 0))
                text_content = elem.text.strip() if elem.text else ''
                
                # Add text as polygon ID if it matches a polygon title
                if text_content in polygon_ids:
                    texts.append((text_content, x, y))  # Store as (ID, x, y)
            except ValueError:
                print(f"Skipping invalid text coordinates: {elem.attrib.get('x', 0)}, {elem.attrib.get('y', 0)}")

    polyline_data = []

    # Match each polyline in the first file with a polygon ID from the second file
    for points_list in polylines:
        polygon = Polygon(points_list)
        polyline_id = None

        # Check if any text or <title> ID is within the polygon in the first file
        for text, tx, ty in texts:
            text_point = Point(tx, ty)
            if polygon.contains(text_point):
                polyline_id = text  # Assign the text or title ID as the polygon ID
                break  # Stop after finding the first match

        if polyline_id:
            print(f"Assigned polygon ID {polyline_id} for polygon with points: {points_list}")
            polyline_data.append(f"{int(polyline_id)} {' '.join([f'{x},{y}' for x, y in points_list])}")
        else:
            print(f"No matching text or title found for polygon with points: {points_list}")

    # Write the results to the output file
    with open(output_file, 'w') as f:
        for line in polyline_data:
            f.write(line + '\n')

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <svg_file1> <svg_file2> <output_file>")
        sys.exit(1)

    svg_file = sys.argv[1].rstrip()
    svg_file2 = sys.argv[2].rstrip()
    output_file = sys.argv[3].rstrip()
    
    extract_polylines_with_ids(svg_file, svg_file2, output_file)
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    main()
