import xml.etree.ElementTree as ET
from shapely.geometry import Polygon, Point
import pygraphviz as pgv
import re


def parse_svg(svg_file):
    """
    Parse the SVG file to extract polygons and their points.
    """
    tree = ET.parse(svg_file)
    root = tree.getroot()

    polygons = []

    # Extract polygons
    for elem in root.iter():
        if elem.tag.endswith('polyline'):
            points = elem.attrib.get('points', '').strip()
            if points:
                points_list = [
                    tuple(map(float, point.split(','))) for point in points.split()
                ]
                if points_list[0] != points_list[-1]:
                    points_list.append(points_list[0])  # Ensure closed polygon
                polygons.append({"polygon": Polygon(points_list), "points": points_list})

    print(f"Extracted {len(polygons)} polygons from SVG.")
    return polygons


def parse_dot(dot_file):
    """
    Parse the DOT file to extract node positions and attributes.
    """
    graph = pgv.AGraph(dot_file)
    nodes = {}

    for node in graph.nodes():
        attributes = node.attr
        pos = attributes.get("pos")
        department = attributes.get("department")

        if pos and department:
            try:
                x, y = map(float, pos.split(','))
                nodes[node] = {"pos": (x, y), "department": int(department)}
            except ValueError:
                print(f"Invalid position or department for node {node}: pos={pos}, department={department}")

    print(f"Extracted {len(nodes)} nodes from DOT file.")
    return nodes


def flip_y(pos):
    """
    Flip the Y-axis for a position (DOT to SVG conversion).
    """
    x, y = pos
    return x, -y


def match_nodes_to_polygons(svg_file, dot_file, output_file):
    """
    Match nodes from the DOT file to polygons in the SVG file and write matches as PolID ppoints.
    """
    polygons = parse_svg(svg_file)
    nodes = parse_dot(dot_file)

    output_lines = []

    for polygon in polygons:
        polygon_points = polygon["points"]
        matched_department = None

        # Match the first DOT node that falls inside the polygon
        for node, node_data in nodes.items():
            flipped_pos = flip_y(node_data["pos"])
            node_point = Point(flipped_pos)
            if polygon["polygon"].contains(node_point):
                matched_department = node_data["department"]
                break  # Stop after finding the first matching node

        if matched_department is not None:
            # Ensure only one entry per polygon
            points_str = " ".join([f"{x},{y}" for x, y in polygon_points])
            output_lines.append(f"{matched_department} {points_str}")
        else:
            print("No matching department found for polygon.")

    # Write the results to the output file
    with open(output_file, "w") as f:
        for line in output_lines:
            f.write(line + "\n")

    print(f"Output written to {output_file}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python script.py <svg_file> <dot_file> <output_file>")
        sys.exit(1)

    svg_file = sys.argv[1].rstrip()
    dot_file = sys.argv[2].rstrip()
    output_file = sys.argv[3].rstrip()

    match_nodes_to_polygons(svg_file, dot_file, output_file)
