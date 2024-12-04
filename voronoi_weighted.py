import random
import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.affinity import translate

# Function to read input file and parse squares
def read_input_file(filename):
    ids = []
    square_data = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            square_id = int(parts[0])
            points = []
            for coord in parts[1:]:
                x, y = map(float, coord.split(','))
                points.append((x, y))
            ids.append(square_id)
            square_data.append(points)
    return ids, square_data

# Function to divide a square into 8 points: 4 corners, 4 midpoints on each edge, and 4 centroids of smaller squares
def get_square_points(polygon_points):
    # Extract corner points
    top_left, top_right, bottom_right, bottom_left = polygon_points
    
    # Calculate midpoints of each edge
    midpoints = [
        ((top_left[0] + top_right[0]) / 2, (top_left[1] + top_right[1]) / 2),  # Top edge midpoint
        ((top_right[0] + bottom_right[0]) / 2, (top_right[1] + bottom_right[1]) / 2),  # Right edge midpoint
        ((bottom_right[0] + bottom_left[0]) / 2, (bottom_right[1] + bottom_left[1]) / 2),  # Bottom edge midpoint
        ((bottom_left[0] + top_left[0]) / 2, (bottom_left[1] + top_left[1]) / 2)   # Left edge midpoint
    ]

    # Calculate centroids of four smaller squares
    mid_x = (top_left[0] + top_right[0]) / 2
    mid_y = (top_left[1] + bottom_left[1]) / 2
    centroids = [
        ((top_left[0] + mid_x) / 2, (top_left[1] + mid_y) / 2),  # Top left smaller square centroid
        ((mid_x + top_right[0]) / 2, (top_right[1] + mid_y) / 2),  # Top right smaller square centroid
        ((bottom_left[0] + mid_x) / 2, (mid_y + bottom_left[1]) / 2),  # Bottom left smaller square centroid
        ((mid_x + bottom_right[0]) / 2, (mid_y + bottom_right[1]) / 2)   # Bottom right smaller square centroid
    ]

    # Return 4 corners, 4 midpoints, and 4 centroids
    return [top_left, top_right, bottom_right, bottom_left] + midpoints + centroids

# Function to combine Voronoi polygons for each original square
from shapely.errors import TopologicalError

def combine_voronoi_polygons_for_square(indices, vor, boundary_box):
    polygons = []
    for i in indices:
        region = vor.regions[vor.point_region[i]]
        if region and -1 not in region:  # Ensure region is finite
            try:
                polygon_points = [vor.vertices[j] for j in region]
                poly = Polygon(polygon_points)
                if poly.is_valid:
                    polygons.append(poly)
                else:
                    print(f"Warning: Polygon {i} is not valid. Skipping.")
            except IndexError:
                print(f"Error: Invalid index in Voronoi region for point {i}. Skipping this point.")
            except TopologicalError:
                print(f"Topological error for region {region}. Skipping.")

    # Attempt to union polygons and clip to the boundary
    if polygons:
        combined_polygon = unary_union(polygons)
        clipped_polygon = combined_polygon.intersection(boundary_box)
        return clipped_polygon
    else:
        print(f"Warning: No valid polygons found for indices {indices}. Returning empty geometry.")
        return Polygon()  # Return an empty Polygon if no valid regions are found


# Function to write output to a file
def write_output_file(filename, ids, combined_polygons):
    with open(filename, 'w') as file:
        for square_id, poly in zip(ids, combined_polygons):
            if isinstance(poly, Polygon):
                points = " ".join([f"{x},{y}" for x, y in poly.exterior.coords])
                file.write(f"{square_id} {points}\n")

# Function to read edge list from a file
def read_edge_list(filename):
    edges = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split('--')
            if len(parts) == 2:
                src = int(parts[0].strip())
                dst = int(parts[1].strip())
                edges.append((src, dst))
    return edges

# Function to adjust polygons to create space between connected polygons
def adjust_polygons(combined_polygons, edges):
    translations = {i: (0, 0) for i in range(len(combined_polygons))}
    distance_factor = 10  # Increased distance factor for better separation
    for src, dst in edges:
        if src - 1 < len(combined_polygons) and dst - 1 < len(combined_polygons):
            poly_src = combined_polygons[src - 1]
            poly_dst = combined_polygons[dst - 1]

            if poly_src.intersects(poly_dst):
                delta_x, delta_y = poly_src.centroid.x - poly_dst.centroid.x, poly_src.centroid.y - poly_dst.centroid.y
                magnitude = np.hypot(delta_x, delta_y)
                if magnitude != 0:
                    delta_x /= magnitude
                    delta_y /= magnitude

                translations[dst - 1] = (translations[dst - 1][0] + delta_x * distance_factor, translations[dst - 1][1] + delta_y * distance_factor)
    
    adjusted_polygons = []
    for i, poly in enumerate(combined_polygons):
        if isinstance(poly, Polygon):
            adjusted_poly = translate(poly, xoff=translations[i][0], yoff=translations[i][1])
            adjusted_polygons.append(adjusted_poly)
        else:
            adjusted_polygons.append(poly)
    
    return adjusted_polygons

# Main function to process the input file, generate the Voronoi polygons, and save the result
def process_voronoi(input_file, output_file, edge_file):
    # Step 1: Read the input file
    ids, square_data = read_input_file(input_file)

    # Step 2: Calculate bounding box from input data
    all_x_coords = [point[0] for square in square_data for point in square]
    all_y_coords = [point[1] for square in square_data for point in square]
    min_x, max_x = min(all_x_coords), max(all_x_coords)
    min_y, max_y = min(all_y_coords), max(all_y_coords)
    width = max_x - min_x
    height = max_y - min_y

    # Step 3: Get 12 points for each square (4 corners + 4 midpoints + 4 centroids)
    points = []
    for square in square_data:
        points.extend(get_square_points(square))

    # Step 4: Convert points to NumPy array for Voronoi and add boundary points
    points = np.array(points)
    boundary_buffer = 500  # Add some buffer around the boundary for better visualization
    boundary_points = [
        (min_x - boundary_buffer, min_y - boundary_buffer),
        (min_x - boundary_buffer, max_y + boundary_buffer),
        (max_x + boundary_buffer, min_y - boundary_buffer),
        (max_x + boundary_buffer, max_y + boundary_buffer)
    ]

    all_points = np.vstack([points, boundary_points])

    # Step 5: Compute Voronoi diagram with extended boundary points
    vor = Voronoi(all_points)

    # Step 6: Create the boundary box
    boundary_box = box(min_x, min_y, max_x, max_y)

    # Step 7: Combine the Voronoi polygons for each original square, clipped to the boundary
    combined_polygons = []
    for i in range(0, len(points), 12):
        combined_polygon = combine_voronoi_polygons_for_square([i + j for j in range(12)], vor, boundary_box)
        combined_polygons.append(combined_polygon)

    # Step 8: Read the edge list and adjust polygons to create space
    edges = read_edge_list(edge_file)
    adjusted_polygons = adjust_polygons(combined_polygons, edges)

    # Step 9: Write the output to the file
    write_output_file(output_file, ids, adjusted_polygons)

    print(f"Processed Voronoi polygons saved to {output_file}")
    return points, square_data, adjusted_polygons, vor, min_x, max_x, min_y, max_y

# Example usage
import networkx as nx
import pygraphviz as pgv
from networkx.drawing.nx_agraph import write_dot
import os
import argparse
from pathlib import Path

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('input_dir', type=str, help='path of input directory containing .dot files')
parser.add_argument('input_dir2', type=str, help='path of input directory containing .dot files')
parser.add_argument('output_path', nargs='?', type=str, help='path for output directory')
args = parser.parse_args()

# Set up input and output directories
INPUT_PATH = Path(args.input_dir)
INPUT_PATH2 = Path(args.input_dir2)
if not INPUT_PATH.is_dir():
    raise ValueError("Input path must be a directory.")

try: 
    OUTPUT_PATH = Path(args.output_path.rstrip())
    OUTPUT_PATH.mkdir(parents=True, exist_ok=True)
except:
    OUTPUT_PATH = Path('tmp/')
    OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

# Process each .dot file in the input directory
for dot_file in INPUT_PATH.glob("G.dot"):
    input_filename = args.input_dir  # Get filename without extension
    
    output_filename = args.output_path # Get filename without extension
    input_filename2 = args.input_dir2


    input_file = args.input_dir2+'/'+'polygon_points.txt'
    output_file = args.output_path+'/'+'Vpolygon_points.txt'
    edge_file = args.input_dir2+'/'+'edgelist.txt'
    points, square_data, combined_polygons, vor, min_x, max_x, min_y, max_y = process_voronoi(input_file, output_file, edge_file)

    # Plot the original boxes, Voronoi polygons, and combined polygons
    fig, ax = plt.subplots()

    # Plot original boxes
    for square in square_data:
        x_coords, y_coords = zip(*square + [square[0]])
        ax.plot(x_coords, y_coords, 'g--', linewidth=1)

    # Plot Voronoi polygons before combining
    for region_index in vor.point_region:
        region = vor.regions[region_index]
        if not -1 in region and len(region) > 0:
            polygon = [vor.vertices[i] for i in region]
            poly = Polygon(polygon)
            x, y = poly.exterior.xy
            ax.plot(x, y, 'r-', linewidth=0.5)

    # Plot combined polygons
    for poly in combined_polygons:
        if isinstance(poly, Polygon):
            x, y = poly.exterior.xy
            ax.plot(x, y, 'k-', linewidth=1)
            color = (random.random(), random.random(), random.random())
            ax.fill(x, y, color=color, alpha=0.6)

    # Plot seed points (corners, midpoints, and centroids)
    for point in points:
        ax.plot(point[0], point[1], 'bo', markersize=5)
    for poly in combined_polygons:
        if isinstance(poly, Polygon):
            x, y = poly.exterior.xy
            ax.plot(x, y, 'r-', linewidth=1)
            #ax.fill(x, y, 'lightblue', alpha=0.4)
    # Set plot limits to match the bounding box
    plt.xlim(min_x - 50, max_x + 50)
    plt.ylim(min_y - 50, max_y + 50)

    # Add labels and title
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Combined Voronoi Polygons for Original Squares with Separation")

    # Save the figure as a PNG file
    plt.savefig(args.output_path+'/'+"combined_voronoi_for_squares.png")
