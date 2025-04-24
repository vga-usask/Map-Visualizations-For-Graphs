import networkx as nx
import pygraphviz as pgv
from networkx.drawing.nx_agraph import write_dot
import os
import argparse
from pathlib import Path
import shutil
import math
from pygraphviz import AGraph
import subprocess

# Helper Functions
def calculate_centroid(positions):
    """Calculate the centroid from a list of positions."""
    x_coords = [pos[0] for pos in positions]
    y_coords = [pos[1] for pos in positions]
    centroid_x = sum(x_coords) / len(x_coords)
    centroid_y = sum(y_coords) / len(y_coords)
    return centroid_x, centroid_y

def calculate_distance(pos, centroid):
    """Calculate the Euclidean distance between a position and the centroid."""
    return math.sqrt((pos[0] - centroid[0])**2 + (pos[1] - centroid[1])**2)

def check_square_overlap(square1, square2):
    """
    Check if two squares overlap.
    Each square is represented as a dictionary with keys:
    - 'x': x-coordinate of the center
    - 'y': y-coordinate of the center
    - 'size': size of the square (width and height are equal)
    """
    half_size1 = square1['size'] / 2
    half_size2 = square2['size'] / 2

    # Calculate the boundaries of the squares
    left1, right1 = square1['x'] - half_size1, square1['x'] + half_size1
    top1, bottom1 = square1['y'] - half_size1, square1['y'] + half_size1

    left2, right2 = square2['x'] - half_size2, square2['x'] + half_size2
    top2, bottom2 = square2['y'] - half_size2, square2['y'] + half_size2

    # Check for overlap
    return not (right1 <= left2 or right2 <= left1 or bottom1 <= top2 or bottom2 <= top1)

def find_overlapping_squares(positions, size):
    """
    Find all overlapping square nodes.
    :param positions: Dictionary with node ID as key and (x, y) as value
    :param size: Size of each square (width and height are equal)
    :return: List of tuples representing overlapping square pairs
    """
    squares = [{'id': node, 'x': pos[0], 'y': pos[1], 'size': size[node]} for node, pos in positions.items()]
    overlaps = []

    for i in range(len(squares)):
        for j in range(i + 1, len(squares)):
            if check_square_overlap(squares[i], squares[j]):
                overlaps.append((squares[i]['id'], squares[j]['id']))
    return overlaps

def move_toward_centroid(node_pos, centroid, step_size):
    """
    Move a node toward the centroid by a step size.
    """
    dx, dy = centroid[0] - node_pos[0], centroid[1] - node_pos[1]
    distance = math.sqrt(dx**2 + dy**2)
    if distance == 0:  # Node is already at the centroid
        return node_pos
    unit_dx, unit_dy = dx / distance, dy / distance
    new_x = node_pos[0] + unit_dx * step_size
    new_y = node_pos[1] + unit_dy * step_size
    return (new_x, new_y)

def adjust_positions(positions, centroid, size, step_size):
    """
    Adjust the positions of squares to avoid overlaps and move closer to the centroid.
    """
    # Sort squares by their distance to the centroid
    sorted_nodes = sorted(positions.items(), key=lambda item: calculate_distance(item[1], centroid))
    updated_positions = positions.copy()
    counter=0
    for node, pos in sorted_nodes:
        last_valid_pos = updated_positions[node]  # Initialize with the current position of the node

        while True:
            # Move the square one step closer to the centroid
            new_pos = move_toward_centroid(last_valid_pos, centroid, step_size)
            updated_positions[node] = new_pos

            # Check if the square has reached the centroid
            if calculate_distance(new_pos, centroid) <= step_size:
                print(f"Node {node} reached the centroid at {new_pos}")
                break

            # Check for overlaps
            overlaps = find_overlapping_squares(updated_positions, size)
            if any(node in overlap for overlap in overlaps):
                # Revert to the last valid position if overlap occurs
                counter+=1
                print(f"Overlap detected for node {node}. Reverting to {last_valid_pos} polygons processed: {counter} ")

                updated_positions[node] = last_valid_pos
                break
            else:
                # Update the last valid position if no overlap occurs
                last_valid_pos = new_pos
                #print(f"Node {node} moved to {new_pos} without overlap.")

    return updated_positions

def save_to_dot_file(graph, updated_positions, output_path):
    """
    Save the graph with updated node positions to a DOT file.
    """
    for node in graph.nodes():
        if node in updated_positions:
            x, y = updated_positions[node]
            node.attr['pos'] = f"{x},{y}"
    graph.write(output_path)

def generate_svg(input_dot, output_svg):
    """Generate an SVG file using Graphviz."""
    os.system(f'neato -Gforcelabels=false -n2 -Tsvg  {input_dot} > {output_svg}')

# Main Function
def process_dot_file(dot_file_path, output_dot_file, step_size=1.1):
    """Process the DOT file and adjust node positions to avoid overlaps."""
    graph = AGraph(dot_file_path)
    positions = {}
    size = {}
    # Extract node positions from the DOT file
    for node in graph.nodes():
        pos = node.attr['pos']
        square_size = node.attr['height']
        if pos:
            x, y = map(float, pos.split(','))
            positions[node] = (x, y)
        if square_size:
            size[node]=float(square_size)*72
    # Calculate centroid
    centroid = calculate_centroid(list(positions.values()))
    print(f"Centroid: {centroid}")

    # Adjust positions to avoid overlaps and move toward centroid
    updated_positions = adjust_positions(positions, centroid, size, step_size)

    # Save the updated graph to a new DOT file
    save_to_dot_file(graph, updated_positions, output_dot_file)
    print(f"Updated DOT file saved to {output_dot_file}")
    return output_dot_file,positions,centroid


# Execution
# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('input_dir', type=str, help='path of input directory containing .dot files')
parser.add_argument('output_path', nargs='?', type=str, help='path for output directory')
args = parser.parse_args()
print("Line 14")
# Set up input and output directories
INPUT_PATH = Path(args.input_dir.rstrip())
if not INPUT_PATH.is_dir():
    raise ValueError("Input path must be a directory.")
print(INPUT_PATH)
try: 
    OUTPUT_PATH = Path(args.output_path.rstrip())
    OUTPUT_PATH.mkdir(parents=True, exist_ok=True)
except:
    OUTPUT_PATH = Path('tmp/')
    OUTPUT_PATH.mkdir(parents=True, exist_ok=True)
print(OUTPUT_PATH)
# Process each .dot file in the input directory
for dot_file in INPUT_PATH.glob("G.dot"):
    input_filename = dot_file.stem  # Get filename without extension
    #shutil.copy(Path(args.input_dir+"/"+input_filename+".dot"), Path(args.input_dir+'/'+dot_file.stem+"/G.dot"))
    #print("line 31")
    # Create the folder
    #os.makedirs(Path(args.input_dir+'/'+dot_file.stem), exist_ok=True)
    #print(Path(args.input_dir+'/'+dot_file.stem))
    # Load graph

    G = nx.Graph(pgv.AGraph(dot_file))

    # Create output paths with input filename prepended
    H_PATH = Path(args.input_dir+'/'+"H.dot")
    H_EMBEDDED_PATH = Path(args.input_dir+'/'+dot_file.stem+"_H_embedded.dot")
    H_EMBEDDED_PATH_svg =Path(args.input_dir+'/'+dot_file.stem+"_H_embedded.svg")

    # Parameterize this later?
    CLUSTER_BY = 'department'

    # Start filling up the nodes of H, corresponding to clusters (departments) of G.
    def generate_H_nodes():
        H = nx.Graph()
        cluster_index = 1
        for node in G.nodes():
            this_cluster_id = G.nodes[node][CLUSTER_BY]
            if this_cluster_id in H.nodes():
                H.nodes[this_cluster_id]['weight'] += 1
            else:
                H.add_node(this_cluster_id, weight=1, cluster=this_cluster_id, shape='rectangle')
                H.nodes[this_cluster_id][CLUSTER_BY] = this_cluster_id
                cluster_index += 1
        return H

    # Multiply node weight by scalars to get rectangle dimensions
    def add_rectangle_dimensions(H):
        width_scalar = 1/3
        height_scalar = 1/3
        for node in H.nodes():
            H.nodes[node]['width'] = H.nodes[node]['weight'] * width_scalar
            H.nodes[node]['height'] = H.nodes[node]['weight'] * height_scalar
            H.nodes[node]['label'] = node
            H.nodes[node]['department'] = node 
        return H

    # Generate weighted edges between departments
    def add_weighted_edges(H):
        for edge in G.edges():
            department_A = G.nodes[edge[0]][CLUSTER_BY]
            department_B = G.nodes[edge[1]][CLUSTER_BY]
            if department_A == department_B: continue
            if department_A in H[department_B]:
                H[department_A][department_B]["weight"] += 1
            else:
                H.add_edge(department_A, department_B)
                H[department_A][department_B]["weight"] = 1
        return H

    # Filter edges by weight


    # Generate layout network for each file
    print(f"Generating layout network for {dot_file} ... ", end="")
    H = generate_H_nodes()
    H = add_rectangle_dimensions(H)
    H = add_weighted_edges(H)

    # Write the file
    write_dot(H, H_PATH)
    os.system(f'neato -Gstart=123 -Goverlap=prism -Tdot -Gscale=0.5 -Gsplines="" {H_PATH} > {H_EMBEDDED_PATH}')
    os.system(f'neato  -Gstart=123 -Goverlap=prism -Tsvg -Gscale=0.5 -Gsplines="" {H_PATH} > {H_EMBEDDED_PATH_svg}')
    print("done")
    dot_file_path = H_EMBEDDED_PATH # Replace with your input DOT file path
    output_dot_file = H_EMBEDDED_PATH  # Output DOT file path
    output_svg = H_EMBEDDED_PATH_svg  # Output SVG file path
    updated_dot_file,positions,centroid = process_dot_file(dot_file_path, output_dot_file)
    generate_svg(updated_dot_file, output_svg)




#plot_positions(positions,7200,centroid)



