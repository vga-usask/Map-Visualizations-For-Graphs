import re
from collections import defaultdict
from math import sqrt
import os
import argparse
from pathlib import Path
import shutil
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

parser = argparse.ArgumentParser()
parser.add_argument('input_filename', type=str, help='drawing from gpu from base3 or....')
parser.add_argument('output_filename', type=str, help='drawing  from base2')
parser.add_argument('output_filename_node_positions', type=str, help='drawing  from base2')
parser.add_argument('real_number_of_nodes', type=int, help='drawing  from base2')
parser.add_argument('min_distance',  type=int, help='final svg')

args = parser.parse_args()


input_filename=Path(args.input_filename.rstrip())
output_filename = Path(args.output_filename.rstrip())
output_filename2 = Path(args.output_filename.rstrip()+"_updated.dot")
output_filename2_svg = Path(args.output_filename.rstrip()+"_updated.svg")
output_filename_node_positions =Path(args.output_filename_node_positions.rstrip())
real_number_of_nodes=int(args.real_number_of_nodes)
min_distance=int(args.min_distance)

#input_filename="base3/bookland/map.dot"
#output_filename = "base3/bookland/CBA_H_embeddedd.dot"
#output_filename_node_positions = "base3/bookland/node_positions.txt"
#real_number_of_nodes=499
#min_distance=7700

#input_filename="base3/musicland/map.dot"
#output_filename = "base3/musicland/CBA_H_embeddedd.dot"
#output_filename_node_positions = "base3/musicland/node_positions.txt"
#real_number_of_nodes=250
#min_distance=6100

#input_filename="base3/recipies/map.dot"
#output_filename = "base3/recipies/CBA_H_embedded.dot"
#output_filename_node_positions = "base3/recipies/node_positions.txt"
#real_number_of_nodes=369
#min_distance=9600

#input_filename="base3/tradeland/map.dot"
#output_filename = "base3/tradeland/CBA_H_embedded.dot"
#output_filename_node_positions = "base3/tradeland/node_positions.txt"
#real_number_of_nodes=150
#min_distance=7000

#input_filename="base3/universities/map.dot"
#output_filename = "base3/universities/CBA_H_embedded.dot"
#output_filename_node_positions = "base3/universities/node_positions.txt"
#real_number_of_nodes=161
#min_distance=3000
new_G_output_filename = "g/G.dot" # replace with your desired output file path


import re

def parse_graph_file(filename):
    nodes = {}
    edges = []


    # Accumulate lines into full entries (nodes or edges) before parsing
    with open(filename, 'r') as file:
        entry = ""
        for line in file:
            entry += line.strip()  # Remove leading/trailing whitespace and join lines

            # Check if we’ve reached the end of a node or edge entry
            if line.strip().endswith("];"):
                # Try matching as a node first
                node_match1 = re.match(
                    r'(\d+)\s+\[cluster=(\d+),\s*clustercolor="([^"]+)",\s*height=([\d.]+),\s*name=([^"]+),\s*pos="([^"]+)",\s*weight=(\d+),\s*width=([\d.]+)\];',
                    entry
                )
                node_match2 = re.match(
                    r'(\d+)\s+\[cluster=(\d+),\s*clustercolor="([^"]+)",\s*height=([\d.]+),\s*name="([^"]+)",\s*pos="([^"]+)",\s*weight=(\d+),\s*width=([\d.]+)\];',
                    entry
                )
                node_match3 = re.match(
                    r'(\d+)\s+\[cluster=(\d+),\s*clustercolor="([^"]+)",\s*height=([\d.]+),\s*name=(\d+),\s*pos="([^"]+)",\s*weight=(\d+),\s*width=([\d.]+)\];',
                    entry
                )
                node_match3 = re.match(
                    r'(\d+)\s+\[cluster=(\d+),\s*clustercolor="([^"]+)",\s*height=([\d.]+),\s*pos="([^"]+)",\s*weight=(\d+),\s*width=([\d.]+)\];',
                    entry
                )
                node_match4 = re.match(
                    r'(\d+)\s+\[cluster=(\d+),\s*clustercolor="([^"]+)",\s*department=(\d+),\s*label=(\d+),\s*pos="([^"]+)",\s*weight=(\d+)\];',
                    entry
                )
                if node_match1 or node_match2 or node_match3 or node_match4:
                    if node_match1:
                        node_id, cluster, clustercolor, height, name, pos, weight, width = node_match1.groups()
                        print(node_id)
                        node_id = int(node_id)
                        if node_id <= real_number_of_nodes:
                            nodes[node_id] = {
                            "department": int(cluster),
                            "cluster": int(cluster),
                            "name": name,
                            "clustercolor": clustercolor,
                            "label": '',
                            "pos":tuple(map(lambda coord: float(coord) / 72, pos.split(','))),
                            "height": float(height),
                            "weight": int(weight),
                            "width": float(width)
                        }
                    if node_match2:
                        node_id, cluster, clustercolor, height, name, pos, weight, width = node_match2.groups()
                        
                        print(node_id)
                        node_id = int(node_id)
                        if node_id <= real_number_of_nodes:
                            nodes[node_id] = {
                            "department": int(cluster),
                            "cluster": int(cluster),
                            "name": name,
                            "clustercolor": clustercolor,
                            "label": '',
                            "pos":tuple(map(lambda coord: float(coord) / 72, pos.split(','))),
                            "height": float(height),
                            "weight": int(weight),
                            "width": float(width)
                        }
                    if node_match3:
                        node_id, cluster, clustercolor, height, pos, weight, width = node_match3.groups()
                        
                        print(node_id)
                        node_id = int(node_id)
                        if node_id <= real_number_of_nodes:
                            nodes[node_id] = {
                            "department": int(cluster),
                            "cluster": int(cluster),
                            "name": node_id,
                            "clustercolor": clustercolor,
                            "label": '',
                            "pos":tuple(map(lambda coord: float(coord) / 72, pos.split(','))),
                            "height": float(height),
                            "weight": int(weight),
                            "width": float(width)
                        }
                    if node_match4:
                        node_id, cluster, clustercolor, department,label, pos, weight = node_match4.groups()
                        
                        print(node_id)
                        node_id = int(node_id)
                        if node_id <= real_number_of_nodes:
                            nodes[node_id] = {
                            "department": int(cluster),
                            "cluster": int(cluster),
                            "name": node_id,
                            "clustercolor": clustercolor,
                            "label": '',
                            "pos":tuple(map(lambda coord: float(coord) / 72, pos.split(','))),
                            "weight": int(weight)
                           
                        }                            
                
                # Try matching as an edge
                else:
                    edge_match = re.match(
                        r'(\d+)\s*--\s*(\d+)\s+\[pos="([^"]+)",\s*type="([^"]+)",\s*weight=(\d+)\];',
                        entry
                    )
                    if edge_match:
                        src, dst, pos, edge_type, weight = edge_match.groups()
                        src = int(src)
                        dst = int(dst)
                        if src <= real_number_of_nodes and dst <= real_number_of_nodes:
                            edge_points = [tuple(map(float, point.split(','))) for point in pos.split()]
                            edges.append({
                                "src": src,
                                "dst": dst,
                                "type": edge_type,
                                "weight": int(weight),
                                "pos": edge_points
                            })

                # Clear entry for the next node or edge
                entry = ""
       
    return nodes, edges


def calculate_bounding_boxes_and_attributes(nodes):
    clusters = defaultdict(list)
    for node_id, data in nodes.items():
        clusters[data["cluster"]].append(data["pos"])

    cluster_attributes = {}
    for cluster_id, positions in clusters.items():
        # Calculate average x and y for the center
        avg_x = sum(pos[0] for pos in positions) / len(positions)
        avg_y = sum(pos[1] for pos in positions) / len(positions)
        
        # Calculate the box size based on the number of nodes
        box_size = len(positions)  # This can be scaled if needed
        width = box_size
        height = box_size
        
        # Define the bounding box's center and dimensions
        center_x = avg_x
        center_y = avg_y

        # Store the attributes
        cluster_attributes[cluster_id] = {
            "height": height,
            "width": width,
            "center": [center_x, center_y],
            "weight": len(positions)
        }
        
    return cluster_attributes

def check_overlap(attr1, attr2):
    # Calculate the center distance between two clusters
    dist_x = attr1["center"][0] - attr2["center"][0]
    dist_y = attr1["center"][1] - attr2["center"][1]
    distance = sqrt(dist_x ** 2 + dist_y ** 2)

    # Calculate the minimum allowed distance
    required_distance = (sqrt((attr1["width"]**2)+(attr1["height"]**2)) + sqrt((attr2["width"]**2)+(attr2["height"]**2)))  + min_distance
    return distance < required_distance, dist_x, dist_y, required_distance

def resolve_overlaps(cluster_attributes, step=10):
    # Continue adjusting positions until no overlaps exist
    overlaps = True
    while overlaps:
        overlaps = False
        for id1, attr1 in cluster_attributes.items():
            for id2, attr2 in cluster_attributes.items():
                if id1 >= id2:
                    continue
                has_overlap, dist_x, dist_y, required_distance = check_overlap(attr1, attr2)
                if has_overlap:
                    # Calculate displacement step
                    step_x = step * dist_x / required_distance
                    step_y = step * dist_y / required_distance
                    
                    # Move both clusters away from each other
                    attr1["center"][0] += step_x
                    attr1["center"][1] += step_y
                    attr2["center"][0] -= step_x
                    attr2["center"][1] -= step_y
                    overlaps = True
        

def calculate_inter_cluster_edges(edges, nodes):
    cluster_edges = defaultdict(lambda: defaultdict(int))

    for edge in edges:
        src_cluster = nodes[edge["src"]]["cluster"]
        dst_cluster = nodes[edge["dst"]]["cluster"]
        
        # Only consider edges between different clusters
        if src_cluster != dst_cluster:
            # Count the edge between the clusters
            cluster_edges[min(src_cluster, dst_cluster)][max(src_cluster, dst_cluster)] += 1

    return cluster_edges

def write_cluster_file(filename, cluster_attributes, cluster_edges):
    with open(filename, 'w') as file:
        file.write("strict graph \"\" {\n")
        file.write("    graph [bb=\"0,0,21605,17647\",\n")
        file.write("        overlap=prism\n")
        file.write("    ];\n")
        file.write("    node [label=\"\\N\"];")
        
        # Write each cluster as a node with bounding box attributes
        for cluster_id, attrs in cluster_attributes.items():
            pos_str = f"{attrs['center'][0]},{attrs['center'][1]}"
            file.write(f'    {cluster_id} [cluster={cluster_id},\n')
            file.write(f'        department={cluster_id},\n')
            file.write(f'        height={attrs["height"]},\n')
            file.write(f'        pos="{pos_str}",\n')
            file.write(f'        shape=rectangle,\n')
            file.write(f'        weight={attrs["weight"]},\n')
            file.write(f'        width={attrs["width"]}];\n')

        # Write inter-cluster edges with weights
        for src_cluster, dst_clusters in cluster_edges.items():
            for dst_cluster, weight in dst_clusters.items():
                file.write(f'    {src_cluster} -- {dst_cluster} [weight={weight}];\n')

        file.write("}\n")
def write_new_G_file(filename, nodes, edges):
    with open(filename, 'w') as file:
        file.write("strict graph \"\" {\n")
        file.write("node [label=\"\\N\"];\n")
        
        # Write each cluster as a node with bounding box attributes
        for node, attrs in nodes.items():
            file.write(f'{node} [department={attrs["department"]},\n')
            file.write(f'name=\"{attrs["department"]}\",\n')
            file.write(f'label=\"{node}\",\n')
            file.write(f'weight={0}];\n')
        

        # Write inter-cluster edges with weights
        for edge in edges:
            file.write(f'    \"{edge["src"]}\" -- \"{edge["dst"]}\" [type=\"t,l\",\n')
            file.write('weight=1];\n')

        file.write("}\n")
def write_node_positions( output_filename_node_positions,nodes, edges):
    with open(output_filename_node_positions, 'w') as file:
        # Write each cluster as a node with bounding box attributes
        for node, attrs in nodes.items():
            pos_str = f"{attrs['pos'][0]},{attrs['pos'][1]}"
            file.write(f'{node} {pos_str}\n')
# Example usage
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



nodes, edges = parse_graph_file(input_filename)
cluster_attributes = calculate_bounding_boxes_and_attributes(nodes)

# Adjust positions to remove overlaps
resolve_overlaps(cluster_attributes)

cluster_edges = calculate_inter_cluster_edges(edges, nodes)
write_cluster_file(output_filename, cluster_attributes, cluster_edges)
write_node_positions( output_filename_node_positions,nodes, edges)
#write_new_G_file(new_G_output_filename, nodes, edges,node_map)
print("Cluster attributes and inter-cluster edges have been written to", output_filename)
os.system(f'neato -Gforcelabels=false -Ecolor=grey28  -Nshape=ellipse -n2 -Tsvg  {output_filename} > {output_filename}".svg"')
updated_dot_file,positions,centroid = process_dot_file(output_filename, output_filename2)
generate_svg(updated_dot_file, output_filename2_svg)