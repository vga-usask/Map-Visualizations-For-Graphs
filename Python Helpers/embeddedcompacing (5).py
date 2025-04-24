import math
from pygraphviz import AGraph
import subprocess
import os
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
    squares = [{'id': node, 'x': pos[0], 'y': pos[1], 'size': size} for node, pos in positions.items()]
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
                print(f"Overlap detected for node {node}. Reverting to {last_valid_pos} ")
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
def process_dot_file(dot_file_path, output_dot_file, square_size=7200, step_size=1.1):
    """Process the DOT file and adjust node positions to avoid overlaps."""
    graph = AGraph(dot_file_path)
    positions = {}

    # Extract node positions from the DOT file
    for node in graph.nodes():
        pos = node.attr['pos']
        if pos:
            x, y = map(float, pos.split(','))
            positions[node] = (x, y)

    # Calculate centroid
    centroid = calculate_centroid(list(positions.values()))
    print(f"Centroid: {centroid}")

    # Adjust positions to avoid overlaps and move toward centroid
    updated_positions = adjust_positions(positions, centroid, square_size, step_size)

    # Save the updated graph to a new DOT file
    save_to_dot_file(graph, updated_positions, output_dot_file)
    print(f"Updated DOT file saved to {output_dot_file}")
    return output_dot_file,positions,centroid

import matplotlib.pyplot as plt

def plot_positions(positions, square_size, centroid=None):
    """
    Plot the positions of squares as boxes on a 2D plane.
    :param positions: Dictionary with node ID as key and (x, y) as value
    :param square_size: Size of each square (width and height are equal)
    :param centroid: Optional tuple (x, y) representing the centroid to plot
    """
    fig, ax = plt.subplots()
    half_size = square_size / 2

    # Plot squares
    for node, (x, y) in positions.items():
        rect = plt.Rectangle((x - half_size, y - half_size), square_size, square_size,
                             edgecolor='blue', facecolor='lightblue', alpha=0.5)
        ax.add_patch(rect)
        ax.text(x, y, str(node), color='black', ha='center', va='center', fontsize=8)

    # Plot centroid if provided
    if centroid:
        ax.plot(centroid[0], centroid[1], 'ro', label='Centroid')
        ax.legend()

    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.title('Node Positions with Squares')
    plt.savefig("embedded.png")
# Execution

dot_file_path = "G_H_embedded.dot"  # Replace with your input DOT file path
output_dot_file = "G_H_embedded_updated.dot"  # Output DOT file path
output_svg = "G_H_embedded_updated.svg"  # Output SVG file path

updated_dot_file,positions,centroid = process_dot_file(dot_file_path, output_dot_file)
generate_svg(updated_dot_file, output_svg)
plot_positions(positions,7200,centroid)


