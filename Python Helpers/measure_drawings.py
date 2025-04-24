'''
    this measures a drawing of a cluster graph with polygon borders
    inputs are H, a graph of the clusters; G, a graph of the individual
    nodes; and drawing, an svg of the drawing of G and H
    inputs are actually folders of the above
    output is a table where rows are of each graph and columns
    are categories we are measuring and the aggregate score
'''
import argparse
import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from Drawing import Drawing
from collections import Counter

# argument stuff
parser = argparse.ArgumentParser()
parser.add_argument('main_directory', nargs=1, type=str, help='directory of drawing directories')
parser.add_argument('output_path', nargs='?', type=str, help='path for output csv')
args = parser.parse_args()
MAIN_DIR = args.main_directory[0]
try: 
    OUTPUT_PATH = Path(args.output_path)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    print(OUTPUT_PATH)
except:
    OUTPUT_PATH = 'drawing_measures_iqbal.csv'

# create Drawing objects 
def process_drawings():
    print("processing graph drawings from filesystem ...")
    drawing_objs = []
    drawing_dirs = [x[0] for x in os.walk(MAIN_DIR)][1:]
    drawing_index = 1
    for dir in drawing_dirs:
        this_name = os.path.basename(os.path.normpath(dir))
        print(f"  {drawing_index} out of {len(drawing_dirs)}: {this_name}")
        G_path = f'{dir}/G.dot'
        H_path = f'{dir}/H.dot'
        drawing_path = f'{dir}/drawing.svg'
        drawing_objs.append(Drawing(this_name, G_path, H_path, drawing_path))
        drawing_index += 1
    return drawing_objs

# get names of drawings for csv
def get_drawing_names(drawings):
    names = []
    for drawing in drawings:
        names.append(drawing.get_name())
    return names

# this makes sure all polygons correspond to a unique cluster id
def all_polygons_unique_cluster(drawings):
    constraint_met_for_drawings = []
    for drawing in tqdm(drawings, desc="seeing if all polygons have only one cluster id: "):
        cumulative_mapping_size = 0
        for polygon in drawing.get_polygons():
            number_of_cluster_ids = polygon.get_number_of_cluster_ids()
            if number_of_cluster_ids>1:
                print("----------------------------")
                print(number_of_cluster_ids)
            cumulative_mapping_size += number_of_cluster_ids
        if cumulative_mapping_size == len(drawing.get_polygons()):
            constraint_met_for_drawings.append(1)
        else:
            constraint_met_for_drawings.append(0)
    return constraint_met_for_drawings

# this makes sure that no cluster occurs in more than one polygon
def drawing_not_fragmented(drawings):
    constraint_met_for_drawings = []
    for drawing in tqdm(drawings, desc="seeing if all clusters have only one polygon (fragmentation): "):
        # we want a list of how many times a cluster occurs in the polygons of the drawing
        polygons = drawing.get_polygons()
        cluster_occurences = []
        for polygon in polygons:
            cluster_occurences.append(polygon.get_cluster_id())
        cluster_count = Counter(cluster_occurences)
        # clusters should occur once and only once; if greater than one, then the map is defragmented
        drawing_defragmented = False
        for cluster in cluster_count:
            if cluster_count[cluster] > 1:
                drawing_defragmented = True
                # print(f"{drawing.get_name()} {cluster} has {cluster_count[cluster]}")
        if drawing_defragmented == True:
            constraint_met_for_drawings.append(0)
        else:
            constraint_met_for_drawings.append(1)  
    return constraint_met_for_drawings

# the above is a constraint, but it's stil useful to see how fragmented a drawing is
# this returns the ratio of the nubmer of clusters in the graph to how many polygons there are
def get_fragmentation_amount(drawings):
    frag_perecentages = []
    for drawing in drawings:
        number_of_clusters = len(drawing.get_clusters())
        number_of_polygons = len(drawing.get_polygons())
        frag_perecentages.append(number_of_clusters / number_of_polygons)
    return frag_perecentages

# measure of the mean convexity of each polygon
def get_mean_convexity_score(drawings):
    mean_convexity_scores = []
    # get convexity score of all polygons
    for drawing in tqdm(drawings, desc="calculating mean convexity scores: "):
        sum_of_ratios = 0.0
        these_polygons = drawing.get_polygons()
        
        
        for polygon in these_polygons:
            print(polygon.get_perimeter())
            #print(polygon.get_convex_hull_perimeter())
            ch_perimeter_to_perimeter = polygon.get_convex_hull_perimeter() / polygon.get_perimeter()
            sum_of_ratios += ch_perimeter_to_perimeter
            #print(sum_of_ratios)
        if len(these_polygons)==0:
            print(f"{drawing.get_name()}**********")
        mean_convexity = sum_of_ratios / len(these_polygons)
        mean_convexity_scores.append(mean_convexity)
    return mean_convexity_scores

# get the space between polygons (M_2); ask about the measure itself tho
# may need to center everything to avoid fp precision errors
def get_space_between_polygons_score(drawings):
    mean_space_scores = []
    for drawing in tqdm(drawings, desc="calculating space between polygon scores: "):
        convex_hull_area = drawing.get_entire_convex_hull_area()
        polygon_area_sum = 0
        for polygon in drawing.get_polygons():
            polygon_area_sum += polygon.get_area()
        this_space_score = polygon_area_sum / convex_hull_area
        mean_space_scores.append(this_space_score)
    return mean_space_scores

# see how closesly the cluster sizes mathch the areas of the polygon (M_3)
def get_cluster_size_to_polygon_area(drawings):
    ratios = []
    for drawing in tqdm(drawings, desc='calculating ratios between cluster size and polygon area: '):
        cumulative_ratio_sum = 0.0
        cumulative_polygon_area = drawing.get_cumulative_polygon_areas()
        cumulative_cluster_size = drawing.get_cumulative_cluster_sizes()
        for polygon in drawing.get_polygons():
            area_ratio = polygon.get_area() / cumulative_polygon_area
            cluster_size_ratio = drawing.get_cluster_size(polygon.get_cluster_id()) / cumulative_cluster_size
            difference = abs(area_ratio - cluster_size_ratio)
            cumulative_ratio_sum += difference
        mean_ratio = cumulative_ratio_sum / len(drawing.get_polygons())
        ratios.append(mean_ratio)
    return ratios

# adjacency measure here
# this is different than what's in the paper
def get_immediate_neighbors_match_score(drawings, drawing_is_not_fragmented):
    match_scores = []
    index = 0
    for drawing in drawings:
        if drawing_is_not_fragmented[index] != 0:
            this_neighbor_mapping = drawing.get_neighbor_mapping()
            this_H_neighbor_mapping = drawing.get_H_neighbor_mapping()
            j_sum = 0
            # ask about iterating over polygons or nodes
            for polygon in this_neighbor_mapping:
                these_neighboring_polygons = this_neighbor_mapping[polygon]
                these_neighboring_H_nodes = this_H_neighbor_mapping[polygon.get_cluster_id()]
                for n_polygon in these_neighboring_polygons:
                    this_polygon_cluster_id = n_polygon.get_cluster_id()
                    if this_polygon_cluster_id in these_neighboring_H_nodes:
                        j_sum += 1
            this_score = j_sum / len(drawing.get_clusters())
            match_scores.append(this_score)
        else:
            match_scores.append(-1)
        index += 1
    return match_scores

# subgraph stress measure here

# get the spread for the interior subgraph of eacah polygon,
# then the average of that for the whole drawing
# this is the perimeter method
def get_mean_spread_score(drawings):
    mean_spreads_for_drawings = []
    for drawing in tqdm(drawings, desc="calculating mean spread scores: "):
        cumulative_spread = 0.0
        mean_denominator = len(drawing.get_polygons())
        for polygon in drawing.get_polygons():
            this_subgraph_ch_perimeter = polygon.get_interior_subgraph_perimeter()
            if this_subgraph_ch_perimeter == -1:
                mean_denominator -= 1
            else:
                this_polygon_perimeter = polygon.get_perimeter()
                this_spread = this_subgraph_ch_perimeter / this_polygon_perimeter
                cumulative_spread += this_spread
        mean_spread = cumulative_spread / mean_denominator
        mean_spreads_for_drawings.append(mean_spread)
    return mean_spreads_for_drawings

import networkx as nx
import numpy as np
# Helper function to calculate Euclidean distance
def euclidean_distance(pos1, pos2):
    return np.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)

# Function to calculate the shortest path distances between nodes in the graph
def get_shortest_path_distances(graph, nodes):
    """
    Calculate the shortest path distances between a set of nodes using the graph's structure.
    :param graph: NetworkX graph (self._G)
    :param nodes: List of node IDs for which to calculate distances
    :return: Dictionary of shortest path distances between node pairs
    """
    distances = {}
    for i, node1 in enumerate(nodes):
        for j, node2 in enumerate(nodes):
            if i < j:
                try:
                    # Compute shortest path distance between node1 and node2
                    distance = nx.shortest_path_length(graph, source=node1, target=node2)
                    distances[(node1, node2)] = distance
                except nx.NetworkXNoPath:
                    # If no path exists, set the distance to infinity (or a very large value)
                    distances[(node1, node2)] = float('inf')
    return distances

# Function to calculate stress for a drawing
def calculate_stress(drawing):
    stress_values = []

    # For each polygon, get the cluster ID and retrieve the corresponding nodes (EmbeddedPoint objects)
    for polygon in drawing.get_polygons():
        cluster_id = polygon.get_cluster_id()  # Get the cluster ID associated with the polygon
        
        # Retrieve nodes associated with the cluster by filtering the points in the drawing
        nodes = [point for point in drawing._points if point.get_cluster_id() == cluster_id]
        
        if not nodes:
            continue  # Skip if no nodes are found for this cluster

        # Get the node IDs to calculate shortest path distances
        node_ids = [point._node for point in nodes]

        # Calculate desired distances (shortest path distances) based on graph G
        desired_distances = get_shortest_path_distances(drawing._G, node_ids)
        actual_distances = {}
        
        # Calculate actual Euclidean distances between nodes
        for i, node1 in enumerate(nodes):
            for j, node2 in enumerate(nodes):
                if i < j:
                    pos1 = node1.get_coord_tuple()  # Position in the layout (SVG)
                    pos2 = node2.get_coord_tuple()
                    actual_distances[(node1._node, node2._node)] = euclidean_distance(pos1, pos2)
        
        # Calculate stress for this polygon
        polygon_stress = 0.0
        for (node1, node2), d_actual in actual_distances.items():
            d_desired = desired_distances.get((node1, node2), float('inf'))  # Use a default of infinity if no path
            if d_desired > 0 and d_desired != float('inf'):
                polygon_stress += ((d_actual - d_desired) ** 2) / d_desired
        
        stress_values.append(polygon_stress)
    
    # Return average stress for the drawing
    return np.mean(stress_values) if stress_values else 0.0

# Integrate stress calculation into the main data dictionary
def get_stress_measure(drawings):
    stress_scores = []
    for drawing in tqdm(drawings, desc="calculating stress for each drawing: "):
        stress_score = calculate_stress(drawing)
        stress_scores.append(stress_score)
    return stress_scores

import numpy as np
from itertools import combinations

# Helper function to check if two line segments (p1, p2) and (q1, q2) intersect
def do_edges_cross(p1, p2, q1, q2):
    """
    Returns True if line segments (p1, p2) and (q1, q2) cross each other.
    :param p1, p2: Endpoints of the first line segment (tuples).
    :param q1, q2: Endpoints of the second line segment (tuples).
    :return: Boolean indicating if the two edges cross.
    """
    def orientation(p, q, r):
        # Return the orientation of the triplet (p, q, r)
        val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
        if val == 0:
            return 0  # Collinear
        elif val > 0:
            return 1  # Clockwise
        else:
            return 2  # Counterclockwise

    def on_segment(p, q, r):
        if min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and min(p[1], r[1]) <= q[1] <= max(p[1], r[1]):
            return True
        return False

    # Calculate orientations needed for general and special cases
    o1 = orientation(p1, p2, q1)
    o2 = orientation(p1, p2, q2)
    o3 = orientation(q1, q2, p1)
    o4 = orientation(q1, q2, p2)

    # General case
    if o1 != o2 and o3 != o4:
        return True

    # Special cases (collinear points)
    if o1 == 0 and on_segment(p1, q1, p2):
        return True
    if o2 == 0 and on_segment(p1, q2, p2):
        return True
    if o3 == 0 and on_segment(q1, p1, q2):
        return True
    if o4 == 0 and on_segment(q1, p2, q2):
        return True

    return False

# Function to calculate edge crossings for a drawing
def calculate_edge_crossings(drawing):
    edge_crossings = 0
    edges = list(drawing._G.edges())  # Get the edges from the graph G
    
    # We must get the positions of the nodes for each edge
    node_positions = {point._node: point.get_coord_tuple() for point in drawing._points}
    
    # Check for crossings between all pairs of edges
    for edge1, edge2 in combinations(edges, 2):
        p1, p2 = node_positions[edge1[0]], node_positions[edge1[1]]
        q1, q2 = node_positions[edge2[0]], node_positions[edge2[1]]
        
        # Check if the two edges cross
        if do_edges_cross(p1, p2, q1, q2):
            edge_crossings += 1

    return edge_crossings

# Integrate edge crossing calculation into the main data dictionary
def get_edge_crossing_measure(drawings):
    crossing_scores = []
    for drawing in tqdm(drawings, desc="calculating edge crossings for each drawing: "):
        edge_crossings = calculate_edge_crossings(drawing)
        crossing_scores.append(edge_crossings)
    return crossing_scores
# Main function (this part remains unchanged from before)
def main():
    data = {}
    # load all drawing objects to conduct tests on 
    drawings = process_drawings()
    drawing_names = get_drawing_names(drawings)
    
    data['drawing'] = drawing_names
    # existing measures
    polygons_have_unique_clusters = all_polygons_unique_cluster(drawings)
    data['C_1'] = polygons_have_unique_clusters
    
    drawing_is_not_fragmented = drawing_not_fragmented(drawings)
    data['C_2'] = drawing_is_not_fragmented
    
    mean_convexity_scores = get_mean_convexity_score(drawings)
    data['M_1'] = mean_convexity_scores
    
    inter_polygon_space_scores = get_space_between_polygons_score(drawings)
    data['M_2'] = inter_polygon_space_scores
    
    immediate_neighbors_match_score = get_immediate_neighbors_match_score(drawings, drawing_is_not_fragmented)
    data['M_3'] = immediate_neighbors_match_score
    
    cluster_size_to_polygon_area = get_cluster_size_to_polygon_area(drawings)
    data['M_4'] = cluster_size_to_polygon_area
    
    mean_spread_scores = get_mean_spread_score(drawings)
    data['M_6'] = mean_spread_scores
    
    fragmentation_percentages = get_fragmentation_amount(drawings)
    data['frag'] = fragmentation_percentages
    
    # Add the stress measure
    stress_scores = get_stress_measure(drawings)
    data['M_stress'] = stress_scores

    # Add the edge crossing measure
    edge_crossing_scores = get_edge_crossing_measure(drawings)
    data['M_edge_crossings'] = edge_crossing_scores
    
    # write out as csv
    print("writing out")
    df = pd.DataFrame(data)
    df.sort_values('drawing')
    df.to_csv(OUTPUT_PATH)

if __name__ == '__main__':
    main()
