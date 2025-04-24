'''
    this embeds nodes of G into the sections generate by embed_H.py
'''
import networkx as nx
import pygraphviz as pgv
from networkx.drawing.nx_agraph import write_dot
import os, argparse
from pathlib import Path
from tqdm import tqdm
import numpy as np

print("hi")
# parse the args
parser = argparse.ArgumentParser()
parser.add_argument('G_path', nargs=1, type=str, help='path of dot file of G')
parser.add_argument('H_embedded_path', nargs=1, type=str, help='path of dot file of H')
parser.add_argument('node_position', nargs=1, type=str, help='path of node_position.txt')
parser.add_argument('output_path', nargs='?', type=str, help='path for output directory')
args = parser.parse_args()
GPOS_MODIFIED_PATH = args.node_position[0]
G_PATH = args.G_path[0]
H_EMBEDDED_PATH = args.H_embedded_path[0]

try: 
    OUTPUT_PATH = Path(args.output_path.rstrip())
    OUTPUT_PATH.mkdir(parents=True, exist_ok=True)
except:
    OUTPUT_PATH = Path('output/')
    OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

G_OUT_PATH = OUTPUT_PATH.joinpath("G_pos.dot")
H_OUT_PATH = OUTPUT_PATH.joinpath("H_pos.dot")
DRAWING_DOT_PATH = OUTPUT_PATH.joinpath("drawing.dot")
DRAWING_PATH = OUTPUT_PATH.joinpath("drawing.svg")

# param this later
CLUSTER_BY = "department"

# map cluster ids -> nodes of that cluster_id (partitions nodes of G)
def generate_cluster_id_to_nodes_obj(G):
    cluster_id_to_nodes_object = {}
    for node in G.nodes():
        this_node = G.nodes[node]
        this_cluster_id = this_node[CLUSTER_BY]
        if this_cluster_id in cluster_id_to_nodes_object:
            cluster_id_to_nodes_object[this_cluster_id].append(node)
        else:
            cluster_id_to_nodes_object[this_cluster_id] = [node]
    return cluster_id_to_nodes_object

# this just maps cluster ids to consecutive integers {1, ..., k}; sfdp works best
# when clustering over these
def generate_cluster_id_to_int_obj(cluster_id_to_nodes_obj):
    cluster_id_to_int_obj = {}
    i = 1
    for cluster_id in cluster_id_to_nodes_obj:
        cluster_id_to_int_obj[cluster_id] = str(i)
        i += 1
    return cluster_id_to_int_obj

def hide_inter_cluster_edges(G):
    for edge in G.edges():
        cluster_A = G.nodes[edge[0]]['cluster']
        cluster_B = G.nodes[edge[1]]['cluster']
        if cluster_A != cluster_B:
            G.edges[edge]['style'] = 'invis'


# This function will extract positions from 'gpos_modified.dot' and apply them to the cluster's dot file
def integrate_positions(cluster_dot_file, node_positions_file):
    # Load the cluster dot file using PyGraphviz
    cluster_graph = pgv.AGraph(cluster_dot_file)

    # Read the node positions from node_positions.txt
    node_positions = {}
    with open(node_positions_file, 'r') as f:
        for line in f:
            # Parse each line: expected format is 'node_id posX,posY'
            parts = line.strip().split()
            if len(parts) == 2:
                node_id = parts[0]
                position = parts[1]  # This is 'posX,posY'
                node_positions[node_id] = position

    # Update the positions in the cluster graph using the parsed node positions
    for node_id, pos in node_positions.items():
        if cluster_graph.has_node(node_id):  # Check if the node exists in the cluster graph
            # Update the position in the cluster graph
            cluster_graph.get_node(node_id).attr['pos'] = pos

    # Write the updated cluster graph back to the cluster_dot_file (tmp_drawing.dot)
    cluster_graph.write(cluster_dot_file)

def main():
    print(f"generating drawing with G from {G_PATH} and H from {H_EMBEDDED_PATH} ... ", end="")
    H = nx.Graph(pgv.AGraph(H_EMBEDDED_PATH))
    G = nx.Graph(pgv.AGraph(G_PATH))
    cluster_id_to_nodes_obj = generate_cluster_id_to_nodes_obj(G)
    cluster_id_to_int_obj = generate_cluster_id_to_int_obj(cluster_id_to_nodes_obj)

    # where majority of generation takes place; also the majority of the math for the map embedding
    for cluster_id in tqdm(cluster_id_to_nodes_obj):
        # for each cluster_id, make a subgraph from G (and dot file)
        cluster_subgraph = G.subgraph(cluster_id_to_nodes_obj[cluster_id])
        os.system('mkdir -p tmp')
        write_dot(cluster_subgraph, 'tmp/tmp.dot')
        os.system('sfdp -Gnodesep=0.001 -Gstart=123 -Goverlap=prism -Tdot -Gsplines="" tmp/tmp.dot > tmp/tmp_drawing.dot')
        integrate_positions('tmp/tmp_drawing.dot', GPOS_MODIFIED_PATH)
        # load this dot -> pygraphviz graph -> networkx graph
        c_sg_drawing = nx.Graph(pgv.AGraph('tmp/tmp_drawing.dot'))
        
        # this is all math for sizing the cluster's representation ...
        x_coord=[ float(c_sg_drawing.nodes[n]['pos'].split(',')[0]) for n in c_sg_drawing.nodes()]
        width=max(x_coord)-min(x_coord)
        y_coord=[ float(c_sg_drawing.nodes[n]['pos'].split(',')[1]) for n in c_sg_drawing.nodes()]
        height=max(y_coord)-min(y_coord)
        x_scale_factor= float(H.nodes[cluster_id]['width'])*72  / width if width>0 else 1
        y_scale_factor= float(H.nodes[cluster_id]['height'])*72  / height if height>0 else 1
        x_center=(min(x_coord) + width/2) * x_scale_factor
        y_center=(min(y_coord) + height/2) * y_scale_factor
        x_translation=float(H.nodes[cluster_id]['pos'].split(',')[0]) - x_center
        y_translation=float(H.nodes[cluster_id]['pos'].split(',')[1]) - y_center
       
        # working in G (the general graph) here
        # for all nodes in this cluster
        #   reposition the ndoe according to the above math
        for node in cluster_subgraph.nodes():
            # this is only needed if there's no label field in the data -- FIX LATER
            try:
                G.nodes[node]['label'] = node#cluster_id#G.nodes[node]['name']
            except:
                pass

            x_s =  x_scale_factor * float(c_sg_drawing.nodes[node]['pos'].split(',')[0]) + x_translation
            y_s = y_scale_factor * float(c_sg_drawing.nodes[node]['pos'].split(',')[1]) + y_translation
            pos=str(x_s) + ","+str(y_s)
            G.nodes[node]['pos'] = pos
            
            # assign the cluster as a string of an integer in {1, ..., k}
            G.nodes[node]['cluster']=cluster_id_to_int_obj[cluster_id]

    # this cleans the final graph output
    # it makes all inter-cluster edges invisible
    #hide_inter_cluster_edges(G)

    # write out G into an svg, ready to be turned into a geojson
    write_dot(G,G_OUT_PATH)
    write_dot(H, H_OUT_PATH)
    #os.system('neato -Gforcelabels=false -Ecolor=grey  -Nshape=ellipse -n2 -Tsvg  tmp/dept_wise_graph.dot > prepolygon.svg')
    os.system(f'gvmap -e -c 3 -r 500  {G_OUT_PATH} > {DRAWING_DOT_PATH}')
    os.system(f'neato -Gforcelabels=false -Nlabel="" -Ecolor=grey28  -Nshape=ellipse -n2 -Tsvg  {DRAWING_DOT_PATH} > {DRAWING_PATH}')
    #os.system(f'neato -Gforcelabels=false -Nlabel="" -Ecolor=grey  -Nshape=ellipse -n2 -Tsvg  {G_OUT_PATH} > {DRAWING_PATH}')

if __name__ == "__main__":
    main()