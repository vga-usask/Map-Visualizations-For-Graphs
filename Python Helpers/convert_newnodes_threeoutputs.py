import re
import xml.etree.ElementTree as ET
from shapely.geometry import Point, Polygon
import argparse
from pathlib import Path

# Define function to create node name file
def create_node_name_file(graph_file_path):
    node_label_pattern = re.compile(r'(\d+)\s+\[department=(\d+),\s*name="([^"]+)"')
    node_names = {}

    with open(graph_file_path, 'r') as file:
        current_line = ""
        for line in file:
            current_line += line.strip()
            node_match = node_label_pattern.search(current_line)
            if node_match:
                node_id = int(node_match.group(1))
                node_name = int(node_match.group(3))
                node_names[node_id] = node_name
                current_line = ""

    with open('node_to_name.txt', 'w') as output_file:
        for node_id, node_name in node_names.items():
            output_file.write(f'{node_id}, {node_name}\n')

    return node_names

def read_graph_file(graph_file_path):
    edge_counter = 0
    graph = {}
    edges = {}

    node_pattern = re.compile(r'(\d+)\s+\[department=(\d+),')
    edge_pattern = re.compile(r'"?(\d+)"?\s+--\s+"?(\d+)"?')

    with open(graph_file_path, 'r') as file:
        for line in file:
            node_match = node_pattern.search(line)
            if node_match:
                node_id = int(node_match.group(1))
                department = int(node_match.group(2))
                graph[node_id] = department
                if node_id not in edges:
                    edges[node_id] = set()

    with open(graph_file_path, 'r') as file:
        for line in file:
            edge_match = edge_pattern.search(line)
            if edge_match:
                edge_counter += 1
                node1 = int(edge_match.group(1))
                node2 = int(edge_match.group(2))
                edges[node1].add(node2)
                edges[node2].add(node1)

    return graph, edges

def extract_bracket_blocks(file_path):
    blocks = []
    current_block = []
    inside_block = False

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if inside_block:
                current_block.append(line)
                if ']' in line:
                    inside_block = False
                    current_block[-1] = current_block[-1].split(']', 1)[0]
                    blocks.append(' '.join(current_block))
                    current_block = []
            if '[' in line and not inside_block:
                inside_block = True
                current_block.append(line.split('[', 1)[1])

    return blocks

def parse_block_to_dict(block):
    data_dict = {}
    for item in block.split(', '):
        if '=' in item:
            key, value = item.split('=', 1)
            value = value.strip('"').strip()

            if key == 'pos':
                value = value.replace(' ', ',')
                value = [float(v) for v in value.split(',')]
            else:
                try:
                    if '.' in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass
            data_dict[key.strip()] = value
    return data_dict

def read_h_embedded_file(embedded_file_path):
    departments_points = {}
    bracket_blocks = extract_bracket_blocks(embedded_file_path)
    
    

    for block in bracket_blocks:
        attributes = parse_block_to_dict(block)
        
        if 'department' in attributes and 'pos' in attributes and 'width' in attributes and 'height' in attributes:
            department_original_id = int(attributes['department'])

            
            department_id = department_original_id
            
            pos_x, pos_y = attributes['pos']
            width = float(attributes['width']) * 72*(2/3)
            height = float(attributes['height']) * 72*(2/3)
            
            half_width = width / 2
            half_height = height / 2
            top_left = (pos_x - half_width, pos_y + half_height)
            top_right = (pos_x + half_width, pos_y + half_height)
            bottom_left = (pos_x - half_width, pos_y - half_height)
            bottom_right = (pos_x + half_width, pos_y - half_height)
            
            departments_points[department_id] = [top_left, top_right, bottom_right, bottom_left]

    return departments_points

def count_nodes_edges(edges):
    num_nodes = len(edges)
    num_edges = sum(len(neighbours) for neighbours in edges.values()) // 2
    return num_nodes, num_edges

def add_extra_nodes(graph, edges):
    max_node_id = max(graph.keys())
    new_edges = {}
    new_graph = {}
    next_node_id = max_node_id + 1

    for node, department in graph.items():
        if not edges[node] or all(graph[neighbour] != department for neighbour in edges[node]):
            new_graph[next_node_id] = department
            new_edges[next_node_id] = {node}
            edges[node].add(next_node_id)
            next_node_id += 1

    graph.update(new_graph)
    edges.update(new_edges)

    return graph, edges, max_node_id

def write_edgelist_file(edges, file_path='edgelist.txt'):
    with open(file_path, 'w') as f:
        for src, neighbours in edges.items():
            for dst in neighbours:
                if src < dst:
                    f.write(f"{src} {dst}\n")

def write_polygon_points_file(departments_points, file_path='polygon_points.txt'):
    with open(file_path, 'w') as f:
        for polygon_id, points in departments_points.items():
            points_str = ' '.join([f"{x},{y}" for (x, y) in points])
            f.write(f"{polygon_id} {points_str}\n")


def write_node_polygon_mapping_file(graph, file_path):
    polygon_to_nodes = {}
    for node, polygon_id in graph.items():
        if polygon_id not in polygon_to_nodes:
            polygon_to_nodes[polygon_id] = []
        polygon_to_nodes[polygon_id].append(node)

    with open(file_path, 'w') as f:
        for polygon_id, nodes in polygon_to_nodes.items():
            nodes_str = ' '.join(str(node) for node in nodes)
            f.write(f"{polygon_id} {nodes_str}\n")

def extract_polylines_with_ids(svg_file, output_file):
    tree = ET.parse(svg_file)
    root = tree.getroot()
    polylines = []
    texts = []

    for elem in root.iter():
        if elem.tag.endswith('polyline'):
            points = elem.attrib.get('points', '').strip()
            if points:
                points_list = [(float(x), float(y)) for x, y in (point.split(',') for point in points.split())]
                if points_list[0] != points_list[-1]:  # Check if the polyline is not closed
                    points_list.append(points_list[0])  # Close the polyline
                polylines.append(points_list)

    # Rest of the function remains unchanged...


    for elem in root.iter():
        if elem.tag.endswith('text'):
            try:
                x = float(elem.attrib.get('x', 0))
                y = float(elem.attrib.get('y', 0))
                text_content = elem.text.strip() if elem.text else ''
                texts.append((text_content, x, y))
            except ValueError:
                print(f"Skipping invalid text coordinates: {elem.attrib.get('x', 0)}, {elem.attrib.get('y', 0)}")

    polyline_data = []

    for points_list in polylines:
        polygon = Polygon(points_list)
        polyline_id = None
        found_text = False

        for text, tx, ty in texts:
            text_point = Point(tx, ty)
            if polygon.contains(text_point):
                polyline_id = text
                found_text = True
                break

        if polyline_id:
            polyline_data.append(f"{polyline_id} {' '.join([f'{x},{y}' for x, y in points_list])}")

    with open(output_file, 'w') as f:
        for line in polyline_data:
            f.write(line + '\n')

parser = argparse.ArgumentParser()
parser.add_argument('input_dir', type=str, help='path of input directory containing .dot files')
parser.add_argument('input_dir2', type=str, help='path of input directory containing .dot files')
parser.add_argument('output_path', nargs='?', type=str, help='path for output directory')
args = parser.parse_args()

max_node_ids = {}

for dot_file in Path(args.input_dir).glob("G.dot"):
    graph_file_path = Path(args.input_dir + '/G.dot')
    graph, edges = read_graph_file(graph_file_path)
    
    departments_points = read_h_embedded_file(Path(args.input_dir + '/' + dot_file.stem + '_H_embedded.dot'))
    graph, edges, max_node_id = add_extra_nodes(graph, edges)
    max_node_ids[dot_file.stem] = max_node_id
    
    write_edgelist_file(edges, Path(args.output_path + '/edgelist.txt'))
    write_polygon_points_file(departments_points, Path(args.output_path + '/polygon_points.txt'))
    write_node_polygon_mapping_file(graph, Path(args.output_path + '/node_polygon_mapping.txt'))
    
    extract_polylines_with_ids(args.input_dir2 + '/' + dot_file.stem + '_H_embedded.svg', Path(args.output_path + '/_polylines_output.txt'))
    print(f"Total number of nodes: {len(edges)}")
    print(f"Total number of edges: {sum(len(neighbors) for neighbors in edges.values()) // 2}")
    print("Files 'edgelist.txt', 'polygon_points.txt', and 'node_polygon_mapping.txt' have been created.")

#with open(Path(args.output_path + '/runcode.sh'), 'w') as f:
 #   f.write(f"mkdir out/{args.input_dir}\n")
  #  for filename, max_node_id in max_node_ids.items():
  #      f.write(f"./graph_viewer ./{args.input_dir}/voronoipoints/{filename}_polygon_points.txt ./{args.input_dir}/outputs_convert/{filename}_edgelist.txt ./{args.input_dir}/{filename} 2000 1 wg 1 1 approximate 1 1 1 png svg gpu {max_node_id} ./{args.input_dir}/outputs_convert/{filename}_node_polygon_mapping.txt 5 13 5 ./{args.input_dir}/outputs_convert/{filename}_node_map.txt ./{args.input_dir}/outputs_convert/{filename}_dep_map.txt\n")
