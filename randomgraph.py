import random
from collections import defaultdict
from pydot import graph_from_dot_file, Dot, Node, Edge
import pandas as pd


def count_nodes_and_edges(graph):
    """Counts the number of nodes and edges in a pydot graph."""
    nodes = set()
    edges = 0
    for node in graph.get_nodes():
        nodes.add(node.get_name())
    for edge in graph.get_edges():
        edges += 1
    return len(nodes), edges


def generate_smaller_graphs(input_dot_file, output_prefix, sizes_per_department, output_excel):
    # Load the original graph
    graphs = graph_from_dot_file(input_dot_file)
    original_graph = graphs[0] if graphs else None

    if not original_graph:
        print("Could not load the DOT file.")
        return

    # Group nodes by department
    department_nodes = defaultdict(list)
    node_attributes = {node.get_name(): node.get_attributes() for node in original_graph.get_nodes()}
    for node_name, attributes in node_attributes.items():
        department = attributes.get("department")
        if department:
            department_nodes[department].append(node_name)

    # Data for Excel
    data = []

    # Precompute edges for efficiency
    original_edges = [(edge.get_source(), edge.get_destination(), edge.get_attributes())
                      for edge in original_graph.get_edges()]

    # Create smaller graphs for each specified size
    for size in sizes_per_department:
        for i in range(10):  # Create 10 graphs per size
            sampled_nodes = set()
            for department, nodes in department_nodes.items():
                sample_size = min(size, len(nodes))
                sampled_nodes.update(random.sample(nodes, sample_size))

            # Create a new graph with sampled nodes
            new_graph = Dot(graph_type="graph")
            for node in sampled_nodes:
                new_node = Node(name=node, **node_attributes[node])
                new_graph.add_node(new_node)

            # Add edges only if both nodes are in the sampled set
            for src, dst, attrs in original_edges:
                if src in sampled_nodes and dst in sampled_nodes:
                    new_graph.add_edge(Edge(src, dst, **attrs))

            # Count nodes and edges in the smaller graph
            num_nodes, num_edges = count_nodes_and_edges(new_graph)

            # Save the smaller graph to a new DOT file
            sample_name = f"{output_prefix}_size_{size}_graph_{i + 1}"
            output_file = f"{sample_name}.dot"
            new_graph.write(output_file)

            # Add data to Excel
            data.append({
                "SampleName": sample_name,
                "NumberOfEdges": num_edges,
                "NumberOfNodes": num_nodes
            })
            print(f"Created: {output_file} - Nodes: {num_nodes}, Edges: {num_edges}")

    # Save the data to an Excel file
    df = pd.DataFrame(data)
    df.to_excel(output_excel, index=False)
    print(f"Excel file created: {output_excel}")


# Input DOT file and sizes
input_dot_file = "G.dot"
output_prefix = "notmappedsampels/small_graph"
sizes_per_department = [100, 300, 500, 700, 900]
output_excel = "notmappedsampels/graph_samples_summary.xlsx"

generate_smaller_graphs(input_dot_file, output_prefix, sizes_per_department, output_excel)
