import pygraphviz as pgv

# Load the DOT file using pygraphviz
graph = pgv.AGraph('notG.dot')
graph2 = pgv.AGraph(strict=False, directed=False)  # Explicitly set the new graph as undirected

# Initialize dictionaries for mapping unique attribute values and node IDs
attribute_map = {}
current_id = 1
node_id_map = {}
new_node_id = 1

# Iterate over each node to map attributes and assign new node IDs
for node in graph.nodes():
    # Map and set the new node ID
    node_id_map[node] = f"{new_node_id}"  # Assign a new unique node ID
    new_node_id += 1

    # Remove 'name' attribute if it exists
    if 'name' in graph.get_node(node).attr:
        del graph.get_node(node).attr['name']

    # Map the specific attribute, e.g., 'department'
    attribute_value = graph.get_node(node).attr.get('department')  # Replace 'department' with your attribute name
    
    if attribute_value:
        # Check if the attribute value is already in the map
        if attribute_value not in attribute_map:
            # Assign a new ID if it's a new value
            attribute_map[attribute_value] = current_id
            current_id += 1

        # Set the mapped ID for the attribute
        graph.get_node(node).attr['department'] = str(attribute_map[attribute_value])

    # Add the new node to the new graph with mapped attributes
    graph2.add_node(node_id_map[node], **graph.get_node(node).attr)

# Copy edges with updated node IDs
for edge in graph.edges():
    src, dst = edge
    graph2.add_edge(node_id_map[src], node_id_map[dst], **edge.attr)

# Write the modified graph back to DOT format
graph2.write('G.dot')
