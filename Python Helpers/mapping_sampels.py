import pygraphviz as pgv
import os

def map_attributes_and_save(input_folder, output_folder):
    """
    Processes all DOT files in a folder, maps attributes, assigns unique IDs,
    and saves the updated graphs to a new folder.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Initialize mapping for attributes and global node IDs


    # Iterate through all DOT files in the input folder
    for file_name in os.listdir(input_folder):
        attribute_map = {}
        global_current_id = 1  # For mapping attribute values globally
        global_new_node_id = 1  # For assigning new node IDs globally
        if file_name.endswith(".dot"):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name)

            # Load the graph
            graph = pgv.AGraph(input_file_path)
            graph2 = pgv.AGraph(strict=False, directed=False)  # New graph

            # Create local mappings for this graph
            node_id_map = {}

            # Map node IDs and attributes
            for node in graph.nodes():
                # Assign a new global node ID
                node_id_map[node] = f"{global_new_node_id}"
                global_new_node_id += 1

                # Remove 'name' attribute if it exists
                #if 'name' in graph.get_node(node).attr:
                 #   del graph.get_node(node).attr['name']

                # Map specific attribute (e.g., 'department')
                attribute_value = graph.get_node(node).attr.get('department')  # Replace 'department' as needed
                if attribute_value:
                    # Map the attribute value globally
                    if attribute_value not in attribute_map:
                        attribute_map[attribute_value] = global_current_id
                        global_current_id += 1
                    # Update the attribute with the mapped value
                    graph.get_node(node).attr['department'] = str(attribute_map[attribute_value])

                # Add the updated node to the new graph
                graph2.add_node(node_id_map[node], **graph.get_node(node).attr)

            # Copy edges with updated node IDs
            for edge in graph.edges():
                src, dst = edge
                graph2.add_edge(node_id_map[src], node_id_map[dst], **edge.attr)

            # Save the modified graph to the output folder
            graph2.write(output_file_path)
            print(f"Processed and saved: {output_file_path}")

# Input and output folders
input_folder = "notmappedsampels"  # Folder containing the 100 DOT files
output_folder = "mappedsamples"  # Folder to save processed DOT files

map_attributes_and_save(input_folder, output_folder)
