import pygraphviz as pgv
import os

def map_attributes_and_save(input_folder):
    """
    Processes all DOT files in a folder, maps attributes, assigns unique IDs,
    and organizes outputs in specific folder structures.
    """
    # Iterate through all DOT files in the input folder
    for file_name in os.listdir(input_folder):
        attribute_map = {}
        global_current_id = 1  # For mapping attribute values globally
        global_new_node_id = 1  # For assigning new node IDs globally

        if file_name.endswith(".dot"):
            base_name = file_name[:-4]  # Remove '.dot' from the filename
            main_folder = os.path.join("mapped_perpolygon_samples", base_name)
            outputs_convert_folder = os.path.join(main_folder, "outputs_convert")

            # Create the main folder and outputs_convert folder
            os.makedirs(outputs_convert_folder, exist_ok=True)

            # Create six additional folders beside the main folder
            additional_folders = [
                f"base1_{base_name}",
                f"base2_{base_name}",
                f"base3_{base_name}",
                f"our1_{base_name}",
                f"our2_{base_name}",
                f"our3_{base_name}",
                f"our4_{base_name}",
            ]
            for folder in additional_folders:
                os.makedirs(os.path.join("mapped_perpolygon_samples", folder), exist_ok=True)

            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(main_folder, "G.dot")

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

            # Save the modified graph to the G.dot file
            graph2.write(output_file_path)
            print(f"Processed and saved: {output_file_path}")
            os.system(f'python3 10_generate_H.py {main_folder} {main_folder}')
            os.system(f'python3 30_drawing_G_with_H_kmp.py {main_folder}/G.dot {main_folder}/G_H_embedded.dot {main_folder}/outputs_convert/')
            os.system(f'python3 convert_newnodes_threeoutputs.py {main_folder} {main_folder} {main_folder}/outputs_convert')
            os.system(f'python3 polylineextraction.py {main_folder}/outputs_convert/drawing.svg {main_folder}/outputs_convert/drawing.dot {main_folder}/outputs_convert/polyline_point.txt')
            os.system(f'python3 voronoi_weighted.py {main_folder} {main_folder}/outputs_convert {main_folder}/outputs_convert')

# Input folder
input_folder = "notmapped_perpolygon_samples"  # Folder containing the DOT files

map_attributes_and_save(input_folder)

