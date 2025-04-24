from collections import defaultdict
  # Mapping of {index: node_id}
  # Start department IDs from 1
      # Start node mapping index from 1
# Read edge list with src,dst format (comma-separated)
def read_edge_list(file_path):
    node_map = {} 
    node_index = 1 
    edges = []
    with open(file_path, 'r') as f:
        for line in f:
            src, dst = map(int, line.strip().split(','))
            if src not in node_map:
                node_map[src]=node_index
                node_index +=1
            if dst not in node_map:
                node_map[dst]=node_index
                node_index +=1
            src_mapped=node_map[src]
            dst_mapped=node_map[dst]
            edges.append((src_mapped, dst_mapped))
    return edges,node_map

# Read department-node associations with node,dep format (comma-separated)
def read_department_nodes(file_path,node_map):
    dept_map = {}
    dept_index = 1
    department_nodes = defaultdict(list)
    with open(file_path, 'r') as f:
        for line in f:
            node_id, dept_id = map(int, line.strip().split(','))
            
            if dept_id not in dept_map:
                dept_map[dept_id]=dept_index
                dept_index +=1
                print(dept_id)
                print(dept_map[dept_id])
            dept_mapped=dept_map[dept_id]
            if node_id in node_map:
                department_nodes[node_map[node_id]].append(dept_mapped)
            else:
                node_map[node_id]=max(node_map.values())+1
                department_nodes[node_map[node_id]].append(dept_mapped)
                #print(node_map[node_id])
                
    max_department_id = dept_index
    return department_nodes,max_department_id

# Generate G.dot file
def generate_dot_file(edges, department_nodes, output_path, max_department_id):
    unique_nodes = {src for src, dst in edges}.union({dst for src, dst in edges})
    dep_size = defaultdict(int)  # Track the number of nodes per department
    main_department = {}  # Store the main department for each node
    node_dept_count = defaultdict(lambda: defaultdict(int))
    current_fake_dept = max_department_id  # Start with the highest department ID

    # Count departments for each node
    for node, depts in department_nodes.items():
        for dept in depts:
            node_dept_count[node][dept] += 1

    # Assign main departments, limiting each department to 500 nodes
    for node in unique_nodes:
        if node in node_dept_count and node_dept_count[node]:
            # Get the department with the highest count
            sorted_depts = sorted(node_dept_count[node].items(), key=lambda x: (-x[1], x[0]))
            for dept, _ in sorted_depts:
                if dep_size[dept] < 100000:
                    main_department[node] = dept
                    dep_size[dept] += 1
                    break

    # Write the DOT file
    with open(output_path, 'w') as f:
        f.write("strict graph \"\" {\n  node [label=\"\\N\"];\n")
        # Write nodes with their departments
        for node in unique_nodes:
            if node in main_department:
                dept_id = main_department[node]
                f.write(f'  {node} [department={dept_id}, name="{node}", weight=0];\n')
        # Write edges
        for src, dst in edges:
            if src in main_department and dst in main_department:
                f.write(f'  {src} -- {dst} [type="t,l", weight=1];\n')
        f.write("}\n")



# Usage
edge_list_file = 'author-author.txt'
department_node_file = 'author-conf.txt'
output_dot_file = edge_list_file+'G_limited_not_final_mapped.dot'
edges,node_map = read_edge_list(edge_list_file)
department_nodes,max_department_id = read_department_nodes(department_node_file,node_map)
generate_dot_file(edges, department_nodes, output_dot_file,max_department_id)

