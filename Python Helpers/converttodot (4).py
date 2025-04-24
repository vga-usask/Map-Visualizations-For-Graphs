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
            dept_mapped=dept_map[dept_id]
            if node_id in node_map:
                department_nodes[node_map[node_id]].append(dept_mapped)
            else:
                node_map[node_id]=max(node_map.values())+1
                department_nodes[node_map[node_id]].append(dept_mapped)
                print(node_map[node_id])
    max_department_id = dept_index
    return department_nodes,max_department_id

# Generate G.dot file
def generate_dot_file(edges, department_nodes, output_path,max_department_id):

    # Build a set of unique nodes from the edge list
    unique_nodes = set()
    
    dep_size={}
    for src, dst in edges:
        unique_nodes.add(src)
        unique_nodes.add(dst)

    # Determine the highest department ID in the input data
    
    fake_department_id = max_department_id  # This will be the fake department ID

    # Prepare department connections for each node in the edge list
    main_department = {}
    node_dept_count = defaultdict(lambda: defaultdict(int))
    for node, depts in department_nodes.items():
        for dept in depts:
            node_dept_count[node][dept] += 1

    # Determine main department for each node
    for node in unique_nodes:
        if node in node_dept_count and node_dept_count[node]:
            # Node has department(s)
            
            main_dept = max(node_dept_count[node], key=node_dept_count[node].get)
            main_department[node] = main_dept
            dep_size[main_dept] +=1
        else:
            # Assign the same fake department ID to all nodes without department or neighbor
            main_department[node] = fake_department_id
            print(f"Assigned fake department {fake_department_id} to node {node}.")    

    # Write to G.dot, only including nodes in the edge list
    with open(output_path, 'w') as f:
        f.write("strict graph \"\" {\nnode [label=\"\\N\"];\n")
        # Write nodes with mapped indices
        for node in unique_nodes:
            dept_id = main_department[node]
            f.write(f'{node} [department={dept_id},\nname="{node}",\nweight=0];\n')
        # Write edges, converting src and dst to their mapped indices
        for src, dst in edges:
            f.write(f'{src} -- {dst} [type="t,l",\nweight=1];\n')
        
        f.write("}")


# Usage
edge_list_file = 'author-author.txt'
department_node_file = 'author-conf.txt'
output_dot_file = edge_list_file+'G_limited_not_final_mapped.dot'
edges,node_map = read_edge_list(edge_list_file)
department_nodes,max_department_id = read_department_nodes(department_node_file,node_map)
generate_dot_file(edges, department_nodes, output_dot_file,max_department_id)

