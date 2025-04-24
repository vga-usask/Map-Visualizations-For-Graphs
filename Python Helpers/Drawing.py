'''
    this is the Drawing class for the objects we want to study
    it is esssentially a container for the the graph of nodes (G), 
    graph of clusters (H), and the drawing of those graphs (drawing), along
    with emthods on those objects

    look into cleaning up the polygon parsing part
'''
import networkx as nx
import pygraphviz as pgv
import svgelements as se
import numpy as np
import matplotlib.patches as mp

from scipy.spatial import ConvexHull
try:
    from scipy.spatial import QhullError
except ImportError:
    QhullError = Exception
from scipy.spatial.distance import pdist
from tqdm import tqdm

# processing -----

# extract EmbeddedPoint objects from the svg
def process_kmap_drawing_points(drawing_path, G):
    #print("   parsing points / nodes: ", end="")
    global_points = []
    with open(drawing_path, 'r') as drawing_file:
        file_lines = drawing_file.readlines()
        line_index = 0
        for line in tqdm(file_lines, desc="   parsing points / nodes: "):
            # check to ssee this is proper type
            type_index = 1
            type_element = line[type_index:type_index + 1]
            if type_element == "g":
                # check to see this is a point (node)
                class_index = line.index("class") + 7
                this_class = line[class_index:class_index + 4]
                if this_class == "node":
                    # we have a valid point, let's get the needed info
                    # and create an EmbeddedPoint
                    # fetch the node
                    this_node_start = file_lines[line_index + 1].index("<title>") + 7
                    this_node_end = file_lines[line_index + 1].index("</title>")
                    this_node = file_lines[line_index + 1][this_node_start:this_node_end]

                    try:
                        # fetch the x coord                        
                        this_x_start = file_lines[line_index + 2].index('cx="') + 4
                        this_x_end = file_lines[line_index + 2].index('" cy')
                        this_x = file_lines[line_index + 2][this_x_start:this_x_end]
                        # fetch the y coord
                        this_y_start = file_lines[line_index + 2].index('cy="') + 4
                        this_y_end = file_lines[line_index + 2].index('" rx')
                        this_y = file_lines[line_index + 2][this_y_start:this_y_end]
                    except:
                        # fetch the x coord                        
                        this_x_start = file_lines[line_index + 2].index('x="') + 3
                        this_x_end = file_lines[line_index + 2].index('" y')
                        this_x = file_lines[line_index + 2][this_x_start:this_x_end]
                        # fetch the y coord
                        this_y_start = file_lines[line_index + 2].index('y="') + 3
                        this_y_end = file_lines[line_index + 2].index('" font')
                        this_y = file_lines[line_index + 2][this_y_start:this_y_end]

                    #import pdb;
                    #pdb.set_trace()

                    # get the label and cluster_id from G
                    # THIS IS VERY KMAP SPECIFIC FIX LATER
                    this_label = G.nodes[this_node]['name']
                    this_cluster_id = G.nodes[this_node]['department']
                    # create an EmbeddedPoint obj and append it to the list
                    global_points.append(EmbeddedPoint(this_x, this_y, this_node, this_label, this_cluster_id))
            line_index += 1
    return global_points
                    
# this turns the svg into a collection of polygons that we can
# work with. it is pretty specific to the kmap data
def process_kmap_drawing_polygons(drawing_path, global_points):
    # fetch the polygons from the svg
    these_polygons = []
    with open(drawing_path, 'r') as drawing_file:
        # we want to loop through and find all polygons that are
        # filled in -- these correspond to departments in the kmap
        file_lines = drawing_file.readlines()
        for line in tqdm(file_lines, "   parsing polygons and interior points: "):
            # check to see if this is a polygon of interest
            type_index = 1
            type_element = line[type_index:8]
            if type_element == "polygon":
                fill_index = line.index("fill") + 6
                fill = line[fill_index:fill_index + 4]
                stroke_index = line.index("stroke") + 8
                stroke = line[stroke_index: stroke_index + 11]
                if fill != "none" and stroke != "transparent":
                    # at this point we have a valid polygon, parse it then add it
                    point_string_index = line.index("points") + 8
                    point_string = line[point_string_index:-4]
                    this_polygon = EmbeddedPolygon(se.Polygon(point_string), global_points)
                    if not this_polygon.is_lake():
                        these_polygons.append(this_polygon)
            else:
                continue
    return these_polygons

# get the interior points from global_points for this polygon
# inputs are this_mp_polygon, which is the matplotlib polygon
# for this polygon object, and global_points, which is a list
# of EmbeddedPoint objects
def get_interior_points(this_mp_polygon, global_points):
    # naive approach?
    #print(f"parsing interior points for {this_mp_polygon}")
    interior_points = []
    for point in global_points:
        this_pt_coord = point.get_coord_tuple()
        if this_mp_polygon.contains_points([this_pt_coord]):
            interior_points.append(point)
    return interior_points

# this gets the cluster_id from the cluster; there might be a 
# quicker way to get this from H, but since H is a bit unreliable
# for graph layouts in general, we use this implementation
def retrieve_cluster_id(interior_points):
    # this maps all cluster_ids from interior regions of the polygon
    # to their count; by C_1 this should be only one cluster_id, but this
    # serves as a measurement for our first hard constraint -- that
    # all points are contained within the cluster's designated polygon
    these_candidate_ids = {}
    for point in interior_points:
        this_cluster_id = point.get_cluster_id()
        if this_cluster_id in these_candidate_ids:
            these_candidate_ids[this_cluster_id] += 1
        else:
            these_candidate_ids[this_cluster_id] = 1
    # now sort for the proper cluster_id
    sorted_cluster_ids = sorted(these_candidate_ids.items(), key=lambda item: item[1], reverse=True)
    #print("HERE")
    #print(interior_points)
    #print(sorted_cluster_ids)
    #print("BYE")
    
    return sorted_cluster_ids[0][0], these_candidate_ids

# classes -----

# this corresponds ot a single polygon in the drawing, along with the corresponding
# graph data
class EmbeddedPolygon:
    def __init__(self, polygon, global_points):
        self._svg_polygon = polygon
        self._mp_polygon = mp.Polygon(polygon.points, closed=True)
        self._interior_points = get_interior_points(self._mp_polygon, global_points)
        if len(self._interior_points) > 0:
            self._cluster_id, self._interior_cluster_ids_w_count = retrieve_cluster_id(self._interior_points)
        else:
            self._cluster_id = None
            self._interior_cluster_ids_w_count = None
   
    # check if the polygon is a lake ie it has no nodes
    # for now this just filters it out for the measure, but work on this later
    def is_lake(self):
        if len(self._interior_points) == 0:
            return True
        return False
   
    # note that this returns the cluster by looking at the interior points in the drawing and
    # then re-mapping them to G gathering the information straight from a dot file would be 
    # quicker, but it's not contained in every output
    def get_cluster_id(self):
        return self._cluster_id
    
    # returns a collection of mpl Point objects that are the polygons vertices
    def get_vertices(self):
        return self._mp_polygon.get_xy()

    # returns a scalar of the polygon's perimeter
    def get_perimeter(self):
        perimeter_points = self._mp_polygon.get_xy()
        all_distances = pdist(perimeter_points, metric="euclidean")
        # we need to jump around the array to retrieve the correct distances
        this_adjacency_jump = len(perimeter_points) - 1
        i = 0
        perimeter = 0.0
        while (i < len(all_distances)):
            this_dist = all_distances[i]
            perimeter += this_dist
            i += this_adjacency_jump
            this_adjacency_jump -= 1
        return perimeter

    # returns a scalar of the polygon's area (shoelace)
    def get_area(self):
        vertices = self.get_vertices()
        vertices_array = np.stack(vertices, axis=1)
        x_vector = vertices_array[0]
        y_vector = vertices_array[1]
        shoelace_area = 0.5*np.abs(np.dot(x_vector,np.roll(y_vector,1))-np.dot(y_vector,np.roll(x_vector,1)))
        return shoelace_area

    # returns a scipy ConevHull object of the vertices of the polygon
    def get_convex_hull_vertices(self):
        return ConvexHull(self._mp_polygon.get_xy())
    
    def get_convex_hull_perimeter(self):
        # this is generalized area, so in 2d it returns the perimeter
        return ConvexHull(self._mp_polygon.get_xy()).area
    
    def get_convex_hull_area(self):
        return ConvexHull(self._mp_polygon.get_xy()).volume

    # returns a list of tuples of the coordinates of all interior points in the polygon
    # (interior meaning within the polygons boundary in the plane) 
    def get_interior_points_coordinates(self):
        # again, look into numpy stuff
        interior_points_coordinates = []
        for point in self._interior_points:
            this_coord_tuple = point.get_coord_tuple()
            interior_points_coordinates.append(this_coord_tuple)
        return interior_points_coordinates

    # returns how many cluster ids (properties of nodes) occur in the interior points of the 
    # the polygon
    def get_number_of_cluster_ids(self):
        return len(self._interior_cluster_ids_w_count)

    # ask about changing the measure to perimeter instead of area
    # it seems to avoid the need to take area-intersections, and will be bound
    # to [0,1] since all points contained in the interior and triangle inequality
    def get_interior_subgraph_perimeter(self):
        if len(self._interior_points) < 3:
            return -1
        else:
            try:
                return ConvexHull(self.get_interior_points_coordinates()).area
            except QhullError:
                return 0.0

# this corresponds to a node of G's embedding on the plane
class EmbeddedPoint:
    def __init__(self, x_str, y_str, node_str, label_str, cluster_str):
        self._x = x_str
        self._y = y_str
        self._node = node_str
        self._label = label_str
        self._cluster_id = cluster_str

    def get_label(self):
        return self._label

    def get_coord_tuple(self):
        return (float(self._x), float(self._y))
    
    # this is still ultimately induced by the plane drawing not the graph
    # update this if that changes
    def get_cluster_id(self):
        return self._cluster_id
    
# this corresponds to any point that is drawn as a boundary point for the polygons. They aren't
# included in the graphs, but serve an important purpose in finding when polygons are adjacent
class BorderPoint:
    def __init__(self, coord_array, curr_polygon):
        self._coords = coord_array
        self._contained_in = {curr_polygon}

    def include_in(self, this_polygon):
        self._contained_in.add(this_polygon)

    def get_all_polygons_contained_in(self):
        return self._contained_in

# this is the main class we use to measure the drawings
class Drawing:
    def __init__(self, name_str, G_path, H_path, drawing_path):
        self._name = name_str
        self._G = nx.Graph(pgv.AGraph(G_path))
        self._H = nx.Graph(pgv.AGraph(H_path))
        self._points = process_kmap_drawing_points(drawing_path, self._G)
        self._polygons = process_kmap_drawing_polygons(drawing_path, self._points)
        
    # returns the name
    def get_name(self):
        return self._name   

    # returns the list of EmbeddedPolygon objects in this drawing
    def get_polygons(self):
        return self._polygons
    
    # returns a set of the clusters in the drawing 
    def get_clusters(self):
        return self._H.nodes()

    # returns a union of all vertices of all polygons in the drawing
    def get_entire_drawing_vertices(self):
        # try to get this to a numpy array creation
        polygon_tuples = set()
        for polygon in self._polygons:
            these_vertices = polygon.get_vertices()
            for vertex in these_vertices:
                polygon_tuples.add(tuple(vertex))
        return list(polygon_tuples)
    
    # returns the convex hull of all polygons (of all polygons' BorderPoints)
    def get_entire_drawing_convex_hull(self):
        return ConvexHull(self.get_entire_drawing_vertices())

    def get_entire_convex_hull_area(self):
        return self.get_entire_drawing_convex_hull().volume
    
    # sums up all areas and returns it; use this if you want the map's area but not that of the lakes / gaps
    def get_cumulative_polygon_areas(self):
        area_sum = 0.0
        for polygon in self._polygons:
            area_sum += polygon.get_area()
        return area_sum
    
    # maps a cluster to its size (number of nodes) and returns the dict
    def get_cluster_to_size_mapping(self):
        cluster_to_size = {}
        for node in self._H.nodes():
            cluster_to_size[node] = int(self._H.nodes[node]['weight'])
        return cluster_to_size

    # probably remove this later, but it might help speed up / ensure the M_3 calculation works
    def get_cumulative_cluster_sizes(self):
        cluster_size_sum = 0.0
        cluster_to_size = self.get_cluster_to_size_mapping()
        for cluster in cluster_to_size:
            cluster_size_sum += cluster_to_size[cluster]
        return cluster_size_sum
    
    # returns the size of a given cluster
    def get_cluster_size(self, cluster):
        return self.get_cluster_to_size_mapping()[cluster]
    
    # border points and neighbor stuff
    def _get_border_points(self):
        border_points = []
        coords_already_visited = {}
        for polygon in self._polygons:
            these_vertices = polygon.get_vertices()
            for vertex_array in these_vertices:
                vertex = tuple(vertex_array)
                if vertex in coords_already_visited:
                    this_border_point = coords_already_visited[vertex]
                    this_border_point.include_in(polygon)
                else:
                    this_border_point = BorderPoint(vertex, polygon)
                    coords_already_visited[vertex] = this_border_point
                    border_points.append(this_border_point)
        return border_points
        
    # go through all the border points, and since a border point being in two polygons 
    # implies they are neighbors, return a dictionary for each polygon to its neighbors
    def get_neighbor_mapping(self):
        border_points = self._get_border_points()
        polygon_to_neighbors = {}
        for point in border_points:
            these_polygons_containing_point = point.get_all_polygons_contained_in()
            indexing_copy = these_polygons_containing_point.copy()
            while len(indexing_copy) > 0:
                this_polygon = indexing_copy.pop()
                these_neighbors = these_polygons_containing_point.difference({this_polygon})
                if this_polygon in polygon_to_neighbors:
                    polygon_to_neighbors[this_polygon].update(these_neighbors)
                else:
                    polygon_to_neighbors[this_polygon] = these_neighbors.copy()
        return polygon_to_neighbors

    # go through H and return a dict of each node mapped to its neighbor
    def get_H_neighbor_mapping(self):
        node_to_neighbors = {}
        for edge in self._H.edges():
            node_A = edge[0]
            node_B = edge[1]
            if node_A in node_to_neighbors:
                node_to_neighbors[node_A].add(node_B)
            else:
                node_to_neighbors[node_A] = {node_B}
            if node_B in node_to_neighbors:
                node_to_neighbors[node_B].add(node_A)
            else:
                node_to_neighbors[node_B] = {node_A}
        # not all edges have nodes, so add the ones that don't here
        for node in self._H.nodes():
            if node not in node_to_neighbors:
                node_to_neighbors[node] = set()
        return node_to_neighbors
