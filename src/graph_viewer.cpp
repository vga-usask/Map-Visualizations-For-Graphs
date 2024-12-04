#include "RPGraphLayout.hpp"
#include "RPCommon.hpp"
#include "RPGraph.hpp"
#include "DataStore.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <cstdlib>
#include <cuda_runtime.h>
#include <thread>
#include <cuda_runtime_api.h>
#include <ctime>
#include <iomanip>
#include "RPGPUForceAtlas2.hpp"

auto program_start = std::chrono::high_resolution_clock::now();
bool png_flag = false;
bool svg_flag = false;
int rep_tune = 0;
typedef std::pair<double, double> Point_2;
std::vector<std::vector<int>> externa_polygonX(100000, std::vector<int>(500));

struct VectorPairEqual
{
    bool operator()(const std::vector<std::pair<double, double>> &lhs, const std::vector<std::pair<double, double>> &rhs) const
    {
        return lhs == rhs;
    }
};
void loadPolygonFile(const std::string &polygon_path, const std::string &node_polygon_path, DataStore &dataStore)
{
    auto start = std::chrono::high_resolution_clock::now();

    // Read the polygon file
    std::ifstream polygonFile(polygon_path);
    if (!polygonFile)
    {
        throw std::runtime_error("Unable to open polygon file: " + polygon_path);
    }
    std::cout << "Reading polygon file: " << polygon_path << std::endl;

    std::string line;
    int polygon_counter = 0;
    // Process the polygon file (polygonId -> points)
    while (std::getline(polygonFile, line))
    {
        std::istringstream iss(line);
        int polygonId;
        iss >> polygonId;
        std::cerr << "Points for Polygon " << polygonId << std::endl;
        std::vector<std::pair<double, double>> polygonPoints;
        std::string point;
        polygon_counter++;
        // Extract x, y points
        while (iss >> point)
        {
            double x, y;
            char comma;
            std::istringstream pointStream(point);
            if (pointStream >> x >> comma >> y && comma == ',')
            {
                polygonPoints.emplace_back(x, y);
                std::cerr << " " << x << "," << y;
            }
            else
            {
                std::cerr << "Error parsing point for polygonId: " << polygonId << std::endl;
                continue;
            }
        }

        dataStore.polygons[polygonId] = polygonPoints; // Store polygon points in the map

        // Convert to Point_2 for center/area calculation
        std::vector<Point_2> points;
        for (const auto &point : polygonPoints)
        {
            points.emplace_back(point.first, point.second);
        }

        // Calculate and store polygon center and area
        dataStore.polygonCenters[polygonId] = dataStore.calculatePolygonCenter(points);
        dataStore.polygonAreas[polygonId] = dataStore.calculatePolygonArea(points);
        if (dataStore.polygonAreas[polygonId] > dataStore.maximum_area)
        {
            dataStore.maximum_area = dataStore.polygonAreas[polygonId];
            dataStore.maximum_area_polygonId = polygonId;
        }
    }
    polygonFile.close();
    dataStore.number_of_polygons = polygon_counter;
    // Read the node-polygon association file
    std::ifstream nodePolygonFile(node_polygon_path);
    if (!nodePolygonFile)
    {
        throw std::runtime_error("Unable to open node-polygon file: " + node_polygon_path);
    }
    std::cout << "Reading node-polygon file: " << node_polygon_path << std::endl;

    // Process the node-polygon file (polygonId -> nodeIds)
    while (std::getline(nodePolygonFile, line))
    {
        std::istringstream iss(line);
        int polygonId;
        iss >> polygonId;
        std::cerr << "Adding nodes for polygon " << polygonId << std::endl;
        // Check if the polygonId exists in the polygon map
        if (dataStore.polygons.find(polygonId) == dataStore.polygons.end())
        {
            std::cerr << "PolygonId " << polygonId << " not found!" << std::endl;
            continue;
        }

        int nodeId;
        while (iss >> nodeId)
        {

            std::cerr << "  " << nodeId;
            // Assign polygon points to node and associate the node with a polygonId
            dataStore.nodeToPolygon[nodeId] = polygonId;
            dataStore.PolygonToNode[polygonId].push_back(nodeId);
        }
    }
    nodePolygonFile.close();

    // Measure the time taken
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to load polygon and node-polygon files: " << elapsed.count() << " seconds" << std::endl;
}

void loadEdgeListFile(const std::string &edgelist_path)
{
    auto start = std::chrono::high_resolution_clock::now();

    std::ifstream edgelist_file(edgelist_path);
    if (!edgelist_file.is_open())
    {
        std::cerr << "Error opening edge list file: " << edgelist_path << std::endl;
        return;
    }
    std::cout << "Reading edge list file: " << edgelist_path << std::endl;

    std::unordered_set<int> uniqueNodes;
    std::string line;
    while (std::getline(edgelist_file, line))
    {
        std::istringstream iss(line);
        int src, dst;
        if (!(iss >> src >> dst))
        {
            continue;
        }

        DataStore::getInstance().addEdge(src, dst);

        uniqueNodes.insert(src);
        uniqueNodes.insert(dst);
    }

    std::cout << "Finished reading edge list file. Number of edges: " << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to load edge list file: " << elapsed.count() << " seconds" << std::endl;
}

void createSubgraphsInMemory(
    const std::unordered_map<int, int> &nodeToPolygon,
    const std::vector<std::pair<int, int>> &edges,
    std::unordered_map<int, std::vector<std::pair<int, int>>> &subgraphs)
{
    std::cout << " Memory is loading..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    DataStore &dataStore = DataStore::getInstance();
    std::cout << " DataStore is loaded!" << std::endl;
    // Process edges to create subgraphs
    for (const auto &edge : edges)
    {
        int src = edge.first;
        int dst = edge.second;

        int itSrc = nodeToPolygon.at(src);
        int itDst = nodeToPolygon.at(dst);

        if (itSrc == itDst)
        {
            subgraphs[itSrc].push_back(edge);
        }
        else
        {
            externa_polygonX[src][itDst]++;
            externa_polygonX[dst][itSrc]++;
        }
    }
dataStore.externa_polygonX=externa_polygonX;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to create subgraphs in memory: " << elapsed.count() << " seconds" << std::endl;

    std::cerr << "createSubgraphsInMemory: " << std::endl;
}

// Utility function to calculate convex hull using Andrew's monotone chain algorithm
std::vector<Point_2> calculateConvexHull(std::vector<Point_2> points)
{
    if (points.size() <= 1)
        return points;

    std::sort(points.begin(), points.end());

    std::vector<Point_2> hull;

    for (int phase = 0; phase < 2; ++phase)
    {
        auto start_size = hull.size();
        for (const auto &point : points)
        {
            while (hull.size() >= start_size + 2 &&
                   (hull[hull.size() - 2].first - hull.back().first) * (point.second - hull.back().second) -
                           (hull[hull.size() - 2].second - hull.back().second) * (point.first - hull.back().first) <=
                       0)
            {
                hull.pop_back();
            }
            hull.push_back(point);
        }
        hull.pop_back();
        std::reverse(points.begin(), points.end());
    }

    return hull;
}

// Utility function to calculate polygon area using Shoelace formula
double calculatePolygonArea(const std::vector<Point_2> &points)
{
    double area = 0.0;
    int n = points.size();
    for (int i = 0; i < n; ++i)
    {
        int j = (i + 1) % n;
        area += points[i].first * points[j].second;
        area -= points[j].first * points[i].second;
    }
    return std::abs(area) / 2.0;
}

void processSubgraph(int polygonId, const std::vector<std::pair<int, int>> &subgraphEdges,
                     RPGraph::UGraph &full_graph, RPGraph::GraphLayout &full_layout, const std::unordered_map<int, int> &nodeToPolygon,
                     const std::unordered_map<int, double> &polygonAreas, bool approximate, bool strong_gravity,
                     float gravity, float scale, int max_iterations, int rep_tune, double &area_difference_percentage, int sizeofsubs, std::string out_filepath, bool cuda_requested, float area_propotion, DataStore &dataStore, int excep_polid, float excep_polid_scale, float progress)
{
    auto start = std::chrono::high_resolution_clock::now();
    double area_difference_percentage2 = area_difference_percentage;

   

    std::cout << "Processing subgraph: " << polygonId << std::endl;

    RPGraph::UGraph sub_graph(subgraphEdges);

    RPGraph::GraphLayout sub_layout(sub_graph, 3000, 3000, polygonId);
    std::cout << "FA2 for the subgraph " << polygonId << " is initializing... with " << sub_graph.num_nodes() << " nodes and " << sub_graph.num_edges() << "edges" << std::endl;

    float scale2 = scale;

    //********** */
    double area_difference_percentage2_maximum_happened = 0;
    float scale2_maximum_happened=0;
    //********** */

    int max_iterations2 = max_iterations;
    if (rep_tune > 0)
    {
        scale2 = (sqrt(sqrt(sqrt((polygonAreas.at(polygonId) / sub_graph.num_nodes())))));
        max_iterations2 = int(60 * sub_graph.num_nodes());
        if (max_iterations2>10000)
        max_iterations2=10000;
    }

    int rep = 0;

    do
    {
       

        //*****************************************************
        RPGraph::CUDAForceAtlas2 *fa2 = new RPGraph::CUDAForceAtlas2(sub_layout, approximate, strong_gravity, gravity, scale2, max_iterations2, externa_polygonX, polygonAreas.at(polygonId));
        for (int iteration = 1; iteration <= max_iterations2; ++iteration)
        {
            fa2->doStep(iteration);
        }

        fa2->sync_layout_void();
        fa2->freeGPUMemory();
        delete fa2;
        //*****************************************************
        if (sub_graph.num_nodes() <= 3)
            break;
        //***************************************************** */
        std::vector<Point_2> graph_nodes_positions;
        for (RPGraph::nid_t n = 0; n < sub_graph.num_nodes(); ++n)
        {
            graph_nodes_positions.push_back(Point_2(sub_layout.getX(n), sub_layout.getY(n)));
        }
        std::vector<Point_2> convex_hull = calculateConvexHull(graph_nodes_positions);
        double convex_hull_area = calculatePolygonArea(convex_hull);
        double original_polygon_area = polygonAreas.at(polygonId);
        area_difference_percentage2 = (  convex_hull_area/original_polygon_area);
        //**********************************************************************
        std::cout << std::endl
                  << "Area difference percentage for graph " << out_filepath << ": " << area_difference_percentage2 << " scale:" << scale2 << "Max Iterations " << max_iterations2 << " Progress: %" << progress<<"."<<rep  << std::endl;
        //*****************************************************
        /*if(area_difference_percentage2>area_difference_percentage2_maximum_happened)
        {
            area_difference_percentage2_maximum_happened=area_difference_percentage2;
            scale2_maximum_happened=scale2;
        }*/
        //*****************************************************
        if (rep_tune > 0)
        {

            if (area_difference_percentage2 < area_propotion)
            {
                scale2 = scale2*(1+(area_propotion-area_difference_percentage2));
                //scale2 = scale2*1.1;
            }
            else
            {
                break;
            }
        }
        //*****************************************************
        rep++;

    } while (area_difference_percentage2 < area_propotion);
    std::cout << " Area difference percentage after the upscaling: " << area_difference_percentage2 << " scale:" << scale2 << std::endl;

    std::vector<Point_2> graph_nodes_positions;

    for (RPGraph::nid_t i = 0; i < sub_graph.num_nodes(); ++i)
    {
        bool no_use_variable = sub_layout.move_node(i,sub_layout.getX(i), sub_layout.getY(i));
        if(!no_use_variable)
        std::cout << " Nodes are moving! "<<std::endl; 
        graph_nodes_positions.push_back(Point_2(sub_layout.getX(i), sub_layout.getY(i)));
        int ii = sub_graph.node_map_r.at(i);
        RPGraph::nid_t full_graph_nid = full_graph.node_map.at(ii);
        full_layout.setX(full_graph_nid, sub_layout.getX(i));
        full_layout.setY(full_graph_nid, sub_layout.getY(i));
    }
    std::vector<Point_2> convex_hull = calculateConvexHull(graph_nodes_positions);
    double convex_hull_area = calculatePolygonArea(convex_hull);
    double original_polygon_area = polygonAreas.at(polygonId);
    area_difference_percentage2 = (original_polygon_area / convex_hull_area);
    std::cout << " **Final Area difference percentage after the merging: " << area_difference_percentage2 << " scale:" << scale2 << std::endl;

    std::cerr << "Integration complete for subgraph " << polygonId << std::endl;
    std::cout << "Processing subgraph " << polygonId << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to process subgraph " << polygonId << ": " << elapsed.count() << " seconds" << std::endl;
}

std::string getUniqueFilename(const std::string &base_path, const std::string &extension)
{
    int counter = 1;
    std::string new_filename = base_path + "_final." + extension;
    std::srand(std::time(0));
    counter = 10 + std::rand() % 900;

    new_filename = base_path + "_final_" + std::to_string(counter) + "." + extension;

    return new_filename;
}

std::string extractFilename(const std::string &path)
{
    size_t lastSlash = path.find_last_of("/\\");
    size_t start = (lastSlash == std::string::npos) ? 0 : lastSlash + 1;

    size_t lastDot = path.find_last_of('.');
    size_t end = (lastDot == std::string::npos || lastDot < start) ? path.length() : lastDot;

    return path.substr(start, end - start);
}

int main(int argc, const char **argv)
{
    if (argc < 10)
    {
        std::cerr << "Usage: " << argv[0] << " <polygon_path> <edgelist_path> <output_path> <max_iterations> <num_screenshots> <strong_gravity> <scale> <gravity> <approximate> <rep_tune> <cornergravity> <Polygonsize> <png> <svg>" << std::endl;
        return EXIT_FAILURE;
    }

    double area_difference_percentage = 100;
    bool cuda_requested = false;
    int max_iterations = std::stoi(argv[4]);
    int num_screenshots = std::stoi(argv[5]);
    bool strong_gravity = std::string(argv[6]) == "sg";
    float scale = std::stof(argv[7]);
    float gravity = std::stof(argv[8]);
    bool approximate = std::string(argv[9]) == "approximate";
    std::string edgelist_path = argv[2];
    std::string out_path = argv[3];
    std::string out_format = "svg";
    rep_tune = std::stoi(argv[10]);
    double polysize = std::stoi(argv[12]);
    double cornergravity = std::stoi(argv[11]); //* polysize;
    std::string polygon_path = argv[1];

    png_flag = std::string(argv[13]) == "png";
    svg_flag = std::string(argv[14]) == "svg";
    cuda_requested = std::string(argv[15]) == "gpu";
    int real_number_nodes = std::stoi(argv[16]);
    std::string node_polygon_path = argv[17];
    float area_propotion = std::stof(argv[18]);
    int excep_polid = std::stoi(argv[19]);
    float excep_polid_scale = std::stof(argv[20]);
    std::string node_map_r_f = std::string(argv[21]);
    std::string dep_map_r = std::string(argv[22]);
    DataStore &dataStore = DataStore::getInstance();
    loadEdgeListFile(edgelist_path);

    loadPolygonFile(polygon_path, node_polygon_path, dataStore);
    std::cout << "loading polygons is done" << std::endl;
    const auto &edges = dataStore.edges;
    const auto &nodeToPolygon = dataStore.nodeToPolygon;
    const auto &polygonAreas = dataStore.polygonAreas;
    std::cout << " polygons variables is done" << std::endl;
    std::unordered_map<int, std::vector<std::pair<int, int>>> subgraphs;

    createSubgraphsInMemory(nodeToPolygon, edges, subgraphs);
    std::cout << " subgraph into memory is done" << std::endl;
    RPGraph::UGraph full_graph(edges);
    RPGraph::GraphLayout full_layout(full_graph);
    std::cout << "Full layout is created." << full_graph.num_nodes() << "**" << full_graph.num_edges() << std::endl;

    std::time_t now = std::time(nullptr);
    // Convert it to local time structure
    std::tm *local_time = std::localtime(&now);

    // Create a stringstream to format the date and time
    std::stringstream ss;
    ss << (local_time->tm_year + 1900) << "-"                                  // Year 
       << std::setw(2) << std::setfill('0') << (local_time->tm_mon + 1) << "-" // Month 
       << std::setw(2) << std::setfill('0') << local_time->tm_mday << "_"      // Day
       << std::setw(2) << std::setfill('0') << local_time->tm_hour << "-"      // Hour
       << std::setw(2) << std::setfill('0') << local_time->tm_min << "-"       // Minute
       << std::setw(2) << std::setfill('0') << local_time->tm_sec;             // Second
    std::srand(std::time(0));
    int randomNumber = 10000 + std::rand() % 90000;
    std::string out_filename = getUniqueFilename(extractFilename(edgelist_path) + "_Rep" + std::to_string(scale) + "_Gra" + std::to_string(gravity) + std::to_string(randomNumber), out_format);
    std::string out_filepath_png = out_path + "/" + "png/" + ss.str() + out_filename;
    std::string out_filepath_svg = out_path + '/' + "drawing.svg"; // + "/" + "svg/" + ss.str() + out_filename;
    std::string out_filepath = out_path + "/" + ss.str() + out_filename;

    std::cout << "Size of the subgraphs: " << subgraphs.size() << std::endl;
    int progress = subgraphs.size();
    for (const auto &[polygonId, subgraphEdges] : subgraphs)
    {

        progress--;

        std::cout << "size of subgraph " << polygonId << " is " << subgraphEdges.size() << std::endl;
        processSubgraph(polygonId, subgraphEdges, full_graph, full_layout, nodeToPolygon, polygonAreas, approximate, strong_gravity,
                        gravity, scale, max_iterations, rep_tune, area_difference_percentage, subgraphs.size(), out_filepath, cuda_requested, area_propotion, dataStore, excep_polid, excep_polid_scale, 100 - ((progress * 100) / subgraphs.size()));
        cudaDeviceReset();
    }

    std::cerr << "Writing into the file:" << std::endl;
    if (svg_flag)
    {
        full_layout.writeToSVG(out_filepath_svg, real_number_nodes, node_map_r_f, dep_map_r);
    }

    std::cerr << "File is saved!" << std::endl;

    auto program_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_elapsed = program_end - program_start;
    std::cout << "Total time taken: " << total_elapsed.count() << " seconds" << std::endl;

    return EXIT_SUCCESS;
}
