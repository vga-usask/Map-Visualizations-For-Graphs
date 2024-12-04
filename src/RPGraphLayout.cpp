#include "RPGraphLayout.hpp"
#include "RPCommon.hpp" // Include for singleton access
#include <fstream>
#include <cmath>
#include <limits>
#include <iostream>
#include <list>
#include <vector>
#include <sstream>
#include <string>
#include <cstdlib>   // For srand() and rand()
#include <ctime>     // For time()
#include <algorithm> // Ensure this is included
#include <iostream>
#include "DataStore.hpp"
#include <cairo/cairo.h>
#include <random>
#include <iomanip>
#include <cfloat>

typedef std::pair<double, double> Point_2;

namespace RPGraph
{
    GraphLayout::GraphLayout(UGraph &graph, float width, float height, int polygonId)
        : graph(graph), width(width), height(height), polygonId(polygonId)
    {

        coordinates = (Coordinate *)malloc(graph.num_nodes() * sizeof(Coordinate));
    }

    GraphLayout::~GraphLayout()
    {
        free(coordinates);
    }
float pointToSegmentDistance(const std::pair<float, float> &point,
                              const std::pair<float, float> &segStart,
                              const std::pair<float, float> &segEnd)
{
    float px = segEnd.first - segStart.first;
    float py = segEnd.second - segStart.second;
    float norm = px * px + py * py;
    float u = ((point.first - segStart.first) * px + (point.second - segStart.second) * py) / norm;

    u = std::clamp(u, 0.0f, 1.0f);
    float x = segStart.first + u * px;
    float y = segStart.second + u * py;

    return std::hypot(x - point.first, y - point.second);
}

std::pair<std::pair<double, double>, std::pair<double, double>> 
findClosestEdge(const std::vector<std::pair<double, double>> &polygonPoints, 
                const std::pair<double, double> &targetCenter)
{
    double minDistance = std::numeric_limits<double>::max();
    std::pair<std::pair<double, double>, std::pair<double, double>> closestEdge;

    for (size_t i = 0; i < polygonPoints.size(); ++i)
    {
        const auto &start = polygonPoints[i];
        const auto &end = polygonPoints[(i + 1) % polygonPoints.size()];

        double distance = pointToSegmentDistance(targetCenter, start, end);
        if (distance < minDistance)
        {
            minDistance = distance;
            closestEdge = {start, end};
        }
    }
    return closestEdge;
}

void GraphLayout::randomizePositions()
{
    DataStore &dataStore = DataStore::getInstance();
    const auto &polygonCenters = dataStore.polygonCenters;
    const auto &nodePolygons = dataStore.polygons;

    std::random_device rd;
    std::mt19937 gen(rd()); // Random number generator

    for (nid_t i = 0; i < graph.num_nodes(); ++i)
    {
        auto nodeId = graph.node_map_r[i];
        const auto &polygonPoints = nodePolygons.at(polygonId);

        if (polygonPoints.empty())
        {
            std::cerr << "Empty polygon for node " << nodeId << std::endl;
            continue;
        }

        float newX = 0, newY = 0;
        bool positioned = false;

        // Step 1: Find the nearest external polygon center
        double minDistance = std::numeric_limits<double>::max();
        std::pair<double, double> targetCenter;

        for (size_t j = 0; j < dataStore.externa_polygonX[nodeId].size(); ++j)
        {
            if (dataStore.externa_polygonX[nodeId][j] > 0) // External link exists
            {
                const auto &[centerXe, centerYe] = polygonCenters.at(j);
                double distance = std::hypot(centerXe - polygonCenters.at(polygonId).first,
                                             centerYe - polygonCenters.at(polygonId).second);

                if (distance < minDistance)
                {
                    minDistance = distance;
                    targetCenter = {centerXe, centerYe};
                }
            }
        }

        // Step 2: If external link exists, start around the external center
        if (minDistance < std::numeric_limits<double>::max())
        {
            auto [externalX, externalY] = targetCenter;
            auto [polygonX, polygonY] = polygonCenters.at(polygonId);

            // Generate a random starting position near the external polygon center
            std::uniform_real_distribution<> dis(-20.0, 20.0);
            newX = externalX + dis(gen);
            newY = externalY + dis(gen);

            // Step toward the polygon center until the position is inside
            while (!isPointInPolygon(newX, newY, polygonPoints))
            {
                // Move step by step toward the polygon center
                newX += (polygonX - newX) * 0.1f; // Move 10% closer to the polygon center
                newY += (polygonY - newY) * 0.1f;

                // Break if movement is too small to avoid infinite loops
                if (std::abs(polygonX - newX) == 0 && std::abs(polygonY - newY) ==0)
                    break;
            }

            positioned = true;
        }

        // Step 3: If no external position, randomize within the polygon
        if (!positioned)
        {
            auto [minX, minY, maxX, maxY] = getBoundingBox(polygonPoints);
            std::uniform_real_distribution<> disX(minX, maxX);
            std::uniform_real_distribution<> disY(minY, maxY);
            newX = disX(gen);
            newY = disY(gen);
        }

        // Step 4: Move the node to the calculated position
        move_node(i, newX, newY);
    }
}






    // Function to get the bounding box of a polygon
    std::tuple<float, float, float, float> GraphLayout::getBoundingBox(const std::vector<std::pair<double, double>> &polygonPoints)
    {
        float minX = std::numeric_limits<float>::max();
        float minY = std::numeric_limits<float>::max();
        float maxX = std::numeric_limits<float>::lowest();
        float maxY = std::numeric_limits<float>::lowest();

        for (const auto &point : polygonPoints)
        {
            if (point.first < minX)
                minX = point.first;
            if (point.second < minY)
                minY = point.second;
            if (point.first > maxX)
                maxX = point.first;
            if (point.second > maxY)
                maxY = point.second;
        }

        return {minX, minY, maxX, maxY};
    }

    float GraphLayout::getX(nid_t node_id)
    {
        return coordinates[node_id].x;
    }

    float GraphLayout::getY(nid_t node_id)
    {
        return coordinates[node_id].y;
    }

    float GraphLayout::minX()
    {
        float minX = std::numeric_limits<float>::max();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            float x = getX(n);
            if (x < minX)
                minX = x;
        }
        return minX;
    }

    float GraphLayout::maxX()
    {
        float maxX = std::numeric_limits<float>::min();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            float x = getX(n);
            if (x > maxX)
                maxX = x;
        }
        return maxX;
    }

    float GraphLayout::minY()
    {
        float minY = std::numeric_limits<float>::max();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            float y = getY(n);
            if (y < minY)
                minY = y;
        }
        return minY;
    }

    float GraphLayout::maxY()
    {
        float maxY = std::numeric_limits<float>::min();
        for (nid_t n = 0; n < graph.num_nodes(); ++n)
        {
            float y = getY(n);
            if (y > maxY)
                maxY = y;
        }
        return maxY;
    }

    float GraphLayout::getXRange()
    {
        return maxX() - minX();
    }

    float GraphLayout::getYRange()
    {
        return maxY() - minY();
    }

    float GraphLayout::getSpan()
    {
        return ceil(fmaxf(getXRange(), getYRange()));
    }

    float GraphLayout::getDistance(nid_t n1, nid_t n2)
    {
        const float dx = getX(n1) - getX(n2);
        const float dy = getY(n1) - getY(n2);
        return std::sqrt(dx * dx + dy * dy);
    }

    Real2DVector GraphLayout::getDistanceVector(nid_t n1, nid_t n2)
    {
        return Real2DVector(getX(n2) - getX(n1), getY(n2) - getY(n1));
    }

    Real2DVector GraphLayout::getNormalizedDistanceVector(nid_t n1, nid_t n2)
    {
        const float x1 = getX(n1);
        const float x2 = getX(n2);
        const float y1 = getY(n1);
        const float y2 = getY(n2);
        const float dx = x2 - x1;
        const float dy = y2 - y1;
        const float len = std::sqrt(dx * dx + dy * dy);

        return Real2DVector(dx / len, dy / len);
    }

    Coordinate GraphLayout::getCoordinate(nid_t node_id)
    {
        return coordinates[node_id];
    }

    Coordinate GraphLayout::getCenter()
    {
        float x = minX() + getXRange() / 2.0;
        float y = minY() + getYRange() / 2.0;
        return Coordinate(x, y);
    }

    void GraphLayout::setX(nid_t node_id, float x_value)
    {
        coordinates[node_id].x = x_value;
    }

    void GraphLayout::setY(nid_t node_id, float y_value)
    {
        coordinates[node_id].y = y_value;
    }

    bool GraphLayout::moveNode(nid_t n, RPGraph::Real2DVector v)
    {
        DataStore &dataStore = DataStore::getInstance();
        const auto &nodePolygons = dataStore.nodeToPolygon;
        const auto &polygonPoints = nodePolygons.at(graph.node_map_r[n]);
        // Calculate potential new position
        float newX = getX(n) + v.x;
        float newY = getY(n) + v.y;

        // Check if the new position is within the polygon
        if (isPointInPolygon(newX, newY, dataStore.polygons[polygonPoints]))
        {
            setX(n, newX);
            setY(n, newY);
            return true;
        }
        else
        {
            return false;
        }
    }

    float GraphLayout::getCx()
    {
        DataStore &dataStore = DataStore::getInstance();
        const auto &polygonCenters = dataStore.polygonCenters[polygonId];
        return polygonCenters.first;
    }
    float GraphLayout::getCy()
    {
        DataStore &dataStore = DataStore::getInstance();

        const auto &polygonCenters = dataStore.polygonCenters[polygonId];
        return polygonCenters.second;
    }

    float *GraphLayout::getPolygonPoints()
    {

        DataStore &dataStore = DataStore::getInstance();
        const auto &PolygonsUnique = dataStore.polygons;

        // Check if polygonId exists in PolygonsUnique
        if (PolygonsUnique.find(polygonId) != PolygonsUnique.end())
        {
            const auto &polygon = PolygonsUnique.at(polygonId);
            int size = polygon.size() * 2;
            float *points = new float[size];

            int index = 0;
            for (const auto &point : polygon)
            {

                points[index] = point.first; // Access x coordinate
                index++;
                points[index] = point.second; // Access y coordinate
                index++;
            }
            return points;
        }

        std::cout << "no points" << std::endl;
        return nullptr; // Return nullptr if polygonId is not found
    }

    // Helper function to calculate the distance between two points
    double distance(double x1, double y1, double x2, double y2)
    {
        return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
    }

    // Helper function to find the closest point on a line segment to a point
    std::pair<double, double> closestPointOnLineSegment(double px, double py, double ax, double ay, double bx, double by)
    {
        double dx = bx - ax;
        double dy = by - ay;
        if (dx == 0 && dy == 0)
        {
            return {ax, ay}; // The segment is a point
        }

        double t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy);
        t = std::max(0.0, std::min(1.0, t));
        return {ax + t * dx, ay + t * dy};
    }

    // Function to find the closest point in the polygon
    std::pair<double, double> closestPointInPolygon(double x, double y, const std::vector<std::pair<double, double>> &polygon, double margin, int polygonId)
    {
        DataStore &dataStore = DataStore::getInstance();
        const auto &polygonCenters = dataStore.polygonCenters[polygonId];

        double minDist = std::numeric_limits<double>::max();
        std::pair<float, float> closestPoint = {x, y};

        for (size_t i = 0; i < polygon.size(); ++i)
        {
            size_t j = (i + 1) % polygon.size();
            auto [closestX, closestY] = closestPointOnLineSegment(x, y, polygon[i].first, polygon[i].second, polygon[j].first, polygon[j].second);
            double dist = distance(x, y, closestX, closestY);
            if (dist < minDist)
            {
                minDist = dist;
                closestPoint = {closestX, closestY};
            }
        }

        // Move the closest point slightly towards the centroid
        auto [centroidX, centroidY] = polygonCenters;
        double adjustedX = closestPoint.first + margin * (centroidX - closestPoint.first) / distance(closestPoint.first, closestPoint.second, centroidX, centroidY);
        double adjustedY = closestPoint.second + margin * (centroidY - closestPoint.second) / distance(closestPoint.first, closestPoint.second, centroidX, centroidY);
        return {adjustedX, adjustedY};
    }

    bool GraphLayout::move_node(nid_t n, float x, float y)
    {
       DataStore &dataStore = DataStore::getInstance();
    const auto &polygonPoints = dataStore.polygons[polygonId];

    // Check if the new position is within the polygon
    if (isPointInPolygon(x, y, polygonPoints))
    {
        setX(n, x);
        setY(n, y);
        return true;
    }

    // Try finding the closest point within the polygon
    float margin = sqrt(dataStore.polygonAreas[polygonId]) / graph.num_nodes();
    std::pair<float, float> closestPoint = closestPointInPolygon(x, y, polygonPoints, margin, polygonId);

    if (isPointInPolygon(closestPoint.first, closestPoint.second, polygonPoints))
    {
        setX(n, closestPoint.first);
        setY(n, closestPoint.second);
        return true;
    }

    // If still outside, fallback to the polygon's center
    const auto &polygonCenter = dataStore.polygonCenters[polygonId];
    setX(n, polygonCenter.first);
    setY(n, polygonCenter.second);
    return isPointInPolygon(polygonCenter.first, polygonCenter.second, polygonPoints);
    }

bool GraphLayout::isPointInPolygon(float x, float y, const std::vector<Point_2> &polygon)
{
     DataStore &dataStore = DataStore::getInstance();
    const auto &polygonPoints = dataStore.polygons[polygonId];
    // Check if point is in polygon with a tolerance for boundary points
    bool inside = false;
    size_t j = polygonPoints.size() - 1;
    for (size_t i = 0; i < polygonPoints.size(); i++)
    {
        if (((polygonPoints[i].second > y) != (polygonPoints[j].second > y)) &&
            (x < (polygonPoints[j].first - polygonPoints[i].first) * (y - polygonPoints[i].second) /
                      (polygonPoints[j].second - polygonPoints[i].second) +
                  polygonPoints[i].first))
        {
            inside = !inside;
        }
        j = i;
    }

    // Allow points on the boundary
    for (size_t i = 0; i < polygonPoints.size(); i++)
    {
        const auto &p1 = polygonPoints[i];
        const auto &p2 = polygonPoints[(i + 1) % polygonPoints.size()];
        float dx = p2.first - p1.first;
        float dy = p2.second - p1.second;
        float lengthSquared = dx * dx + dy * dy;

        // Projection calculation
        float t = ((x - p1.first) * dx + (y - p1.second) * dy) / lengthSquared;
        t = std::max(0.0f, std::min(1.0f, t));

        float projX = p1.first + t * dx;
        float projY = p1.second + t * dy;
float tolerance = 1e-5;
        // Check if point is close to the edge
        if (std::hypot(projX - x, projY - y) <= tolerance)
        {
            return true; // On the edge
        }
    }

    return inside;
}


    void GraphLayout::setCoordinates(nid_t node_id, Coordinate c)
    {
        setX(node_id, c.x);
        setY(node_id, c.y);
    }

    inline float distance(float x1, float y1, float x2, float y2)
    {
        return std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    }

    // Function to apply repulsion between nodes
    void apply_repulsion(std::vector<float> &x_pos, std::vector<float> &y_pos, int num_nodes, float min_distance, float repulsion_strength)
    {
        for (int i = 0; i < num_nodes; ++i)
        {
            for (int j = i + 1; j < num_nodes; ++j)
            {
                float dx = x_pos[j] - x_pos[i];
                float dy = y_pos[j] - y_pos[i];
                float dist = distance(x_pos[i], y_pos[i], x_pos[j], y_pos[j]);

                if (dist < min_distance && dist > 0)
                {
                    float force = repulsion_strength / dist;

                    // Apply repulsion
                    x_pos[i] -= force * dx;
                    y_pos[i] -= force * dy;
                    x_pos[j] += force * dx;
                    y_pos[j] += force * dy;
                }
            }
        }
    }

    void GraphLayout::writeToSVG(std::string path, int real_number_nodes, std::string path_node_map, std::string path_dep_map)
    {
        std::cerr << "************************: " << real_number_nodes << std::endl;
        std::ifstream file_node_map(path_node_map);

        if (!file_node_map.is_open())
        {
            std::cerr << "Error opening file: " << path_node_map << std::endl;
        }

        std::vector<int> node_map_r_f;
        int index_node_map, value_node_map;

        // Read each line of the file
        std::string line_node_map;
        while (std::getline(file_node_map, line_node_map))
        {
            std::istringstream iss(line_node_map);
            if (iss >> value_node_map >> index_node_map)
            {
                // Ensure the vector is large enough to hold the value at 'index'
                if (index_node_map >= node_map_r_f.size())
                {
                    node_map_r_f.resize(index_node_map + 1);
                }
                node_map_r_f[index_node_map] = value_node_map;
            }
        }

        file_node_map.close();

        std::ifstream file_dep_map(path_dep_map);

        if (!file_dep_map.is_open())
        {
            std::cerr << "Error opening file: " << path_dep_map << std::endl;
        }

        std::vector<int> dep_map_r;
        int index_dep_map, value_dep_map;

        // Read each line of the file
        std::string line_dep_map;
        while (std::getline(file_dep_map, line_dep_map))
        {
            std::istringstream iss(line_dep_map);
            if (iss >> value_dep_map >> index_dep_map)
            {
                // Ensure the vector is large enough to hold the value at 'index'
                if (index_dep_map >= dep_map_r.size())
                {
                    dep_map_r.resize(index_dep_map + 1);
                }
                dep_map_r[index_dep_map] = value_dep_map;
            }
        }

        file_dep_map.close();

        DataStore &dataStore = DataStore::getInstance();
        const auto &nodePolygons = dataStore.nodeToPolygon;
        const auto &polygons = dataStore.polygons;
        const float node_opacity = 10000000.0 / graph.num_nodes();
        const float edge_opacity = 100.0 / graph.num_edges();

        // Open the SVG file
        std::ofstream svg_file(path);
        if (!svg_file.is_open())
        {
            std::cerr << "Error: Unable to open SVG file at " << path << std::endl;
            return;
        }
        float minX = FLT_MAX;
        float minY = FLT_MAX;
        float maxX = FLT_MIN;
        float maxY = FLT_MIN;

        // Update minX, minY, maxX, maxY based on node positions
        for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
        {
            double x = getX(n1);
            double y = getY(n1);
            if (x < minX)
                minX = x;
            if (y < minY)
                minY = y;
            if (x > maxX)
                maxX = x;
            if (y > maxY)
                maxY = y;
        }

        // Similarly, update minX, minY, maxX, maxY based on polygon points
        for (const auto &[polygonId, polygonPoints] : dataStore.polygons)
        {
            for (const auto &point : polygonPoints)
            {
                double x = point.first;
                double y = point.second;
                if (x < minX)
                    minX = x;
                if (y < minY)
                    minY = y;
                if (x > maxX)
                    maxX = x;
                if (y > maxY)
                    maxY = y;
            }
        }

        float width = maxX - minX;
        float height = maxY - minY;

        svg_file << "<svg width=\"" << width * 0.17912 << "pt\" height=\"" << height * 0.15221 << "pt\" viewBox=\""
                 << minX << " " << minY << " " << width << " " << height
                 << "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
                 << std::endl;
        svg_file << "<g id=\"graph0\" class=\"graph\" transform=\"scale(1 1) rotate(0) \" >" << std::endl;
        for (const auto &[polygonId, polygonpoints] : dataStore.polygons)
        {
            try
            {
                const auto &polygonPoints = polygonpoints; // Correct way to get polygon points
                std::string points_str;
                for (const auto &point : polygonPoints)
                {
                    if (!points_str.empty())
                        points_str += " ";
                    points_str += std::to_string(point.first) + "," + std::to_string(point.second);
                }
                svg_file << "<title>" << polygonId << "</title>" << std::endl;
                svg_file << "<polygon fill=\"#dae2ff\" stroke=\"#000000\" points=\"" << points_str << "\"/>" << std::endl;
                svg_file << "<polyline fill=\"none\" stroke=\"#000000\" points=\"" << points_str << "\"/>" << std::endl;
                svg_file << "<text x=\"" << dataStore.polygonCenters[polygonId].first << "\" y=\"" << dataStore.polygonCenters[polygonId].second << "\" ext-anchor=\"middle\" fill=\"black\" font-size=\"" << 3 << "\" opacity=\"0.3\" >" << dep_map_r[polygonId] << "</text>" << std::endl;
            }
            catch (const std::out_of_range &e)
            {
                std::cout << "Error: Polygon ID " << polygonId << " not found." << std::endl;
            }
        }
        // Draw edges
        /*for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1) {
            for (nid_t n2 : graph.neighbors_with_geq_id(n1)) {

                svg_file << "<line x1=\"" << getX(n1)
                         << "\" y1=\"" << getY(n1)
                         << "\" x2=\"" << getX(n2)
                         << "\" y2=\"" << getY(n2)
                         << "\" style=\"stroke:rgb(0,0,0);stroke-opacity:.5" << 1
                         << ";stroke-width:.5\"/>" << std::endl;
            }
        }*/

        for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
        {
            for (nid_t n2 : graph.neighbors_with_geq_id(n1))
            {
                if (graph.node_map_r[n1] <= real_number_nodes && graph.node_map_r[n2] <= real_number_nodes)
                {

                    // Get coordinates for nodes n1 and n2
                    double x1 = getX(n1);
                    double y1 = getY(n1);
                    double x2 = getX(n2);
                    double y2 = getY(n2);

                    // Calculate control points to create a curved path between n1 and n2
                    // These are simple approximations; you can adjust them for a more specific curve
                    double cx1 = (x1 + x2) / 2; // First control point x-coordinate (midpoint for simplicity)
                    double cy1 = y1 - 50;       // First control point y-coordinate (above the line)
                    double cx2 = (x1 + x2) / 2; // Second control point x-coordinate (midpoint)
                    double cy2 = y2 + 50;       // Second control point y-coordinate (below the line)

                    // Write the <path> element to the SVG file
                    svg_file << "<!-- " << node_map_r_f[graph.node_map_r[n1]] << "&#45;&#45;" << node_map_r_f[graph.node_map_r[n2]] << " -->" << std::endl;
                    svg_file << "<g id=\"edge1\" class=\"edge\">" << std::endl;
                    svg_file << "<title>" << node_map_r_f[graph.node_map_r[n1]] << "&#45;&#45;" << node_map_r_f[graph.node_map_r[n2]] << "</title>" << std::endl;

                    /*svg_file << "<path fill=\"none\" stroke=\"#c0c0c0\" d=\""
                            << "M" << x1 << "," << y1 << " "      // Move to the start point (n1)
                            << "C" << cx1 << "," << cy1 << " "    // First control point
                            << cx2 << "," << cy2 << " "           // Second control point
                            << x2 << "," << y2 << "\" />"         // End at the second node (n2)
                            << std::endl;
                    svg_file<<"</g>"<< std::endl;*/
                    if (dataStore.nodeToPolygon[node_map_r_f[graph.node_map_r[n2]]] != dataStore.nodeToPolygon[node_map_r_f[graph.node_map_r[n1]]])
                    {
                       
                        svg_file << "<path fill=\"none\" stroke=\"#c0c0c0\" opacity=\"0.4\" stroke-width=\"2\" d=\""

                                 << "M" << x1 << "," << y1 << " "     // Move to the start point (n1)
                                 << "L" << x2 << "," << y2 << "\" />" // Draw a straight line to the end point (n2)
                                 << std::endl;
                    }
                    svg_file << "</g>" << std::endl;
                }
            }
        }

     for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
        {
            for (nid_t n2 : graph.neighbors_with_geq_id(n1))
            {
                if (graph.node_map_r[n1] <= real_number_nodes && graph.node_map_r[n2] <= real_number_nodes)
                {

                    // Get coordinates for nodes n1 and n2
                    double x1 = getX(n1);
                    double y1 = getY(n1);
                    double x2 = getX(n2);
                    double y2 = getY(n2);

                    // Calculate control points to create a curved path between n1 and n2
                    // These are simple approximations; you can adjust them for a more specific curve
                    double cx1 = (x1 + x2) / 2; // First control point x-coordinate (midpoint for simplicity)
                    double cy1 = y1 - 50;       // First control point y-coordinate (above the line)
                    double cx2 = (x1 + x2) / 2; // Second control point x-coordinate (midpoint)
                    double cy2 = y2 + 50;       // Second control point y-coordinate (below the line)

                    // Write the <path> element to the SVG file
                    svg_file << "<!-- " << node_map_r_f[graph.node_map_r[n1]] << "&#45;&#45;" << node_map_r_f[graph.node_map_r[n2]] << " -->" << std::endl;
                    svg_file << "<g id=\"edge1\" class=\"edge\">" << std::endl;
                    svg_file << "<title>" << node_map_r_f[graph.node_map_r[n1]] << "&#45;&#45;" << node_map_r_f[graph.node_map_r[n2]] << "</title>" << std::endl;

                    /*svg_file << "<path fill=\"none\" stroke=\"#c0c0c0\" d=\""
                            << "M" << x1 << "," << y1 << " "      // Move to the start point (n1)
                            << "C" << cx1 << "," << cy1 << " "    // First control point
                            << cx2 << "," << cy2 << " "           // Second control point
                            << x2 << "," << y2 << "\" />"         // End at the second node (n2)
                            << std::endl;
                    svg_file<<"</g>"<< std::endl;*/
                    if (dataStore.nodeToPolygon[node_map_r_f[graph.node_map_r[n2]]] == dataStore.nodeToPolygon[node_map_r_f[graph.node_map_r[n1]]])
                    {
                        svg_file << "<path fill=\"none\" stroke=\"#000000\" stroke-width=\"5\" d=\""

                                 << "M" << x1 << "," << y1 << " "     // Move to the start point (n1)
                                 << "L" << x2 << "," << y2 << "\" />" // Draw a straight line to the end point (n2)
                                 << std::endl;
                    }
                 
                    svg_file << "</g>" << std::endl;
                }
            }
        }
        // Draw nodes
        /*for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1) {
            svg_file << "<circle cx=\"" << getX(n1)
                     << "\" cy=\"" << getY(n1)
                     << "\" r=\"" << 5.0 << "\" style=\"fill:rgb(255,0,0);fill-opacity:" << 100 << "\">"
                     << "<title>Node ID: " << graph.node_map_r[n1] <<" for polygon" <<dataStore.nodeToPolygon[graph.node_map_r[n1]]<< "</title>"
                     << "</circle>" << std::endl;
        }*/
        for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
        {
            if (graph.node_map_r[n1] <= real_number_nodes)
            {

                svg_file << "<!-- " << node_map_r_f[graph.node_map_r[n1]] << " -->" << std::endl;
                svg_file << "<g id=\"node" << node_map_r_f[graph.node_map_r[n1]] << "\" class=\"node\">" << std::endl;
                svg_file << "<title>" << node_map_r_f[graph.node_map_r[n1]] << "</title>" << std::endl;
                svg_file << "<ellipse fill=\"transparent\" fill-opacity=\"0.8\" stroke=\"#000000\" cx=\"" << getX(n1) << "\" cy=\"" << getY(n1) << "\" rx=\"10\" ry=\"10\"/>" << std::endl;
                svg_file << "<text text-anchor=\"middle\" x=\"" << getX(n1) << "\" y=\"" << getY(n1) << "\" font-family=\"Helvetica,sans-Serif\" font-weight=\"bold\" font-size=\"" << 14 << "\"  fill=\"#000000\">" << node_map_r_f[graph.node_map_r[n1]] << "</text>" << std::endl;
                svg_file << "</g>" << std::endl;
            }
        }

        // Draw polygons

        // Write SVG footer and close the file
        svg_file << "</g>" << std::endl;
        svg_file << "</svg>" << std::endl;
        svg_file.close();

        std::cout << "SVG file successfully written to: " << path << std::endl;

        // Now create the dot file
        std::string dot_path = path + ".dot";
        std::ofstream dot_file(dot_path);
        if (!dot_file.is_open())
        {
            std::cerr << "Error: Unable to open dot file at " << dot_path << std::endl;
            return;
        }

        // Write the header for the dot file
        dot_file << "digraph G {" << std::endl;
        dot_file << "    graph [layout=neato, overlap=false];" << std::endl;

        // Iterate over all nodes and write their attributes in the dot format
        for (nid_t i = 0; i < graph.num_nodes(); ++i)
        {
            auto nodeId = graph.node_map_r[i];

            int polygonId = dep_map_r[dataStore.nodeToPolygon[nodeId]]; // Get the department ID (polygon ID)
            Point_2 center = dataStore.polygonCenters[polygonId];
            // Get node position using the GraphLayout methods
            float posX = getX(i);
            float posY = getY(i);

            // Write node attributes
            dot_file << "    " << node_map_r_f[nodeId] << " ["
                     << "cluster=" << 100 << ", \n"
                     << "department=" << 100 << ", \n"
                     << "height=0.5, \n"
                     << "name=\"" << node_map_r_f[nodeId] << "\", \n"
                     << "pos=\"" << posX << "," << posY << "\", \n"
                     << "weight=0, \n"
                     << "width=0.83788 \n"
                     << "];" << std::endl;
        }

        // Close the dot file
        dot_file << "}" << std::endl;

        dot_file.close();
        std::cout << "Dot file successfully written to: " << dot_path << std::endl;

        std::string node_positions_path = path + ".txt";
        std::ofstream node_position_file(node_positions_path);
        if (!node_position_file.is_open())
        {
            std::cerr << "Error: Unable to open file at " << node_positions_path << std::endl;
            return;
        }

        // Collect node positions
        std::vector<float> x_pos(graph.num_nodes());
        std::vector<float> y_pos(graph.num_nodes());

        for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
        {
            if (graph.node_map_r[n1] <= real_number_nodes)
            {
                x_pos[n1] = getX(n1);
                y_pos[n1] = getY(n1);
            }
        }

        // Apply repulsion to avoid overlaps
        float min_distance = 100.0f;     // Minimum allowable distance between nodes
        float repulsion_strength = 1.0f; // Strength of repulsion
        // for(int i=0;i<100;i++)
        // apply_repulsion(x_pos, y_pos, graph.num_nodes(), min_distance, repulsion_strength);

        // Write final node positions to the file
        for (nid_t n1 = 0; n1 < graph.num_nodes(); ++n1)
        {
            if (graph.node_map_r[n1] <= real_number_nodes)
            {
                node_position_file << graph.node_map_r[n1] << " " << x_pos[n1] << "," << y_pos[n1] << std::endl;
            }
        }

        node_position_file.close();
        std::cout << "Node positions successfully written to: " << node_positions_path << std::endl;
    }

}
