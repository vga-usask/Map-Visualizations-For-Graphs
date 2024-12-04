#ifndef DATASTORE_HPP
#define DATASTORE_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath> // For mathematical functions
#include <boost/functional/hash.hpp> // Include Boost hash

typedef std::pair<double, double> Point_2;

class DataStore {
public:
    // Singleton instance method
    static DataStore& getInstance() {
        static DataStore instance;
        return instance;
    }

    // Disable copy constructor and assignment operator for singleton
    DataStore(const DataStore&) = delete;
    DataStore& operator=(const DataStore&) = delete;

    // Methods
    void addEdge(int src, int dst);
    void addNode(int node);
    double calculatePolygonArea(const std::vector<Point_2>& points);
    std::pair<double, double> calculatePolygonCenter(const std::vector<Point_2>& points);

    // Member variables
    std::unordered_map<int, int> nodeToPolygon;
    std::unordered_map<int, std::vector<int>> PolygonToNode;
    std::unordered_map<int, std::pair<double, double>> polygonCenters;
    std::unordered_map<int, double> polygonAreas;
    std::vector<std::pair<int, int>> edges;
    std::unordered_set<int> nodes;
    std::unordered_map<int, std::vector<std::pair<double, double>>> polygons;  // Map polygonId -> polygonPoints
    float maximum_area = 0.0f;
    int maximum_area_polygonId = -1;
    int number_of_polygons = 0;
    int number_of_points = 0;

    // Member initialization
    std::vector<std::vector<int>> externa_polygonX;

private:
    // Private constructor for singleton
    DataStore() : externa_polygonX(100000, std::vector<int>(500)) {}
};

#endif // DATASTORE_HPP
