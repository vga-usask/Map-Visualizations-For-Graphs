#include "DataStore.hpp"
#include <cmath> // For mathematical functions
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <algorithm> 
float cornergravity = 1.0f;
int polysize = 1;



double DataStore::calculatePolygonArea(const std::vector<Point_2>& points) {
    // Using the Shoelace formula for polygon area
    double area = 0.0;
    size_t n = points.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        area += points[i].first * points[j].second;
        area -= points[j].first * points[i].second;
    }
    return std::abs(area) / 2.0;
}

std::pair<double, double> DataStore::calculatePolygonCenter(const std::vector<Point_2>& points) {
    double centroidX = 0.0, centroidY = 0.0;
    for (const auto& point : points) {
        centroidX += point.first;
        centroidY += point.second;
    }
    centroidX /= points.size();
    centroidY /= points.size();
    return {centroidX, centroidY};
}




void DataStore::addEdge(int src, int dst) {
    edges.emplace_back(src, dst);
  
}

void DataStore::addNode(int node) {
    nodes.insert(node);
}



