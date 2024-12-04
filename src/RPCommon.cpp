
#include "RPCommon.hpp"
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include <algorithm> // Include this for std::find_if

bool is_file_exists(std::string filepath)
{
    std::ifstream infile(filepath);
    return infile.good();
}

std::string basename(std::string filepath)
{
    char *result_p = new char[filepath.size() + 1];
    strcpy(result_p, filepath.c_str());
    std::string result = basename(result_p);
    delete[] result_p;
    return result;
}

namespace RPGraph
{
    float get_random(float lowerbound, float upperbound)
    {
        return lowerbound + (upperbound - lowerbound) * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    }

    Real2DVector::Real2DVector(float x, float y): x(x), y(y) {};

    float Real2DVector::magnitude()
    {
        return std::sqrt(x*x + y*y);
    }

    float Real2DVector::distance(RPGraph::Real2DVector to)
    {
        const float dx = (x - to.x)*(x - to.x);
        const float dy = (y - to.y)*(y - to.y);
        return std::sqrt(dx*dx + dy*dy);
    }

    Real2DVector Real2DVector::operator*(float b)
    {
        return Real2DVector(this->x * b, this->y * b);
    }

    Real2DVector Real2DVector::operator/(float b)
    {
        return Real2DVector(this->x / b, this->y / b);
    }

    Real2DVector Real2DVector::operator+(Real2DVector b)
    {
        return Real2DVector(this->x + b.x, this->y + b.y);
    }

    Real2DVector Real2DVector::operator-(Real2DVector b)
    {
        return Real2DVector(this->x - b.x, this->y - b.y);
    }

    void Real2DVector::operator+=(Real2DVector b)
    {
        this->x += b.x;
        this->y += b.y;
    }

    Real2DVector Real2DVector::getNormalized()
    {
        const float length = magnitude();
        return Real2DVector(this->x / length, this->y / length);
    }

    Real2DVector Real2DVector::normalize()
    {
        const float length = magnitude();
        this->x /= length;
        this->y /= length;
        return *this;
    }

    Coordinate::Coordinate(float x, float y): x(x), y(y) {};

    Coordinate Coordinate::operator+(float b)
    {
        return Coordinate(this->x + b, this->y + b);
    }

    Coordinate Coordinate::operator*(float b)
    {
        return Coordinate(this->x * b, this->y * b);
    }

    Coordinate Coordinate::operator/(float b)
    {
        return Coordinate(this->x / b, this->y / b);
    }

    Coordinate Coordinate::operator+(Real2DVector b)
    {
        return Coordinate(this->x + b.x, this->y + b.y);
    }

    Coordinate Coordinate::operator-(Coordinate b)
    {
        return Coordinate(this->x - b.x, this->y - b.y);
    }

    bool Coordinate::operator==(Coordinate b)
    {
        return (this->x == b.x && this->y == b.y);
    }

    float Coordinate::distance(RPGraph::Coordinate to)
    {
        return std::sqrt((x - to.x)*(x - to.x) + (y - to.y)*(y - to.y));
    }

    float Coordinate::distance2(RPGraph::Coordinate to)
    {
        return (x - to.x)*(x - to.x) + (y - to.y)*(y - to.y);
    }

    void Coordinate::operator/=(float b)
    {
        this->x /= b;
        this->y /= b;
    }

    void Coordinate::operator+=(RPGraph::Coordinate b)
    {
        this->x += b.x;
        this->y += b.y;
    }

    void Coordinate::operator+=(RPGraph::Real2DVector b)
    {
        this->x += b.x;
        this->y += b.y;
    }

    int Coordinate::quadrant()
    {
        if (x <= 0)
        {
            if (y >= 0) return 0;
            else        return 3;

        }
        else
        {
            if (y >= 0) return 1;
            else        return 2;
        }
    }

    float distance(Coordinate from, Coordinate to)
    {
        const float dx = from.x - to.x;
        const float dy = from.y - to.y;
        return std::sqrt(dx*dx + dy*dy);
    }

    float distance2(Coordinate from, Coordinate to)
    {
        const float dx = from.x - to.x;
        const float dy = from.y - to.y;
        return dx*dx + dy*dy;
    }

    Real2DVector normalizedDirection(Coordinate from, Coordinate to)
    {
        const float dx = from.x - to.x;
        const float dy = from.y - to.y;
        const float len = std::sqrt(dx*dx + dy*dy);
        return Real2DVector(dx/len, dy/len);
    }

    Real2DVector direction(Coordinate from, Coordinate to)
    {
        const float dx = from.x - to.x;
        const float dy = from.y - to.y;
        return Real2DVector(dx, dy);
    }
}









const std::vector<int>& RPCommon::getNodes() const {
    return nodes;
}

const std::vector<std::pair<int, int>>& RPCommon::getEdges() const {
    return edges;
}

// Custom hash function for vector of pairs
struct VectorPairHash {
    size_t operator()(const std::vector<std::pair<double, double>>& v) const {
        std::hash<double> hash_fn;
        size_t hash = 0;
        for (const auto& p : v) {
            hash ^= hash_fn(p.first) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            hash ^= hash_fn(p.second) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// Custom equality function for vector of pairs
struct VectorPairEqual {
    bool operator()(const std::vector<std::pair<double, double>>& v1, const std::vector<std::pair<double, double>>& v2) const {
        return v1 == v2;
    }
};

