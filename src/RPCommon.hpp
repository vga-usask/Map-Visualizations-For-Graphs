
#ifndef RPCommonUtils_hpp
#define RPCommonUtils_hpp
#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

#ifdef __NVCC__
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <stdlib.h>

#define cudaCatchError(ans) { assert_d((ans), __FILE__, __LINE__); }
inline void assert_d(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"error: (GPUassert) %s (error %d). %s:%d\n", cudaGetErrorString(code), code, file, line);
        if (abort) exit(code);
    }
}
#endif
bool is_file_exists(std::string filepath);
std::string basename(std::string filepath);

namespace RPGraph
{
    float get_random(float lowerbound, float upperbound);

    class Real2DVector
    {
    public:
        Real2DVector(float x, float y);
        float x, y;
        float magnitude();
        float distance(Real2DVector to);

        Real2DVector operator*(float b);
        Real2DVector operator/(float b);
        Real2DVector operator+(Real2DVector b);
        Real2DVector operator-(Real2DVector b);
        void operator+=(Real2DVector b);
        float dot(const Real2DVector& b);
        float length();
        Real2DVector getNormalized();
        Real2DVector normalize();
        float angleBetween(const Real2DVector& b);
    };

    class Coordinate
    {
    public:
        float x, y;
        Coordinate(float x, float y);

        Coordinate operator+(float b);
        Coordinate operator*(float b);
        Coordinate operator/(float b);
        Coordinate operator+(Real2DVector b);
        Coordinate operator-(Coordinate b);
        bool operator==(Coordinate b);
        void operator/=(float b);
        void operator+=(Coordinate b);
        void operator+=(RPGraph::Real2DVector b);

        int quadrant();
        float distance(Coordinate to);
        float distance2(Coordinate to);
    };

    float distance(Coordinate from, Coordinate to);
    float distance2(Coordinate from, Coordinate to);
    Real2DVector normalizedDirection(Coordinate from, Coordinate to);
    Real2DVector direction(Coordinate from, Coordinate to);
}

// New additions for reading updated_nodes_with_polygons.txt
class RPCommon {
public:
    RPCommon(const RPCommon&) = delete;
    void operator=(const RPCommon&) = delete;
    static RPCommon& getInstance();
    const std::unordered_map<int, std::vector<std::pair<double, double>>>& getNodePolygons() const;
    const std::vector<std::vector<std::pair<double, double>>>& getPolygonsUnique() const;
    const std::unordered_map<int, int>& getNodeToPolygon() const;

    const std::vector<int>& getNodes() const;
    const std::vector<std::pair<int, int>>& getEdges() const;

private:
    RPCommon();

    std::unordered_map<int, std::vector<std::pair<double, double>>> nodePolygons;
    std::unordered_map<int, int> nodeToPolygon;
    std::vector<std::vector<std::pair<double, double>>> polygonsUnique;
    std::vector<int> nodes;
    std::vector<std::pair<int, int>> edges;
};

bool is_file_exists(std::string filepath);
std::string basename(std::string filepath);
#endif /* RPCommonUtils_hpp */

