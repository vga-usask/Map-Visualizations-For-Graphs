
#ifndef RPGraphLayout_hpp
#define RPGraphLayout_hpp

#include "RPGraph.hpp"
#include "RPCommon.hpp"
#include <string>
#include <vector>

typedef std::pair<double, double> Point_2;
extern bool area_difference_flag;

namespace RPGraph
{
    class GraphLayout
    {
    private:
        Coordinate *coordinates;

    protected:
        float width, height;
        float minX(), minY(), maxX(), maxY();

    public:
        GraphLayout(RPGraph::UGraph &graph,
                    float width=10000, float height=10000,int polygonId=-1);
        ~GraphLayout();

        UGraph &graph; // to lay-out
		int polygonId;
        std::vector<std::vector<int>> externa_polygonX;
        // randomize the layout position of all nodes.
        void randomizePositions();
float* getPolygonPoints(); 
        float getX(nid_t node_id), getY(nid_t node_id);
        float getXRange(), getYRange(), getSpan();
        float getDistance(nid_t n1, nid_t n2);
        Real2DVector getDistanceVector(nid_t n1, nid_t n2);
        Real2DVector getNormalizedDistanceVector(nid_t n1, nid_t n2);
        Coordinate getCoordinate(nid_t node_id);
        Coordinate getCenter();
 std::pair<float, float> getCenterf(int n);
        void setX(nid_t node_id, float x_value), setY(nid_t node_id, float y_value);
        bool isPointInPolygon(float x, float y, const std::vector<std::pair<double, double>>& polygon);
        bool moveNode(nid_t, Real2DVector v);
        bool  move_node(nid_t n, float x, float y);
        float getCx();
        float getCy();
        void setCoordinates(nid_t node_id, Coordinate c);

        void writeToSVG(std::string path, int real_number_nodes,std::string node_map_r_f,std::string dep_map_r);
        std::vector<Point_2> generateBoundaryNodes(const std::vector<Point_2>& polygonPoints, int count);
        double calculatePolygonArea(const std::vector<Point_2>& points);
        std::tuple<float, float, float, float> getBoundingBox(const std::vector<std::pair<double, double>>& polygonPoints) ;
    };
}

#endif /* RPGraphLayout_hpp */

