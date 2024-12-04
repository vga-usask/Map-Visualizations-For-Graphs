
#ifndef RPBarnesHutApproximator_hpp
#define RPBarnesHutApproximator_hpp

#include "RPGraph.hpp"
#include "RPCommon.hpp"

namespace RPGraph
{
    class BarnesHutCell
    {
    public:
        void add_leafcell(int quadrant, float mass, Coordinate pos);
        float lb, rb, ub, bb;

        // BarnesHutCell always contain either a single particle, or subcells (at most 4).
        BarnesHutCell(Coordinate position, float length, Coordinate particle_position, float particle_mass);
        ~BarnesHutCell();

        Coordinate cell_center, mass_center;
        nid_t num_subparticles = 0;
        float total_mass;
        const float length;   // length of a cell = width = height
        BarnesHutCell *sub_cells[4] = {nullptr, nullptr, nullptr, nullptr}; // per quadrant.

        void insertParticle(Coordinate particle_position, float particle_mass);
    };

    class BarnesHutApproximator
    {
    public:
        BarnesHutApproximator(Coordinate root_center, float root_length, float theta);
        Real2DVector approximateForce(Coordinate particle_pos, float particle_mass, float theta);
        void insertParticle(Coordinate particle_position, float particle_mass);

        void reset(Coordinate root_center, float root_length);
        void setTheta(float theta);

    private:
        BarnesHutCell *root_cell = nullptr;
        const float theta;
        Coordinate root_center;
        float root_length;

    };
}

#endif /* RPBarnesHutApproximator_hpp */
