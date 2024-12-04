

#include "RPForceAtlas2.hpp"
#include "RPCommon.hpp" // Include for singleton access

namespace RPGraph
{
    ForceAtlas2::ForceAtlas2(GraphLayout &layout, bool use_barneshut,
                             bool strong_gravity, float gravity, float scale,int max_iterations)
    : LayoutAlgorithm(layout), use_barneshut{use_barneshut},
      strong_gravity{strong_gravity}
    {
        iteration = 0;

        setGravity(gravity);
        setScale(scale);

        global_speed = 1.0;
        speed_efficiency = 1.0;
        jitter_tolerance = 1.0;

        k_s = 0.1;
        k_s_max = 10.0;

        theta = 1.0;
        epssq  = 0.05 * 0.05;
        itolsq = 1.0f / (theta * theta);

        delta = 0.0;

        prevent_overlap = false;
        use_linlog = false;

        layout.randomizePositions();

    
    }

    ForceAtlas2::~ForceAtlas2(){};

    void ForceAtlas2::doSteps(int n)
    {
        for (int i = 0; i < n; ++i) doStep(i);
    }

    void ForceAtlas2::setScale(float s)
    {
        k_r = s;
    }

    void ForceAtlas2::setGravity(float g)
    {
        k_g = g;
    }

    float ForceAtlas2::mass(nid_t n)
    {
        return layout.graph.degree(n) + 1.0;
    }
}
