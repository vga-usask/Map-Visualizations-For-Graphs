
#ifndef RPForceAtlas2_hpp
#define RPForceAtlas2_hpp

#include "RPLayoutAlgorithm.hpp"
#include "RPBarnesHutApproximator.hpp"
#include "RPCommon.hpp"  // Include for singleton access

namespace RPGraph
{
    class ForceAtlas2 : public LayoutAlgorithm
    {
        public:
            ForceAtlas2(GraphLayout &layout, bool use_barneshut,
                        bool strong_gravity, float gravity, float scale,int max_iterations);
             ~ForceAtlas2();

            virtual void doStep(int iter) = 0;
        
            void doSteps(int n);
            void setScale(float s);
            void setGravity(float s);
            float mass(nid_t n);
            bool prevent_overlap, use_barneshut, use_linlog, strong_gravity,terminate;
           
        protected:
            int iteration;
            int max_iterations2;
            float k_r, k_g; // scalars for repulsive and gravitational force.
            float delta; // edgeweight influence.
            float global_speed;

            // Parameters used in adaptive temperature
            float speed_efficiency, jitter_tolerance;
            float k_s, k_s_max; // magic constants related to swinging.

            // Barnes-Hut parameters
            float theta;   // Accuracy
            float epssq;   // Softening (Epsilon, squared)
            float itolsq;  // Inverse tolerance, squared
    };
}
#endif
