/*
 ==============================================================================

 RPGPUForceAtlas2.hpp
 Copyright © 2016, 2017, 2018  G. Brinkmann

 This file is part of graph_viewer.

 graph_viewer is free software: you can redistribute it and/or modify
 it under the terms of version 3 of the GNU Affero General Public License as
 published by the Free Software Foundation.

 graph_viewer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with graph_viewer.  If not, see <https://www.gnu.org/licenses/>.

 ==============================================================================
*/

#ifndef RPGPUForceAtlas2_hpp
#define RPGPUForceAtlas2_hpp
#include "RPForceAtlas2.hpp"
#include "DataStore.hpp"
#include <vector>
namespace RPGraph
{
    class CUDAForceAtlas2: public ForceAtlas2
    {
    public:
        CUDAForceAtlas2(GraphLayout &layout, bool use_barneshut,
                        bool strong_gravity, float gravity, float scale, int max_iterations,std::vector<std::vector<int>> &externa_polygonX,float maximum_area);
        ~CUDAForceAtlas2();
        void doStep(int iter) override;
      
        bool sync_layout() override;
        void sync_layout_void() override;
void freeGPUMemory();

    private:
        float* d_points;
        int* d_externa_polygonX;
        float* d_externa_polygonX_x;
        float* d_externa_polygonX_y;
        float* d_effective_d;
        float* h_points;
        float Cy;
        float Cx;
        bool external_edge_existtence=false;
        /// CUDA Specific stuff.
        // Host storage.
        float *body_mass;
        float2 *body_pos;
        float *fx, *fy, *fx_prev, *fy_prev;

        // Quick way to represent a graph on the GPU
        int *sources, *targets;

        // Pointers to device memory (all suffixed with 'l').
        int   *errl,  *sortl, *childl, *countl, *startl;
        int   *sourcesl, *targetsl;
        float *body_massl, *node_massl;
        float2 *body_posl, *node_posl;
        float *minxl, *minyl, *maxxl, *maxyl;
        float *fxl, *fyl, *fx_prevl, *fy_prevl;
        float *swgl, *etral;
        int max_degree;
        int mp_count; // Number of multiprocessors on GPU.
        int max_threads_per_block;
        int nnodes;
        int num_points;
        int num_polygons;
        float scaled;
        int nbodies;
        int nedges;
        int max_iterations2;
        float maximum_aread;
        void sendGraphToGPU();
        void sendLayoutToGPU();
        void retrieveLayoutFromGPU();
         float max_distance_corner;
        
    };
};


#endif /* RPGPUForceAtlas2_hpp */
