/*
 ==============================================================================

 RPFA2Kernels.cuh
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

#ifndef RPFA2Kernels_cuh
#define RPFA2Kernels_cuh

#include "RPBHFA2LaunchParameters.cuh"

__global__
__launch_bounds__(THREADS6, FACTOR6)
void GravityKernel(int nbodiesd, const float k_g, const bool strong_gravity,
                   volatile float * __restrict body_massd,
                   volatile float2 * __restrict body_posd,
                   volatile float * __restrict fxd, volatile float * __restrict fyd,float Cx,float Cy,float* d_points,int* d_externa_polygonX,float* d_externa_polygonX_x,float* d_externa_polygonX_y,int max_degree,const float k_rd,int num_points);
__global__
__launch_bounds__(THREADS6, FACTOR6)
void exGravityKernel2(int nbodiesd, const float k_g, const bool strong_gravity,
                   volatile float* __restrict body_massd,
                   volatile float2* __restrict body_posd,
                   volatile float* __restrict fxd, volatile float* __restrict fyd,
                   float Cx, float Cy, float* d_points, int* d_externa_polygonX,
                   float* d_externa_polygonX_x, float* d_externa_polygonX_y,
                   int max_degree, int num_polygons, int num_points,float* d_effective_d);
__global__
__launch_bounds__(THREADS6, FACTOR6)
void cornerKernel(int nbodiesd, const float k_g, const bool strong_gravity,
                   volatile float* __restrict body_massd,
                   volatile float2* __restrict body_posd,
                   volatile float* __restrict fxd, volatile float* __restrict fyd,
                   float Cx, float Cy, float* d_points, int* d_externa_polygonX,
                   float* d_externa_polygonX_x, float* d_externa_polygonX_y,
                   int max_degree, int num_polygons, int num_points,float* d_effective_d,float max_distance_corner);
  
  
__global__
__launch_bounds__(THREADS6, FACTOR6)
void AttractiveForceKernel(int nedgesd,int nbodiesd,volatile float * __restrict body_massd,float maximum_area,
                           volatile float2 * __restrict body_posd,
                           volatile float * __restrict fxd, volatile float * __restrict fyd,
                           volatile int * __restrict sourcesd, volatile int * __restrict targetsd,float scale);
__global__
__launch_bounds__(THREADS6, FACTOR6)
void AttractiveForceKernel2(int nedgesd,int nbodiesd,volatile float * __restrict body_massd,float maximum_area,
                           volatile float2 * __restrict body_posd,
                           volatile float * __restrict fxd, volatile float * __restrict fyd,
                           volatile int * __restrict sourcesd, volatile int * __restrict targetsd,float scale, int* d_externa_polygonX,float* d_externa_polygonX_x,float* d_externa_polygonX_y);

__global__
__launch_bounds__(THREADS1, FACTOR1)
void SpeedKernel(int nbodiesd,
                 volatile float * __restrict fxd , volatile float * __restrict fyd,
                 volatile float * __restrict fx_prevd , volatile float * __restrict fy_prevd,
                 volatile float * __restrict body_massd, volatile float * __restrict swgd, volatile float * __restrict etrad);

__global__
__launch_bounds__(THREADS6, FACTOR6)
void DisplacementKernel(int nbodiesd,
                       volatile float2 * __restrict body_posd,
                       volatile float * __restrict fxd, volatile float * __restrict fyd,
                       volatile float * __restrict fx_prevd, volatile float * __restrict fy_prevd,float* d_points);

#endif
