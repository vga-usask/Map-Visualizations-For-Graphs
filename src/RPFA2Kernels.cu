
#include <cstdio>
#include <stdio.h>
#include "RPFA2Kernels.cuh"
#include "RPBHFA2LaunchParameters.cuh"
#include <cfloat>
/// Some variables for FA2 related to `speed'
static __device__ float k_s_maxd = 10.0;
static __device__ float global_speedd = 1.0;
static __device__ float speed_efficiencyd = 1.0;
static __device__ float jitter_toleranced = 1.0;
static __device__ unsigned int blkcntd_speed_kernel = 0;
__global__
__launch_bounds__(THREADS6, FACTOR6) void GravityKernel(int nbodiesd, const float k_g, const bool strong_gravity,
                                                        volatile float *__restrict body_massd,
                                                        volatile float2 *__restrict body_posd,
                                                        volatile float *__restrict fxd, volatile float *__restrict fyd,
                                                        float Cx, float Cy, float *d_points, int *d_externa_polygonX, float *d_externa_polygonX_x, float *d_externa_polygonX_y, int max_degree, const float k_r, int num_points)
{
    register int i, inc;

    // iterate over all bodies assigned to thread
    inc = blockDim.x * gridDim.x;
    for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc)
    {
        if (i >= nbodiesd)
            continue;

        const float px = body_posd[i].x;
        const float py = body_posd[i].y;
        const float epsilon = 0.0000001;
        // Distance from centroid to node (d_n)
        const float dx_c = px - Cx;
        const float dy_c = py - Cy;
        const float d_n = sqrtf(dx_c * dx_c + dy_c * dy_c);

        // Find d_c: the distance from the centroid to the polygon edge along the line passing through the node
        float d_c = FLT_MAX; // Initialize with a large number

        for (int j = 0; j < num_points; ++j)
        {
            float x1 = d_points[2 * j];
            float y1 = d_points[2 * j + 1];
            float x2 = d_points[2 * ((j + 1) % num_points)];
            float y2 = d_points[2 * ((j + 1) % num_points) + 1];

            // Vector from centroid to node (line direction)
            float dx_cn = dx_c;
            float dy_cn = dy_c;

            // Vector from edge start to end
            float edge_dx = x2 - x1;
            float edge_dy = y2 - y1;

            // Find intersection along the line from centroid to node with the polygon edge
            float denom = dx_cn * edge_dy - dy_cn * edge_dx;

            if (fabsf(denom) > epsilon)
            { // Ensure no division by zero (parallel check)
                float t = ((Cx - x1) * edge_dy - (Cy - y1) * edge_dx) / denom;

                // Intersection point along the centroid-node line (d_c)
                float intersect_x = Cx + t * dx_cn;
                float intersect_y = Cy + t * dy_cn;

                // Check if the intersection point is within the edge bounds
                float edge_t = ((intersect_x - x1) * edge_dx + (intersect_y - y1) * edge_dy) / ((edge_dx * edge_dx + edge_dy * edge_dy)+epsilon);
                if (edge_t >= 0.0f && edge_t <= 1.0f)
                {
                    float dist_to_edge = sqrtf((Cx - intersect_x) * (Cx - intersect_x) +
                                               (Cy - intersect_y) * (Cy - intersect_y));
                    d_c = fminf(d_c, dist_to_edge); // Keep the minimum distance to any edge
                }
            }
        }

        // Calculate d_n' as the difference: d_n' = d_c - d_n
        float d_n_prime = fabsf(d_c - d_n);

        // Compute gravitational force using the modified equation
        float f_g = 0.0f;

        if (d_n > d_c)
            f_g = k_g * body_massd[i] / (d_n - d_c);
        else if (d_n > (d_c * .7) && d_n != d_c)
            f_g = k_g * body_massd[i] / (d_c*.5);
                else if (d_n != 0)
                    f_g = k_g * body_massd[i] / (d_n);

        // f_g = k_g * body_massd[i]  / (d_n);

        // Apply gravitational force in the direction toward the centroid
        fxd[i] += (-dx_c * f_g);
        fyd[i] += (-dy_c * f_g);
    }
}
__global__ __launch_bounds__(THREADS6, FACTOR6) void exGravityKernel2(int nbodiesd, const float k_g, const bool strong_gravity,
                                                                      volatile float *__restrict body_massd,
                                                                      volatile float2 *__restrict body_posd,
                                                                      volatile float *__restrict fxd, volatile float *__restrict fyd,
                                                                      float Cx, float Cy, float *d_points, int *d_externa_polygonX,
                                                                      float *d_externa_polygonX_x, float *d_externa_polygonX_y,
                                                                      int max_degree, int num_polygons, int num_points, float *d_effective_d)
{
    register int i, inc;
    inc = blockDim.x * gridDim.x;
    const float epsilon = 1e-6f;

    for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc)
    {
        if (i >= nbodiesd)
            continue;

        const float px = body_posd[i].x;
        const float py = body_posd[i].y;
        float total_fx = 0.0f;
        float total_fy = 0.0f;
        float Cdx = px - Cx;
        float Cdy = py - Cy;
        // Check for external score and apply external forces if present
        for (int j = 0; j < 500; j++)
        {

            if (d_externa_polygonX[i * 500 + j] > 0)
            {
                float ex = d_externa_polygonX_x[j];
                float ey = d_externa_polygonX_y[j];
                float edx = px - ex;
                float edy = py - ey;
                if (sqrt(edx * edx + edy * edy)  > 0)
                {
                    float f_g = k_g *body_massd[i]*2  / ((sqrtf(edx * edx + edy * edy) *(d_externa_polygonX[i * 500 + j]) ));
                    if (body_massd[i] <= 1)
                        f_g = f_g * sqrtf(max_degree);
                    total_fx += -edx * f_g;
                    total_fy += -edy * f_g;
                }
            }
        }

        // Apply accumulated forces
        fxd[i] += total_fx;
        fyd[i] += total_fy;
    }
}
__global__ __launch_bounds__(THREADS6, FACTOR6) void cornerKernel(int nbodiesd, const float k_g, const bool strong_gravity,
                                                                  volatile float *__restrict body_massd,
                                                                  volatile float2 *__restrict body_posd,
                                                                  volatile float *__restrict fxd, volatile float *__restrict fyd,
                                                                  float Cx, float Cy, float *d_points, int *d_externa_polygonX,
                                                                  float *d_externa_polygonX_x, float *d_externa_polygonX_y,
                                                                  int max_degree, int num_polygons, int num_points, float *d_effective_d, float max_distance_corner)
{
    register int i, inc;
    inc = blockDim.x * gridDim.x;
    const float epsilon = 1e-6f;

    for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc)
    {
        if (i >= nbodiesd)
            continue;

        const float px = body_posd[i].x;
        const float py = body_posd[i].y;
        float total_fx = 0.0f;
        float total_fy = 0.0f;
        float Cdx = px - Cx;
        float Cdy = py - Cy;
        bool has_external_score = false;
        // Check for external score and apply external forces if present

        // Apply additional gravity forces toward sorted corner points
        // if(!has_external_score)
        for (int j = 0; j < num_points; j++)
        {

            int corner_idx = j; // centroid_indices[j];
            float effective_d = d_effective_d[j];
            float corner_x = d_points[2 * corner_idx];
            float corner_y = d_points[2 * corner_idx + 1];
            float dx = px - corner_x;
            float dy = py - corner_y;
            float dist = sqrtf((dx * dx + dy * dy));
            if (dist > effective_d&&dist>0)
            {
                float f_g_extra = (k_g*body_massd[i]) / (dist);
                total_fx += -dx * f_g_extra;
                total_fy += -dy * f_g_extra;
            }
        }

        // Apply accumulated forces
        fxd[i] += total_fx;
        fyd[i] += total_fy;
    }
}
__global__
__launch_bounds__(THREADS6, FACTOR6) void AttractiveForceKernel2(int nedgesd, int nbodiesd, volatile float *__restrict body_massd, float maximum_area,
                                                                 volatile float2 *__restrict body_posd,
                                                                 volatile float *__restrict fxd, volatile float *__restrict fyd,
                                                                 volatile int *__restrict sourcesd, volatile int *__restrict targetsd, float scale, int *d_externa_polygonX, float *d_externa_polygonX_x, float *d_externa_polygonX_y)
{
    int i, inc;
    inc = blockDim.x * gridDim.x;

    for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nedgesd; i += inc)
    {
        if (i >= nedgesd)
            continue;
        if(nbodiesd>0){
        const float min_distance = sqrt((maximum_area) / nbodiesd);
        int source = sourcesd[i];
        int target = targetsd[i];

        // dx and dy are distance between source and target nodes.
        float dx = body_posd[target].x - body_posd[source].x;
        float dy = body_posd[target].y - body_posd[source].y;

        // Current distance between the nodes
        float distance = sqrtf(dx * dx + dy * dy);

        // Desired distance after scaling
        float desired_distance = distance + (min_distance );

        // Adjust positions proportionally to achieve the desired distance
        if (distance > 0)
        { // Avoid division by zero
            float adjustment_factor = (desired_distance - distance) / distance;
            float force_x = dx * adjustment_factor;
            float force_y = dy * adjustment_factor;

            atomicAdd((float *)fxd + source, -force_x);
            atomicAdd((float *)fyd + source, -force_y);
            atomicAdd((float *)fxd + target, force_x);
            atomicAdd((float *)fyd + target, force_y);
        }
    }}
}
__global__
__launch_bounds__(THREADS6, FACTOR6) void AttractiveForceKernel(int nedgesd, int nbodiesd, volatile float *__restrict body_massd, float maximum_area,
                                                                volatile float2 *__restrict body_posd,
                                                                volatile float *__restrict fxd, volatile float *__restrict fyd,
                                                                volatile int *__restrict sourcesd, volatile int *__restrict targetsd, float scale)
{
    register int i, inc, source, target;
    // iterate over all edges assigned to thread
    inc = blockDim.x * gridDim.x;
    for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nedgesd; i += inc)
    {
        if (i >= nedgesd)
            continue;

        source = sourcesd[i];
        target = targetsd[i];

        // dx and dy are distance to between the neighbors.
        const float dx = body_posd[target].x - body_posd[source].x;
        const float dy = body_posd[target].y - body_posd[source].y;

        // Calculate the Euclidean distance between the source and target nodes

        const float fsx = dx;
        const float fsy = dy;
        const float ftx = -dx;
        const float fty = -dy;

        atomicAdd((float *)fxd + source, fsx);
        atomicAdd((float *)fyd + source, fsy);
        atomicAdd((float *)fxd + target, ftx);
        atomicAdd((float *)fyd + target, fty);
    }
}

__global__
__launch_bounds__(THREADS1, FACTOR1) void SpeedKernel(int nbodiesd,
                                                      volatile float *__restrict fxd, volatile float *__restrict fyd,
                                                      volatile float *__restrict fx_prevd, volatile float *__restrict fy_prevd,
                                                      volatile float *__restrict body_massd, volatile float *__restrict swgd, volatile float *__restrict etrad)
{
    register int i, j, k, inc;
    register float swg_thread, swg_body, etra_thread, etra_body, dx, dy, mass;
    // setra: effective_traction (in shared mem.)
    // sswg: swing per node (in shared mem.)
    __shared__ volatile float sswg[THREADS1], setra[THREADS1];

    // initialize with valid data (in case #bodies < #threads)
    swg_thread = 0;
    etra_thread = 0;

    // scan all bodies
    i = threadIdx.x;
    inc = THREADS1 * gridDim.x;

    for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc)
    {
        if (i >= nbodiesd)
            continue;
        mass = body_massd[j];

        dx = fxd[j] - fx_prevd[j];
        dy = fyd[j] - fy_prevd[j];
        swg_body = sqrtf(dx * dx + dy * dy);
        swg_thread += mass * swg_body;

        dx = fxd[j] + fx_prevd[j];
        dy = fyd[j] + fy_prevd[j];
        etra_body = sqrtf(dx * dx + dy * dy) / 2.0;
        etra_thread += mass * etra_body;
    }

    // reduction in shared memory
    sswg[i] = swg_thread;
    setra[i] = etra_thread;

    for (j = THREADS1 / 2; j > 0; j /= 2)
    {
        __syncthreads();
        if (i < j)
        {
            k = i + j;
            sswg[i] = swg_thread = sswg[i] + sswg[k];
            setra[i] = etra_thread = setra[i] + setra[k];
        }
    }

    // swg_thread and etra_thread are now the total swinging
    // and the total effective traction (accross all threads)

    // write block result to global memory
    if (i == 0)
    {
        k = blockIdx.x;
        swgd[k] = swg_thread;
        etrad[k] = etra_thread;
        __threadfence();

        inc = gridDim.x - 1;
        if (inc == atomicInc(&blkcntd_speed_kernel, inc))
        {
            swg_thread = 0;
            etra_thread = 0;

            for (j = 0; j <= inc; j++)
            {
                swg_thread += swgd[j];
                etra_thread += etrad[j];
            }
            // we need to do some calculations to derive
            // from this the new global speed
            float estimated_optimal_jitter_tollerance = 0.05 * sqrtf(nbodiesd);
            float minJT = sqrtf(estimated_optimal_jitter_tollerance);
            float jt = jitter_toleranced * fmaxf(minJT, fminf(k_s_maxd, estimated_optimal_jitter_tollerance * etra_thread / powf(nbodiesd, 2.0)));
            float min_speed_efficiency = 0.05;

            // `Protect against erratic behavior'
            if (swg_thread / etra_thread > 2.0)
            {
                if (speed_efficiencyd > min_speed_efficiency)
                    speed_efficiencyd *= 0.5;
                jt = fmaxf(jt, jitter_toleranced);
            }

            // `Speed efficiency is how the speed really corrosponds to the swinging vs. convergence tradeoff.'
            // `We adjust it slowly and carefully'
            float targetSpeed = jt * speed_efficiencyd * etra_thread / swg_thread;

            if (swg_thread > jt * etra_thread)
            {
                if (speed_efficiencyd > min_speed_efficiency)
                {
                    speed_efficiencyd *= 0.7;
                }
            }
            else if (global_speedd < 1000)
            {
                speed_efficiencyd *= 1.3;
            }

            // `But the speed shouldn't rise much too quickly, ... would make convergence drop dramatically'.
            float max_rise = 0.5;
            global_speedd += fminf(targetSpeed - global_speedd, max_rise * global_speedd);
        }
    }
}

__global__
__launch_bounds__(THREADS6, FACTOR6) void DisplacementKernel(int nbodiesd,
                                                             volatile float2 *__restrict body_posd,
                                                             volatile float *__restrict fxd, volatile float *__restrict fyd,
                                                             volatile float *__restrict fx_prevd, volatile float *__restrict fy_prevd,
                                                             float *d_points)
{
    register int i, inc;
    register float factor, swg, dx, dy, fx, fy;
    register float global_speed = global_speedd;

    // iterate over all bodies assigned to thread
    inc = blockDim.x * gridDim.x;
    for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc)
    {
        if (i >= nbodiesd)
            continue;

        fx = fxd[i];
        fy = fyd[i];
        dx = fx - fx_prevd[i];
        dy = fy - fy_prevd[i];
        swg = sqrtf(dx * dx + dy * dy);
        factor = global_speed / (1.0 + sqrtf(global_speed * swg));

        // Update the position
        float new_x = body_posd[i].x + fx * factor;
        float new_y = body_posd[i].y + fy * factor;

        // Only update position if inside the polygon

        body_posd[i].x = new_x;
        body_posd[i].y = new_y;

        // Save previous forces
        fx_prevd[i] = fx;
        fy_prevd[i] = fy;
        fxd[i] = 0.0;
        fyd[i] = 0.0;
    }
}
