#include <stdio.h>
#include <fstream>
#include <chrono>
#include <algorithm>
#include "time.h"
#include "DataStore.hpp"
#include <vector>
#include <cmath>
#include "RPGPUForceAtlas2.hpp"
#include "RPBHFA2LaunchParameters.cuh"
#include "RPBHKernels.cuh"
#include "RPFA2Kernels.cuh"
#include <iostream>

namespace RPGraph
{
    CUDAForceAtlas2::CUDAForceAtlas2(GraphLayout &layout, bool use_barneshut,
                                     bool strong_gravity, float gravity,
                                     float scale, int max_iterations,std::vector<std::vector<int>> &externa_polygonX,float maximum_area)
    : ForceAtlas2(layout, use_barneshut, strong_gravity, gravity, scale,  max_iterations)
    {
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        if (deviceCount == 0)
        {
            fprintf(stderr, "error: No CUDA devices found.\n");
            exit(EXIT_FAILURE);
        }

        // Host initialization and setup //
        nbodies = layout.graph.num_nodes();
        nedges  = layout.graph.num_edges();
        h_points = layout.getPolygonPoints();
        body_pos = (float2 *)malloc(sizeof(float2) * layout.graph.num_nodes());
        body_mass = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        sources  = (int *)  malloc(sizeof(int)   * layout.graph.num_edges());
        targets  = (int *)  malloc(sizeof(int)   * layout.graph.num_edges());
        fx       = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        fy       = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        fx_prev  = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        fy_prev  = (float *)malloc(sizeof(float) * layout.graph.num_nodes());
        Cy=layout.getCy();
        Cx=layout.getCx();
        max_iterations2= max_iterations;
        DataStore& dataStore = DataStore::getInstance();
        num_points=dataStore.polygons[layout.polygonId].size();
        num_polygons=dataStore.number_of_polygons;
        std::cout<<"number of the points"<<num_points<<std::endl;
        std::cout<<"number of the polygons"<<num_polygons<<std::endl;
        std::cout<<"number of the nodes"<<nbodies<<std::endl;
         std::cout<<"polygon ID"<<layout.polygonId<<std::endl;
    // Calculate effective distances for each corner
    float effective_d[num_points];
   
    for (int k = 0; k < num_points; k++) {
        float corner_x = h_points[ 2 * k];
        float corner_y = h_points[1 + 2 * k];
        float corner_x_next = h_points[1+( 1+2 * k)];
        float corner_y_next = h_points[1+1+(1 + 2 * k)];
        float dx_next = corner_x_next - corner_x;
        float dy_next = corner_y_next- corner_y;
        float d_next=sqrtf(dx_next * dx_next + dy_next * dy_next) ;
        float corner_x_rev = 0;
        float corner_y_rev = 0;
         if(k==0)
        {
             corner_x_rev = h_points[ num_points-2];
             corner_y_rev = h_points[num_points-1];
        }
        else
        {
             corner_x_rev = h_points[ 2-2 * k];
             corner_y_rev = h_points[2-1 + 2 * k];
        }
        float dx_rev = corner_x_rev - corner_x;
        float dy_rev = corner_y_rev - corner_y;
        float d_rev=sqrtf(dx_rev * dx_rev + dy_rev * dy_rev) ;
        float dxC=corner_x-Cx ;
        float dyC=corner_y-Cy ;
        float dC=sqrtf(dxC * dxC + dyC * dyC) ;
        if(dC>max_distance_corner)
        max_distance_corner=dC;
        effective_d[k] = std::min({d_next,d_rev})/2 ;//for now I need to think more about the Dc
    }
    //*************** */
        int numRows = 100000;
        int numCols = 500;
        std::vector<int> flattened(numRows * numCols);
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {

                flattened[i * numCols + j] = externa_polygonX[layout.graph.node_map_r[i]][j];
                
                if(externa_polygonX[layout.graph.node_map_r[i]][j]>0)
                external_edge_existtence=true;
            }
        }        
        float externa_polygonX_x[500];
        float externa_polygonX_y[500];
        scaled=scale;
        maximum_aread=maximum_area;
        
       for(int i=0;i<500;i++)
       {
            externa_polygonX_x[i]=dataStore.polygonCenters[i].first;
            externa_polygonX_y[i]=dataStore.polygonCenters[i].second;
            
    
       }

      std::cout<<"number of the polygons**"<<num_polygons<<std::endl;
        cudaCatchError(cudaMalloc((void**)&d_externa_polygonX_x, 500 * sizeof(float)));
        cudaMemset(d_externa_polygonX_x, 0, 500 * sizeof(float));
        cudaCatchError(cudaMemcpy(d_externa_polygonX_x, externa_polygonX_x, 500 * sizeof(float), cudaMemcpyHostToDevice));
        
        cudaCatchError(cudaMalloc((void**)&d_externa_polygonX_y, 500 * sizeof(float)));
        cudaMemset(d_externa_polygonX_y, 0, 500 * sizeof(float));
        cudaCatchError(cudaMemcpy(d_externa_polygonX_y, externa_polygonX_y, 500 * sizeof(float), cudaMemcpyHostToDevice));

        cudaCatchError(cudaMalloc((void**)&d_externa_polygonX, numRows * numCols * sizeof(int)));
        cudaMemset(d_externa_polygonX, 0, numRows * numCols * sizeof(float));
        cudaCatchError(cudaMemcpy(d_externa_polygonX, flattened.data(), numRows * numCols * sizeof(int), cudaMemcpyHostToDevice));
        
        cudaCatchError(cudaMalloc((void**)&d_effective_d, num_points * sizeof(float)));
        cudaMemset(d_effective_d, 0, num_points * sizeof(float));
        cudaCatchError(cudaMemcpy(d_effective_d, effective_d, num_points * sizeof(float), cudaMemcpyHostToDevice));
        
        cudaCatchError(cudaMalloc((void**)&d_points, num_points*2 * sizeof(float)));
        cudaMemset(d_points, 0, num_points*2 * sizeof(float));
        cudaCatchError(cudaMemcpy(d_points, h_points, num_points*2 * sizeof(float), cudaMemcpyHostToDevice));



      

        max_degree=0;

        for (nid_t n = 0; n < layout.graph.num_nodes(); ++n)
        {
            body_pos[n] = {layout.getX(n), layout.getY(n)};
            body_mass[n] = ForceAtlas2::mass(n);
      
            if( body_mass[n] >max_degree)
            max_degree= body_mass[n] ;
            fx[n] = 0.0;
            fy[n] = 0.0;
            fx_prev[n] = 0.0;
            fy_prev[n] = 0.0;
        }
printf("MAX DEGREE IS %d\n", max_degree);
        int cur_sources_idx = 0;
        int cur_targets_idx = 0;

        // Initialize the sources and targets arrays with edge-data.
        for (nid_t source_id = 0; source_id < layout.graph.num_nodes(); ++source_id)
        {
            for (nid_t target_id : layout.graph.neighbors_with_geq_id(source_id))
            {
                sources[cur_sources_idx++] = source_id;
                targets[cur_targets_idx++] = target_id;
            }
        }

        // GPU initialization and setup //
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, 0);

        if (deviceProp.warpSize != WARPSIZE)
        {
            printf("Warpsize of device is %d, but we anticipated %d\n", deviceProp.warpSize, WARPSIZE);
            exit(EXIT_FAILURE);

        }
        cudaFuncSetCacheConfig(BoundingBoxKernel, cudaFuncCachePreferShared);
        cudaFuncSetCacheConfig(TreeBuildingKernel, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(ClearKernel1, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
        cudaFuncSetCacheConfig(SummarizationKernel, cudaFuncCachePreferShared);
        cudaFuncSetCacheConfig(SortKernel, cudaFuncCachePreferL1);
#if __CUDA_ARCH__ < 300
        cudaFuncSetCacheConfig(ForceCalculationKernel, cudaFuncCachePreferL1);
#endif
        cudaFuncSetCacheConfig(DisplacementKernel, cudaFuncCachePreferL1);

        cudaGetLastError();  // reset error value

        // Allocate space on device.
        mp_count = deviceProp.multiProcessorCount;
        max_threads_per_block = deviceProp.maxThreadsPerBlock;

        nnodes = std::max(2 * nbodies, mp_count * max_threads_per_block);

        // Round up to next multiple of WARPSIZE
        while ((nnodes & (WARPSIZE-1)) != 0) nnodes++;
        nnodes--;

        // child stores structure of the quadtree. values point to IDs.
        cudaCatchError(cudaMalloc((void **)&childl,  sizeof(int)   * (nnodes+1) * 4));

        // the following properties, for each node in the quadtree (both internal and leaf)
        cudaCatchError(cudaMalloc((void **)&body_massl,   sizeof(float) * nbodies));
        cudaCatchError(cudaMalloc((void **)&node_massl,   sizeof(float) * (nnodes+1)));
        cudaCatchError(cudaMalloc((void **)&body_posl,sizeof(float2) * nbodies));
        cudaCatchError(cudaMalloc((void **)&node_posl,    sizeof(float2) * (nnodes+1)));
        // count contains the number of nested nodes for each node in quadtree
        cudaCatchError(cudaMalloc((void **)&countl,  sizeof(int)   * (nnodes+1)));
        // start contains ...
        cudaCatchError(cudaMalloc((void **)&startl,  sizeof(int)   * (nnodes+1)));
        cudaCatchError(cudaMalloc((void **)&sortl,   sizeof(int)   * (nnodes+1)));


        cudaCatchError(cudaMalloc((void **)&sourcesl,sizeof(int)   * (nedges)));
        cudaCatchError(cudaMalloc((void **)&targetsl,sizeof(int)   * (nedges)));
        cudaCatchError(cudaMalloc((void **)&fxl,     sizeof(float) * (nbodies)));
        cudaCatchError(cudaMalloc((void **)&fyl,     sizeof(float) * (nbodies)));
        cudaCatchError(cudaMalloc((void **)&fx_prevl,sizeof(float) * (nbodies)));
        cudaCatchError(cudaMalloc((void **)&fy_prevl,sizeof(float) * (nbodies)));

        // Used for reduction in BoundingBoxKernel
        cudaCatchError(cudaMalloc((void **)&maxxl,   sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&maxyl,   sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&minxl,   sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&minyl,   sizeof(float) * mp_count * FACTOR1));

        // Used for reduction in SpeedKernel
        cudaCatchError(cudaMalloc((void **)&swgl,    sizeof(float) * mp_count * FACTOR1));
        cudaCatchError(cudaMalloc((void **)&etral,   sizeof(float) * mp_count * FACTOR1));

        // Copy host data to device.
        cudaCatchError(cudaMemcpy(body_massl, body_mass, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(body_posl,  body_pos,  sizeof(float2) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(sourcesl, sources, sizeof(int) * nedges, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(targetsl, targets, sizeof(int) * nedges, cudaMemcpyHostToDevice));

        // cpy fx, fy , fx_prevl, fy_prevl so they are all initialized to 0 in device memory.
        cudaCatchError(cudaMemcpy(fxl, fx,           sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(fyl, fy,           sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(fx_prevl, fx_prev, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(fy_prevl, fy_prev, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
    }

    void CUDAForceAtlas2::freeGPUMemory()
    {
        cudaDeviceSynchronize();
        cudaFree(childl);

        cudaFree(body_massl);
        cudaFree(node_massl);
        cudaFree(body_posl);
        cudaFree(node_posl);
        cudaFree(sourcesl);
        cudaFree(targetsl);
        cudaFree(countl);
        cudaFree(startl);
        cudaFree(sortl);

        cudaFree(fxl);
        cudaFree(fx_prevl);
        cudaFree(fyl);
        cudaFree(fy_prevl);
        cudaFree(fy_prevl);


        cudaFree(maxxl);
        cudaFree(maxyl);
        cudaFree(minxl);
        cudaFree(minyl);

        cudaFree(swgl);
        cudaFree(etral);
        cudaFree(d_points);
        cudaFree(d_externa_polygonX);
        cudaFree(d_externa_polygonX_x);
        cudaFree(d_externa_polygonX_y);
        cudaFree(d_effective_d);
        

   
    }

    CUDAForceAtlas2::~CUDAForceAtlas2()
    {

      
        free(body_mass);
        free(body_pos);
        free(sources);
        free(targets);
        free(fx);
        free(fy);
        free(fx_prev);
        free(fy_prev);
        free(h_points);

    
    }

 void CUDAForceAtlas2::doStep(int inter)
{
 
    cudaGetLastError(); // clear any errors
        exGravityKernel2<<<mp_count * FACTOR6, THREADS6>>>(nbodies, k_g,strong_gravity, body_massl, body_posl,fxl, fyl,Cx,  Cy,  d_points, d_externa_polygonX,d_externa_polygonX_x,  d_externa_polygonX_y,max_degree,  num_polygons,  num_points, d_effective_d);
        cudaCatchError(cudaGetLastError());
         cornerKernel<<<mp_count * FACTOR6, THREADS6>>>(nbodies, k_g,strong_gravity, body_massl, body_posl,fxl, fyl,Cx,  Cy,  d_points, d_externa_polygonX,d_externa_polygonX_x,  d_externa_polygonX_y,max_degree,  num_polygons,  num_points, d_effective_d,max_distance_corner);
        cudaCatchError(cudaGetLastError());
        GravityKernel<<<mp_count * FACTOR6, THREADS6>>>(nbodies, k_g, strong_gravity, body_massl, body_posl, fxl, fyl, Cx, Cy,d_points,d_externa_polygonX,d_externa_polygonX_x,d_externa_polygonX_y,max_degree,k_r,num_points);
        cudaCatchError(cudaGetLastError());
        
       
       AttractiveForceKernel<<<mp_count * FACTOR6, THREADS6>>>(nedges, nbodies,body_massl, maximum_aread, body_posl, fxl, fyl, sourcesl, targetsl,scaled);
        cudaCatchError(cudaGetLastError());
         AttractiveForceKernel2<<<mp_count * FACTOR6, THREADS6>>>(nedges, nbodies,body_massl, maximum_aread, body_posl, fxl, fyl, sourcesl, targetsl,scaled,d_externa_polygonX,d_externa_polygonX_x,d_externa_polygonX_y);
        cudaCatchError(cudaGetLastError());
        BoundingBoxKernel<<<mp_count * FACTOR1, THREADS1>>>(nnodes, nbodies, startl, childl, node_massl, body_posl, node_posl, maxxl, maxyl, minxl, minyl);
        cudaCatchError(cudaGetLastError());
		
        // Build Barnes-Hut Tree
        ClearKernel1<<<mp_count, 1024>>>(nnodes, nbodies, childl);
        cudaCatchError(cudaGetLastError());
        TreeBuildingKernel<<<mp_count * FACTOR2, THREADS2>>>(nnodes, nbodies, childl, body_posl, node_posl);
        cudaCatchError(cudaGetLastError());
        ClearKernel2<<<mp_count, 1024>>>(nnodes, startl, node_massl);
        cudaCatchError(cudaGetLastError());
        SummarizationKernel<<<mp_count * FACTOR3, THREADS3>>>(nnodes, nbodies, countl, childl, body_massl, node_massl, body_posl, node_posl);
        cudaCatchError(cudaGetLastError());
        SortKernel<<<mp_count * FACTOR4, THREADS4>>>(nnodes, nbodies, sortl, countl, startl, childl);
        cudaCatchError(cudaGetLastError());
        ForceCalculationKernel<<<mp_count * FACTOR5, THREADS5>>>(nnodes, nbodies, itolsq, epssq, sortl, childl, body_massl, node_massl, body_posl, node_posl, fxl, fyl, k_r);
        cudaCatchError(cudaGetLastError());
   
                
       
        

        SpeedKernel<<<mp_count * FACTOR1, THREADS1>>>(nbodies, fxl, fyl, fx_prevl, fy_prevl, body_massl, swgl, etral);
        cudaCatchError(cudaGetLastError());

        DisplacementKernel<<<mp_count * FACTOR6, THREADS6>>>(nbodies, body_posl, fxl, fyl, fx_prevl, fy_prevl,d_points);
        cudaCatchError(cudaGetLastError());

    cudaDeviceSynchronize();
    
     
}



    void CUDAForceAtlas2::retrieveLayoutFromGPU()
    {
        cudaCatchError(cudaMemcpy(body_pos, body_posl, sizeof(float2) * nbodies, cudaMemcpyDeviceToHost));
        cudaDeviceSynchronize();
    }

    void CUDAForceAtlas2::sendLayoutToGPU()
    {
        cudaCatchError(cudaMemcpy(body_posl, body_pos, sizeof(float2) * nbodies, cudaMemcpyHostToDevice));
        cudaDeviceSynchronize();
    }

    void CUDAForceAtlas2::sendGraphToGPU()
    {
        cudaCatchError(cudaMemcpy(body_massl, body_mass, sizeof(float) * nbodies, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(sourcesl, sources, sizeof(int) * nedges, cudaMemcpyHostToDevice));
        cudaCatchError(cudaMemcpy(targetsl, targets, sizeof(int) * nedges, cudaMemcpyHostToDevice));
        //*************************
 


//**********************
        cudaDeviceSynchronize();
    }

    bool CUDAForceAtlas2::sync_layout()
    {
        bool return_flag=true;
        retrieveLayoutFromGPU();
       		 for(nid_t n = 0; n < layout.graph.num_nodes(); ++n)
       	 {
            bool result=layout.move_node(n,body_pos[n].x,body_pos[n].y);
       	 	if(result==false)
            {
       	 		return_flag= false;
                
            		//layout.setX(n, body_pos[n].x);
            		//layout.setY(n, body_pos[n].y);
            }
       	 }

        	return return_flag;
    }
     void CUDAForceAtlas2::sync_layout_void()
    {
        retrieveLayoutFromGPU();
       		 for(nid_t n = 0; n < layout.graph.num_nodes(); ++n)
       	 {
       	 			//layout.move_node_void(n,body_pos[n].x,body_pos[n].y);
       	 					
            		layout.setX(n, body_pos[n].x);
            		layout.setY(n, body_pos[n].y);
       	 }
        	
    }
}
