/* This code accompanies
 *   The Lattice Boltzmann Method: Principles and Practice
 *   T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
 *   ISBN 978-3-319-44649-3 (Electronic) 
 *        978-3-319-44647-9 (Print)
 *   http://www.springer.com/978-3-319-44647-9
 *
 * This code is provided under the MIT license. See LICENSE.txt.
 *
 * Author: Orest Shardt
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <version>


#  include <mdspan>
#include <ostream>

#include "seconds.h"
#include "LBM.h"

int main(int argc, char* argv[])
{
    // Example use of mdspan (C++23).

    /*
    auto ptr = std::unique_ptr<double[]>(new double[]{1,2,3,4,5,6,7,8,9,10,11,12});
    double a1[12] = {1,2,3,4,5,6,7,8,9,10,11,12};

    auto test= std::mdspan<double, std::extents<std::size_t, 2,6>> (ptr.get());

    for (std::size_t i = 0; i < test.extent(0);i++) {
        for (std::size_t j= 0; j< test.extent(1);j++) {
            std::cout << test[i,j] << " ";
        }
        std::cout << std::endl;
    }
    */

    printf("Simulating Taylor-Green vortex decay\n");
    printf("      domain size: %ux%u\n",NX,NY);
    printf("               nu: %g\n",nu);
    printf("              tau: %g\n",tau);
    printf("            u_max: %g\n",u_max);
    printf("             rho0: %g\n",rho0);
    printf("        timesteps: %u\n",NSTEPS);
    printf("       save every: %u\n",NSAVE);
    printf("    message every: %u\n",NMSG);
    printf("\n");
    
    double bytesPerMiB = 1024.0*1024.0;
    double bytesPerGiB = 1024.0*1024.0*1024.0;

    /*
    double *f0  = (double*) malloc(mem_size_0dir);
    double *f1  = (double*) malloc(mem_size_n0dir);
    double *f2  = (double*) malloc(mem_size_n0dir);
    double *rho = (double*) malloc(mem_size_scalar);
    double *ux  = (double*) malloc(mem_size_scalar);
    double *uy  = (double*) malloc(mem_size_scalar);
    */
    auto ptr_f = std::make_unique<double[]>(mem_size_0dir+2*mem_size_n0dir);
    auto f0 = std::mdspan<double,std::extent<size_t,2,3,4>>(ptr_f.get());
    auto ptr_f1 = std::make_unique<double[]>(mem_size_n0dir);
    auto ptr_f2 = std::make_unique<double[]>(mem_size_n0dir);
    auto ptr_rho = std::make_unique<double[]>(mem_size_scalar);
    auto ptr_ux = std::make_unique<double[]>(mem_size_scalar);
    auto ptr_uy = std::make_unique<double[]>(mem_size_scalar);

    size_t total_mem_bytes = mem_size_0dir + 2*mem_size_n0dir + 3*mem_size_scalar;
    
    if(f0 == NULL || f1 == NULL || f2 == NULL || rho == NULL || ux == NULL || uy == NULL)
    {
        fprintf(stderr,"Error: unable to allocate required memory (%.1f MiB).\n",total_mem_bytes/bytesPerMiB);
        exit(-1);
    }
    
    // compute Taylor-Green flow at t=0 
    // to initialise rho, ux, uy fields.
    taylor_green(0,rho,ux,uy);
    
    // initialise f1 as equilibrium for rho, ux, uy
    init_equilibrium(f0,f1,rho,ux,uy);
    
    save_scalar("rho",rho,0);
    save_scalar("ux", ux, 0);
    save_scalar("uy", uy, 0);
    if(computeFlowProperties)
    {
        report_flow_properties(0,rho,ux,uy);
    }
    
    double start = seconds();
    
    // main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < NSTEPS; ++n)
    {
        bool save = (n+1)%NSAVE == 0;
        bool msg  = (n+1)%NMSG == 0;
        bool need_scalars = save || (msg && computeFlowProperties);
        
        // stream and collide from f1 storing to f2
        // optionally compute and save moments
        stream_collide_save(f0,f1,f2,rho,ux,uy,need_scalars);
        
        if(save)
        {
            save_scalar("rho",rho,n+1);
            save_scalar("ux", ux, n+1);
            save_scalar("uy", uy, n+1);
        }
        
        // swap pointers
        double *temp = f1;
        f1 = f2;
        f2 = temp;
        
        if(msg)
        {
            if(computeFlowProperties)
            {
                report_flow_properties(n+1,rho,ux,uy);
            }
            
            if(!quiet)
                printf("completed timestep %d\n",n+1);
        }
    }
    double end = seconds();
    double runtime = end-start;

    size_t doubles_read = ndir; // per node every time step
    size_t doubles_written = ndir;
    size_t doubles_saved = 3; // per node every NSAVE time steps
    
    // note NX*NY overflows when NX=NY=65536
    size_t nodes_updated = NSTEPS*size_t(NX*NY);
    size_t nodes_saved   = (NSTEPS/NSAVE)*size_t(NX*NY);
    double speed = nodes_updated/(1e6*runtime);
    
    double bandwidth = (nodes_updated*(doubles_read + doubles_written)+nodes_saved*(doubles_saved))*sizeof(double)/(runtime*bytesPerGiB);
    
    printf(" ----- performance information -----\n");
    printf(" memory allocated: %.1f (MiB)\n",total_mem_bytes/bytesPerMiB);
    printf("        timesteps: %u\n",NSTEPS);
    printf("          runtime: %.3f (s)\n",runtime);
    printf("            speed: %.2f (Mlups)\n",speed);
    printf("        bandwidth: %.1f (GiB/s)\n",bandwidth);
    
    // deallocate memory
    free(f0);  free(f1); free(f2);
    free(rho); free(ux); free(uy);
    
    return 0;
}

