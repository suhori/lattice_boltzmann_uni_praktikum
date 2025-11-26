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

#include <cstdlib>
#include <iostream>

using namespace std;

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

    auto test1= std::mdspan<double, std::extents<std::size_t, 2,6>> (ptr.get());
    auto test = std::mdspan(ptr.get(),2,6);
    for (std::size_t i = 0; i < test.extent(0);i++) {
        for (std::size_t j= 0; j< test.extent(1);j++) {
            std::cout << test[i,j] << " ";
        }
        std::cout << std::endl;
    }
    */
    auto lbm = LBM();
    printf("Simulating Taylor-Green vortex decay\n");
    printf("      domain size: %ux%u\n",lbm.NX,lbm.NY);
    printf("               nu: %g\n",lbm.nu);
    printf("              tau: %g\n",lbm.tau);
    printf("            u_max: %g\n",lbm.u_max);
    printf("             rho0: %g\n",lbm.rho0);
    printf("        timesteps: %u\n",lbm.NSTEPS);
    printf("       save every: %u\n",lbm.NSAVE);
    printf("    message every: %u\n",lbm.NMSG);
    printf("\n");

    double bytesPerMiB = 1024.0*1024.0;
    double bytesPerGiB = 1024.0*1024.0*1024.0;

    /*
    double *f0_alt  = (double*) malloc(mem_size_0dir);
    double *f1_alt  = (double*) malloc(mem_size_n0dir);
    double *f2_alt  = (double*) malloc(mem_size_n0dir);
    double *rho_alt = (double*) malloc(mem_size_scalar);
    double *ux_alt  = (double*) malloc(mem_size_scalar);
    double *uy_alt  = (double*) malloc(mem_size_scalar);
    */
    // rho is a two dimensional field
    // ux and uy are two dimensional fields respectivly
    // the field f is of the form f[N_x][N_y][q]

    auto ptr_f0 = make_unique<double[]>(lbm.mem_size_0dir);
    auto ptr_f1 = make_unique<double[]>(lbm.mem_size_n0dir);
    auto ptr_f2 = make_unique<double[]>(lbm.mem_size_n0dir);
    auto ptr_rho =make_unique<double[]>(lbm.mem_size_scalar);
    auto ptr_ux = make_unique<double[]>(lbm.mem_size_scalar);
    auto ptr_uy = make_unique<double[]>(lbm.mem_size_scalar);


    size_t total_mem_bytes = lbm.mem_size_0dir + 2*lbm.mem_size_n0dir + 3*lbm.mem_size_scalar;

    /*if(f0 == NULL || f1 == NULL || f2 == NULL || rho == NULL || ux == NULL || uy == NULL)
    {
        fprintf(stderr,"Error: unable to allocate required memory (%.1f MiB).\n",total_mem_bytes/bytesPerMiB);
        exit(-1);
    }*/
// TODO init with static extend<> ?
    auto m = lbm.ndir-1;
    auto f0 = mdspan(ptr_f0.get(),lbm.NX,lbm.NY);
    auto f1 = mdspan(ptr_f1.get(),lbm.NX,lbm.NY,m);
    auto f2 = mdspan(ptr_f2.get(),lbm.NX,lbm.NY,m);
    auto rho = mdspan(ptr_rho.get(),lbm.NX,lbm.NY);
    auto ux = mdspan(ptr_ux.get(),lbm.NX,lbm.NY);
    auto uy = mdspan(ptr_uy.get(),lbm.NX,lbm.NY);
    // compute Taylor-Green flow at t=0 
    // to initialise rho, ux, uy fields.
    lbm.taylor_green(0,rho,ux,uy);
    
    // initialise f1 as equilibrium for rho, ux, uy
    lbm.init_equilibrium(f0,f1,rho,ux,uy);


    lbm.save_scalar("rho",rho,0);
    lbm.save_scalar("ux", ux, 0);
    lbm.save_scalar("uy", uy, 0);

    if(lbm.computeFlowProperties)
    {
        lbm.report_flow_properties(0,rho,ux,uy);
    }
    
    double start = seconds();
    
    // main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < lbm.NSTEPS; ++n)
    {
        bool save = (n+1)%lbm.NSAVE == 0;
        bool msg  = (n+1)%lbm.NMSG == 0;
        bool need_scalars = save || (msg && lbm.computeFlowProperties);
        
        // stream and collide from f1 storing to f2
        // optionally compute and save moments
        lbm.stream_collide_save(f0,f1,f2,rho,ux,uy,need_scalars);

        if(save)
        {
            lbm.save_scalar("rho",rho,n+1);
            lbm.save_scalar("ux", ux, n+1);
            lbm.save_scalar("uy", uy, n+1);
        }
        //TODO
        // swap pointers
        // double *temp = f1;
        // f1 = f2;
        // f2 = temp;
        //
        f1 = mdspan(ptr_f2.get(),lbm.NX,lbm.NY,m);
        f2 = mdspan(ptr_f1.get(),lbm.NX,lbm.NY,m);
        //ISSUE what happens with the old mdspan objects?
        if(msg)
        {
            if(lbm.computeFlowProperties)
            {
                lbm.report_flow_properties(n+1,rho,ux,uy);
            }
            
            if(!lbm.quiet)
                printf("completed timestep %d\n",n+1);
        }
    }
    double end = seconds();
    double runtime = end-start;

    size_t doubles_read = lbm.ndir; // per node every time step
    size_t doubles_written = lbm.ndir;
    size_t doubles_saved = 3; // per node every NSAVE time steps
    
    // note NX*NY overflows when NX=NY=65536
    size_t nodes_updated = lbm.NSTEPS*size_t(lbm.NX*lbm.NY);
    size_t nodes_saved   = (lbm.NSTEPS/lbm.NSAVE)*size_t(lbm.NX*lbm.NY);
    double speed = nodes_updated/(1e6*runtime);
    
    double bandwidth = (nodes_updated*(doubles_read + doubles_written)+nodes_saved*(doubles_saved))*sizeof(double)/(runtime*bytesPerGiB);
    
    printf(" ----- performance information -----\n");
    printf(" memory allocated: %.1f (MiB)\n",total_mem_bytes/bytesPerMiB);
    printf("        timesteps: %u\n",lbm.NSTEPS);
    printf("          runtime: %.3f (s)\n",runtime);
    printf("            speed: %.2f (Mlups)\n",speed);
    printf("        bandwidth: %.1f (GiB/s)\n",bandwidth);
    
    // deallocate memory
    // free(f0);  free(f1); free(f2);
    // free(rho); free(ux); free(uy);
    return 0;
}

