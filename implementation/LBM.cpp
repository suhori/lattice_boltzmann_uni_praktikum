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
#include <cstdio>
#define _USE_MATH_DEFINES
#include <iostream>
//#include <math.h>
#include <cmath>
#include <mdspan>
#include <vector>
#include <memory>
#include "LBM.h"
using namespace std;

void LBM::taylor_green(unsigned int t, unsigned int x, unsigned int y,mdspan<double, dextents<size_t, 2>> r,mdspan<double, dextents<size_t, 2>> u,mdspan<double, dextents<size_t, 2>> v)
{
    double kx = 2.0*M_PI/NX;
    double ky = 2.0*M_PI/NY;
    double td = 1.0/(nu*(kx*kx+ky*ky));
    
    double X = x+0.5;
    double Y = y+0.5;
    double ux = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1.0*t/td);
    double uy =  u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1.0*t/td);
    double P = -0.25*rho0*u_max*u_max*((ky/kx)*cos(2.0*kx*X)+(kx/ky)*cos(2.0*ky*Y))*exp(-2.0*t/td);
    double rho = rho0+3.0*P;
    
    r[x,y] = rho;
    u[x,y] = ux;
    v[x,y] = uy;
}

void LBM::taylor_green(unsigned int t, mdspan<double, dextents<size_t, 2>> r,mdspan<double, dextents<size_t, 2>> u,mdspan<double, dextents<size_t, 2>> v)
{
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            size_t sidx = scalar_index(x,y);
//TODO slice arrays
            this->taylor_green(t,x,y,r,u,v);
        }
    }
}
void LBM::taylor_green_cfp(unsigned int t, unsigned int x, unsigned int y, double *r, double *u, double *v)
{
    double kx = 2.0*M_PI/NX;
    double ky = 2.0*M_PI/NY;
    double td = 1.0/(nu*(kx*kx+ky*ky));

    double X = x+0.5;
    double Y = y+0.5;
    double ux = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1.0*t/td);
    double uy =  u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1.0*t/td);
    double P = -0.25*rho0*u_max*u_max*((ky/kx)*cos(2.0*kx*X)+(kx/ky)*cos(2.0*ky*Y))*exp(-2.0*t/td);
    double rho = rho0+3.0*P;

    *r = rho;
    *u = ux;
    *v = uy;
}
void LBM::init_equilibrium(mdspan<double, dextents<size_t, 2>> f0,mdspan<double, dextents<size_t, 3>> f1, mdspan<double, dextents<size_t, 2>> r,mdspan<double, dextents<size_t, 2>> u,mdspan<double, dextents<size_t, 2>> v)
{
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[x,y];
            double ux  = u[x,y];
            double uy  = v[x,y];
            
            // load equilibrium
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]
            
            // temporary variables
            double w0r = w0*rho;
            double wsr = ws*rho;
            double wdr = wd*rho;
            double omusq = 1.0 - 1.5*(ux*ux+uy*uy);
            
            double tux = 3.0*ux;
            double tuy = 3.0*uy;
            
            f0[x,y]    = w0r*(omusq);
            
            double cidot3u = tux;
            f1[x,y,1]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            f1[x,y,2]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            f1[x,y,3]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            f1[x,y,4]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            
            cidot3u = tux+tuy;
            f1[x,y,5]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            f1[x,y,6]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            f1[x,y,7]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            f1[x,y,8]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
        }
    }
}

void LBM::stream_collide_save(mdspan<double, dextents<size_t, 2>> f0, mdspan<double, dextents<size_t, 3>> f1, mdspan<double, dextents<size_t, 3>> f2, mdspan<double, dextents<size_t, 2>> r, mdspan<double, dextents<size_t, 2>> u, mdspan<double, dextents<size_t, 2>> v, bool save)
{
    // useful constants
    const double tauinv = 2.0/(6.0*nu+1.0); // 1/tau
    const double omtauinv = 1.0-tauinv;     // 1 - 1/tau

    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            unsigned int xp1 = (x+1)%NX;
            unsigned int yp1 = (y+1)%NY;
            unsigned int xm1 = (NX+x-1)%NX;
            unsigned int ym1 = (NY+y-1)%NY;
            
            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8
            
            double ft0 = f0[x,y];
            
            // load populations from adjacent nodes
            double ft1 = f1[xm1,y,  1];
            double ft2 = f1[x,  ym1,2];
            double ft3 = f1[xp1,y,  3];
            double ft4 = f1[x,  yp1,4];
            double ft5 = f1[xm1,ym1,5];
            double ft6 = f1[xp1,ym1,6];
            double ft7 = f1[xp1,yp1,7];
            double ft8 = f1[xm1,yp1,8];
            
            // compute moments
            double rho = ft0+ft1+ft2+ft3+ft4+ft5+ft6+ft7+ft8;
            double rhoinv = 1.0/rho;
            
            double ux = rhoinv*(ft1+ft5+ft8-(ft3+ft6+ft7));
            double uy = rhoinv*(ft2+ft5+ft6-(ft4+ft7+ft8));
            
            // only write to memory when needed
            if(save)
            {
                r[x,y] = rho;
                u[x,y] = ux;
                v[x,y] = uy;
            }
            
            // now compute and relax to equilibrium
            // note that
            // feq_i  = w_i rho [1 + 3(ci . u) + (9/2) (ci . u)^2 - (3/2) (u.u)]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u) + (1/2) (ci . 3u)^2]
            // feq_i  = w_i rho [1 - 3/2 (u.u) + (ci . 3u){ 1 + (1/2) (ci . 3u) }]
            
            // temporary variables
            double tw0r = tauinv*w0*rho; //   w[0]*rho/tau 
            double twsr = tauinv*ws*rho; // w[1-4]*rho/tau
            double twdr = tauinv*wd*rho; // w[5-8]*rho/tau
            double omusq = 1.0 - 1.5*(ux*ux+uy*uy); // 1-(3/2)u.u
            
            double tux = 3.0*ux;
            double tuy = 3.0*uy;
            
            
            f0[x,y]    = omtauinv*ft0  + tw0r*(omusq);
            
            double cidot3u = tux;
            f2[x,y,1]  = omtauinv*ft1  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            f2[x,y,2]  = omtauinv*ft2  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            f2[x,y,3]  = omtauinv*ft3  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            f2[x,y,4]  = omtauinv*ft4  + twsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            
            cidot3u = tux+tuy;
            f2[x,y,5]  = omtauinv*ft5  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            f2[x,y,6]  = omtauinv*ft6  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            f2[x,y,7]  = omtauinv*ft7  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            f2[x,y,8]  = omtauinv*ft8  + twdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
        }
    }
}


void LBM::compute_flow_properties(unsigned int t, mdspan<double, dextents<size_t, 2>> r,mdspan<double, dextents<size_t, 2>> u,mdspan<double, dextents<size_t, 2>> v, double *prop)
{
    // prop must point to space for 4 doubles:
    // 0: energy
    // 1: L2 error in rho
    // 2: L2 error in ux
    // 3: L2 error in uy
    
    double E = 0.0; // kinetic energy
    
    double sumrhoe2 = 0.0; // sum of error squared in rho
    double sumuxe2 = 0.0;  //                         ux
    double sumuye2 = 0.0;  //                         uy
    
    double sumrhoa2 = 0.0; // sum of analytical rho squared
    double sumuxa2 = 0.0;  //                   ux
    double sumuya2 = 0.0;  //                   uy
    
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[x,y];
            double ux  = u[x,y];
            double uy  = v[x,y];
            E += rho*(ux*ux + uy*uy);
            
            double rhoa, uxa, uya;
            taylor_green_cfp(t,x,y,&rhoa,&uxa,&uya);
            
            sumrhoe2 += (rho-rhoa)*(rho-rhoa);
            sumuxe2  += (ux-uxa)*(ux-uxa);
            sumuye2  += (uy-uya)*(uy-uya);

            sumrhoa2 += (rhoa-rho0)*(rhoa-rho0);
            sumuxa2  += uxa*uxa;
            sumuya2  += uya*uya;
        }
    }
    
    prop[0] = E;
    prop[1] = sqrt(sumrhoe2/sumrhoa2);
    prop[2] = sqrt(sumuxe2/sumuxa2);
    prop[3] = sqrt(sumuye2/sumuya2);
}

void LBM::report_flow_properties(unsigned int t, mdspan<double, dextents<size_t, 2>> rho,mdspan<double, dextents<size_t, 2>> ux,mdspan<double, dextents<size_t, 2>> uy)
{
    double prop[4];
    compute_flow_properties(t,rho,ux,uy,prop);
    printf("%u,%g,%g,%g,%g\n",t,prop[0],prop[1],prop[2],prop[3]);
}

void LBM::save_scalar(const char* name, mdspan<double, dextents<size_t, 2>> scalar, unsigned int n)
{
    // assume reasonably-sized file names
    char filename[128];
    char format[16];
    
    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS)+1.0);
    
    // generate format string
    // file name format is name0000nnn.bin
    sprintf(format,"%%s%%0%dd.bin",ndigits);
    sprintf(filename,format,name,n);
    
    // open file for writing
    FILE *fout = fopen(filename,"wb+");

    double tmp[scalar.size()];
    // write data
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            tmp[i+j]=scalar[i,j];
            //cout<<tmp[i+j]<<" ";
        }
        //std::cout<<endl;
    }

    fwrite(tmp,1,mem_size_scalar,fout);
    //TODO write data correctly
    // close file
    fclose(fout);
    
    if(ferror(fout))
    {
        fprintf(stderr,"Error saving to %s\n",filename);
        perror("");
    }
    else
    {
        if(!quiet)
            printf("Saved to %s\n",filename);
    }
}

