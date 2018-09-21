/**
 * Created on: Sep 21, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/
 
#include "cudaBCS.cuh"

__global__ void cacl_bssn_bcs_x(double * dev_var_out, double * dev_var_in, int u_offset, double * dxf, double * dyf, double * dzf, double pmin_x, double pmin_y, double pmin_z, double pmax_x, double pmax_y, double pmax_z, const double f_falloff, const double f_asymptotic, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) 
{
    int j = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int k = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int nx = host_sz_x;
    int ny = host_sz_y;
    int nz = host_sz_z;

    if(j >= ny-3 || k >= nz-3) return;

    double inv_r;
    double hx = (pmax_x - pmin_x) / (nx - 1);
    double hy = (pmax_y - pmin_y) / (ny - 1);
    double hz = (pmax_z - pmin_z) / (nz - 1);
    double x, y, z;
    int pp;

    if (bflag & (1u<<OCT_DIR_LEFT)) {
        
        x = pmin_x + 3*hx;
        z = pmin_z + k*hz;
        y = pmin_y + j*hy;
        pp = IDX(3,j,k);
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        dev_var_out[u_offset + pp] = -  inv_r * (
                        x * dxf[pp]
                        + y * dyf[pp]
                        + z * dzf[pp]
                        + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
        }
    
    if (bflag & (1u<<OCT_DIR_RIGHT)) {
        x = pmin_x + (nx - 3)*hx;
        z = pmin_z + k*hz;
        y = pmin_y + j*hy;
        pp = IDX((nx - 3),j,k);
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        dev_var_out[u_offset + pp] = -  inv_r * (
                        x * dxf[pp]
                    + y * dyf[pp]
                    + z * dzf[pp]
                    + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
    }
}

__global__ void cacl_bssn_bcs_y(double * dev_var_out, double * dev_var_in, int u_offset, double * dxf, double * dyf, double * dzf, double pmin_x, double pmin_y, double pmin_z, double pmax_x, double pmax_y, double pmax_z, const double f_falloff, const double f_asymptotic, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) 
{
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int k = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int nx = host_sz_x;
    int ny = host_sz_y;
    int nz = host_sz_z;

    if(i >= nx-3 || k >= nz-3) return;

    double inv_r;
    double hx = (pmax_x - pmin_x) / (nx - 1);
    double hy = (pmax_y - pmin_y) / (ny - 1);
    double hz = (pmax_z - pmin_z) / (nz - 1);
    double x, y, z;
    int pp;

    if (bflag & (1u<<OCT_DIR_DOWN)) {
        
        y = pmin_y + 3*hy;
        z = pmin_z + k*hz;
        x = pmin_x + i*hx;
        pp = IDX(i,3,k);
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        dev_var_out[u_offset + pp] = -  inv_r * (
                        x * dxf[pp]
                        + y * dyf[pp]
                        + z * dzf[pp]
                        + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
        
        }
    
    if (bflag & (1u<<OCT_DIR_UP)) {
        x = pmin_x + i*hx;
        z = pmin_z + k*hz;
        y = pmin_y + (ny-3)*hy;
        pp = IDX(i,(ny - 3),k);
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        dev_var_out[u_offset + pp] = -  inv_r * (
                        x * dxf[pp]
                    + y * dyf[pp]
                    + z * dzf[pp]
                    + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
        
    }
}

__global__ void cacl_bssn_bcs_z(double * dev_var_out, double * dev_var_in, int u_offset, double * dxf, double * dyf, double * dzf, double pmin_x, double pmin_y, double pmin_z, double pmax_x, double pmax_y, double pmax_z, const double f_falloff, const double f_asymptotic, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag) 
{
    int i = 3 + threadIdx.x + blockIdx.x * blockDim.x;
    int j = 3 + threadIdx.y + blockIdx.y * blockDim.y;
    int nx = host_sz_x;
    int ny = host_sz_y;
    int nz = host_sz_z;

    if(i >= nx-3 || j >= ny-3) return;

    double inv_r;
    double hx = (pmax_x - pmin_x) / (nx - 1);
    double hy = (pmax_y - pmin_y) / (ny - 1);
    double hz = (pmax_z - pmin_z) / (nz - 1);
    double x, y, z;
    int pp;

    if (bflag & (1u<<OCT_DIR_BACK)) {
        
        y = pmin_y + j*hy;
        z = pmin_z + 3*hz;
        x = pmin_x + i*hx;
        pp = IDX(i,j,3);
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        dev_var_out[u_offset + pp] = -  inv_r * (
                        x * dxf[pp]
                        + y * dyf[pp]
                        + z * dzf[pp]
                        + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
    
    }
    
    if (bflag & (1u<<OCT_DIR_FRONT)) {
        x = pmin_x + i*hx;
        z = pmin_z + (nz-3)*hz;
        y = pmin_y + j*hy;
        pp = IDX(i,j,3);
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        dev_var_out[u_offset + pp] = -  inv_r * (
                            x * dxf[pp]
                        + y * dyf[pp]
                        + z * dzf[pp]
                        + f_falloff * (   dev_var_in[u_offset + pp] - f_asymptotic ) );
        
    }
}


void bssn_bcs(double * dev_var_out, double * dev_var_in, int u_offset, double * dxf, double * dyf, double * dzf, const double * pmin, const double * pmax, const double f_falloff, const double f_asymptotic, const unsigned int * host_sz, int bflag, cudaStream_t stream) 
{
    const unsigned int nx = host_sz[0];
    const unsigned int ny = host_sz[1];
    const unsigned int nz = host_sz[2];

    const int ie = nx - 3;//x direction
    const int je = ny - 3;//y direction
    const int ke = nz - 3;//z direction

    double pmin_x = pmin[0];
    double pmin_y = pmin[1];
    double pmin_z = pmin[2];

    double pmax_x = pmax[0];
    double pmax_y = pmax[1];
    double pmax_z = pmax[2];

    const unsigned int host_sz_x = host_sz[0];
    const unsigned int host_sz_y = host_sz[1];
    const unsigned int host_sz_z = host_sz[2];

    int maximumIterations = (je>ke) ? je: ke;
    
    int requiredBlocks = (9 + maximumIterations) / 10;
    
    int threads_y = (requiredBlocks-1+je) / requiredBlocks;
    int threads_z = (requiredBlocks-1+ke) / requiredBlocks;
    
    cacl_bssn_bcs_x <<< dim3(threads_y,threads_z), dim3(threads_y,threads_z), 0, stream >>> (
        dev_var_out, dev_var_in, 
        u_offset, dxf, dyf, dzf, 
        pmin_x, pmin_y, pmin_z, pmax_x, pmax_y, pmax_z, 
        f_falloff, f_asymptotic, 
        host_sz_x, host_sz_y, host_sz_z, 
        bflag );
    
    CHECK_ERROR(cudaGetLastError(), "cacl_bssn_bcs_x Kernel launch failed");
        
    maximumIterations = (ke>ie) ? ke : ie ;
    requiredBlocks = (9 + maximumIterations)/10;
    int threads_x = (requiredBlocks-1+ie) / requiredBlocks;
    threads_z = (requiredBlocks-1+ke) / requiredBlocks;
    cacl_bssn_bcs_y <<< dim3(threads_x,threads_z), dim3(threads_x,threads_z), 0, stream >>> (
        dev_var_out, dev_var_in, 
        u_offset, dxf, dyf, dzf, 
        pmin_x, pmin_y, pmin_z, pmax_x, pmax_y, pmax_z, 
        f_falloff, f_asymptotic, 
        host_sz_x, host_sz_y, host_sz_z, 
        bflag );

    CHECK_ERROR(cudaGetLastError(), "cacl_bssn_bcs_y Kernel launch failed");

    maximumIterations = (je>ie) ? je : ie ;
    requiredBlocks = (9 + maximumIterations)/10;
    threads_x = (requiredBlocks-1+ie) / requiredBlocks;
    threads_y = (requiredBlocks-1+je) / requiredBlocks;
    cacl_bssn_bcs_z <<< dim3(threads_x,threads_y), dim3(threads_x,threads_y), 0, stream >>> (
        dev_var_out, dev_var_in, 
        u_offset, dxf, dyf, dzf, 
        pmin_x, pmin_y, pmin_z, pmax_x, pmax_y, pmax_z, 
        f_falloff, f_asymptotic, 
        host_sz_x, host_sz_y, host_sz_z, 
        bflag );

    CHECK_ERROR(cudaGetLastError(), "cacl_bssn_bcs_z Kernel launch failed");
}