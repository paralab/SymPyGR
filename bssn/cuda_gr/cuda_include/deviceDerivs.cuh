/**
 * Created on: March 15, 2018
 * 		Author: Akila, Eranga, Eminda, Ruwan
 **/

#ifndef DEVICE_DERIV_CUH
#define DEVICE_DERIV_CUH

#include "cuda_runtime.h"
#include "def.h"

#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

__device__ void calc_deriv42_x_bflag(int id, double * output, double * dev_var_in, const int u_offset, double dx, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_y_bflag(int id, double * output, double * dev_var_in, const int u_offset, double dy, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_z_bflag(int id, double * output, double * dev_var_in, const int u_offset, double dz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);

__device__ void calc_deriv42_xx_bflag(int id, double * output, double * dev_var_in, const int u_offset, double dx, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_yy_bflag(int id, double * output, double * dev_var_in, const int u_offset, double dy, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_zz_bflag(int id, double * output, double * dev_var_in, const int u_offset, double dz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);

__device__ void calc_deriv42_adv_x_bflag(int id, double * output, double * dev_var_in, int u_offset, double dx, int betax, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_adv_y_bflag(int id, double * output, double * dev_var_in, int u_offset, double dy, int betay, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_adv_z_bflag(int id, double * output, double * dev_var_in, int u_offset, double dz, int betaz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);

__device__ void calc_ko_deriv42_x_bflag(int id, double * output, double * dev_var_in, int u_offset, double dx, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_ko_deriv42_y_bflag(int id, double * output, double * dev_var_in, int u_offset, double dy, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_ko_deriv42_z_bflag(int id, double * output, double * dev_var_in, int u_offset, double dz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);


__device__ void calc_deriv42_x(int tile_size, double * output, double * dev_var_in, const int u_offset, double dx, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_y(int tile_size, double * output, double * dev_var_in, const int u_offset, double dy, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_z(int tile_size, double * output, double * dev_var_in, const int u_offset, double dz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);

__device__ void calc_deriv42_xx(int tile_size, double * output, double * dev_var_in, const int u_offset, double dx, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_yy(int tile_size, double * output, double * dev_var_in, const int u_offset, double dy, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_zz(int tile_size, double * output, double * dev_var_in, const int u_offset, double dz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);

__device__ void calc_deriv42_adv_x(int tile_size, double * output, double * dev_global_var_in, double * dev_var_in, int u_offset, double dx, int betax, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_adv_y(int tile_size, double * output, double * dev_global_var_in, double * dev_var_in, int u_offset, double dy, int betay, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_deriv42_adv_z(int tile_size, double * output, double * dev_global_var_in, double * dev_var_in, int u_offset, double dz, int betaz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);

__device__ void calc_ko_deriv42_x(int id, double * output, double * dev_var_in, int u_offset, double dx, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_ko_deriv42_y(int id, double * output, double * dev_var_in, int u_offset, double dy, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);
__device__ void calc_ko_deriv42_z(int id, double * output, double * dev_var_in, int u_offset, double dz, const unsigned int host_sz_x, const unsigned int host_sz_y, const unsigned int host_sz_z, int bflag);

__device__ void globaltoshared(double * dev_var_in, double * shared_var_in, int glb_offset, int halo_length, int x_offset, int y_offset, int z_offset, int nx, int ny);

#endif
