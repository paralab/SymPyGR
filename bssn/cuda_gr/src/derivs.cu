//
// Created by milinda on 8/9/18.
//

/**
 * @brief Contians cuda derivs for bssn computation
 *
 * */

#include "derivs.cuh"
#include "dendro.h"

namespace cuda
{



/*----------------------------------------------------------------------;
 *
 * compute first derivative in x direction
 *
 *----------------------------------------------------------------------*/

        __device__ void deriv42_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag) {

            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)1);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-1);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(1));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-1);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];


            const double idx = 1.0 / dx;
            const double idx_by_2 = 0.5 * idx;
            const double idx_by_12 = idx / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 1;
            const int kb = 1;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 1;
            const int ke = sz[2] - 1;

            //printf("dx threadid (%d,%d,%d) loop begin: (%d,%d,%d) loop end: (%d,%d,%d)  tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \n", threadIdx.x,threadIdx.y,threadIdx.z,ix_b,jy_b,kz_b,ix_e,jy_e,kz_e,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);
                        Dxu[pp] = (u[pp - 2] - 8.0 * u[pp - 1] + 8.0 * u[pp + 1] - u[pp + 2]) * idx_by_12;

                    }



        if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && ix_b==ib)  ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        Dxu[IDX(3, j, k)] = (-3.0 * u[IDX(3, j, k)]
                                     + 4.0 * u[IDX(4, j, k)]
                                     - u[IDX(5, j, k)]
                                    ) * idx_by_2;

                        Dxu[IDX(4, j, k)] = (-u[IDX(3, j, k)]
                                     + u[IDX(5, j, k)]
                                    ) * idx_by_2;

                    }

                
            }

        if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && ix_e==ie) ) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        Dxu[IDX(ie - 2, j, k)] = (-u[IDX(ie - 3, j, k)]
                                          + u[IDX(ie - 1, j, k)]
                                         ) * idx_by_2;

                        Dxu[IDX(ie - 1, j, k)] = (u[IDX(ie - 3, j, k)]
                                          - 4.0 * u[IDX(ie - 2, j, k)]
                                          + 3.0 * u[IDX(ie - 1, j, k)]
                                         ) * idx_by_2;
                    }

                

            }


#ifdef DEBUG_DERIVS_COMP
            if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }

/*----------------------------------------------------------------------;
 *
 * compute first derivative in y direction
 *
 *----------------------------------------------------------------------*/

        __device__ void deriv42_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(1));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-1);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idy = 1.0 / dy;
            const double idy_by_2 = 0.5 * idy;
            const double idy_by_12 = idy / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 1;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 1;

            //printf("dy threadid (%d,%d,%d) loop begin: (%d,%d,%d) loop end: (%d,%d,%d)  tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \n", threadIdx.x,threadIdx.y,threadIdx.z,ix_b,jy_b,kz_b,ix_e,jy_e,kz_e,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);
                        Dyu[pp] = (u[pp - 2 * nx] - 8.0 * u[pp - nx] + 8.0 * u[pp + nx] - u[pp + 2 * nx]) * idy_by_12;
                    }




        if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && jy_b==jb) ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Dyu[IDX(i, 3, k)] = (-3.0 * u[IDX(i, 3, k)]
                                     + 4.0 * u[IDX(i, 4, k)]
                                     - u[IDX(i, 5, k)]
                                    ) * idy_by_2;

                        Dyu[IDX(i, 4, k)] = (-u[IDX(i, 3, k)]
                                     + u[IDX(i, 5, k)]
                                    ) * idy_by_2;


                    }

                
            }

        if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && jy_e==je)) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        Dyu[IDX(i, je - 2, k)] = (-u[IDX(i, je - 3, k)]
                                          + u[IDX(i, je - 1, k)]
                                         ) * idy_by_2;

                        Dyu[IDX(i, je - 1, k)] = (u[IDX(i, je - 3, k)]
                                          - 4.0 * u[IDX(i, je - 2, k)]
                                          + 3.0 * u[IDX(i, je - 1, k)]
                                         ) * idy_by_2;

                    }

                

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif

        }

/*----------------------------------------------------------------------;
 *
 * compute first derivative in z direction
 *
 *----------------------------------------------------------------------*/


        __device__ void deriv42_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idz = 1.0 / dz;
            const double idz_by_2 = 0.5 * idz;
            const double idz_by_12 = idz / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;

            const int n = nx * ny;

            //printf("dz threadid (%d,%d,%d) loop begin: (%d,%d,%d) loop end: (%d,%d,%d)  tile begin: (%d,%d,%d) tile end: (%d,%d,%d) \n", threadIdx.x,threadIdx.y,threadIdx.z,ix_b,jy_b,kz_b,ix_e,jy_e,kz_e,ijk_lm[0],ijk_lm[2],ijk_lm[4],ijk_lm[1],ijk_lm[3],ijk_lm[5]);

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);
                        Dzu[pp] = (u[pp - 2 * n] - 8.0 * u[pp - n] + 8.0 * u[pp + n] - u[pp + 2 * n]) * idz_by_12;
                    }


        if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && kz_b==kb) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Dzu[IDX(i, j, 3)] = (-3.0 * u[IDX(i, j, 3)]
                                     + 4.0 * u[IDX(i, j, 4)]
                                     - u[IDX(i, j, 5)]
                                    ) * idz_by_2;

                       Dzu[IDX(i, j, 4)] = (-u[IDX(i, j, 3)]
                                     + u[IDX(i, j, 5)]
                                    ) * idz_by_2;

                    }

            }

        if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && kz_e==ke) ) {


                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Dzu[IDX(i, j, ke - 2)] = (-u[IDX(i, j, ke - 3)]
                                          + u[IDX(i, j, ke - 1)]
                                         ) * idz_by_2;

                        Dzu[IDX(i, j, ke - 1)] = (u[IDX(i, j, ke - 3)]
                                          - 4.0 * u[IDX(i, j, ke - 2)]
                                          + 3.0 * u[IDX(i, j, ke - 1)]
                                         ) * idz_by_2;
                    }

                

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute second derivative in x direction
 *
 *----------------------------------------------------------------------*/

        __device__ void deriv42_xx(double * const  DxDxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw,unsigned bflag) {


            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idx_sqrd = 1.0 / (dx * dx);
            const double idx_sqrd_by_12 = idx_sqrd / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;


            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        const int pp = IDX(i, j, k);
                        DxDxu[pp] = (-u[pp - 2]
                         + 16.0 * u[pp - 1]
                         - 30.0 * u[pp]
                         + 16.0 * u[pp + 1]
                         - u[pp + 2]
                        ) * idx_sqrd_by_12;

                    }




        if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && ix_b==ib)  ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {

                        DxDxu[IDX(3, j, k)] = (2.0 * u[IDX(3, j, k)]
                                       - 5.0 * u[IDX(4, j, k)]
                                       + 4.0 * u[IDX(5, j, k)]
                                       - u[IDX(6, j, k)]
                                      ) * idx_sqrd;

                        DxDxu[IDX(4, j, k)] = (u[IDX(3, j, k)]
                                       - 2.0 * u[IDX(4, j, k)]
                                       + u[IDX(5, j, k)]
                                      ) * idx_sqrd;


                    }


                

            }

        if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && ix_e==ie) ) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        DxDxu[IDX(ie - 2, j, k)] = (u[IDX(ie - 3, j, k)]
                                            - 2.0 * u[IDX(ie - 2, j, k)]
                                            + u[IDX(ie - 1, j, k)]
                                           ) * idx_sqrd;

                        DxDxu[IDX(ie - 1, j, k)] = (-u[IDX(ie - 4, j, k)]
                                            + 4.0 * u[IDX(ie - 3, j, k)]
                                            - 5.0 * u[IDX(ie - 2, j, k)]
                                            + 2.0 * u[IDX(ie - 1, j, k)]
                                           ) * idx_sqrd;
                    }



                
            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif

        }


/*----------------------------------------------------------------------;
 *
 * compute second derivative in y direction
 *
 *----------------------------------------------------------------------*/



        __device__ void deriv42_yy(double * const  DyDyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {


            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idy_sqrd = 1.0 / (dy * dy);
            const double idy_sqrd_by_12 = idy_sqrd / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];


            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);
                        DyDyu[pp] = (-u[pp - 2 * nx] + 16.0 * u[pp - nx] - 30.0 * u[pp]
                            + 16.0 * u[pp + nx] - u[pp + 2 * nx]
                            ) * idy_sqrd_by_12;

                    }





        if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && jy_b==jb) ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        DyDyu[IDX(i, 3, k)] = (2.0 * u[IDX(i, 3, k)]
                                       - 5.0 * u[IDX(i, 4, k)]
                                       + 4.0 * u[IDX(i, 5, k)]
                                       - u[IDX(i, 6, k)]
                                      ) * idy_sqrd;

                        DyDyu[IDX(i, 4, k)] = (u[IDX(i, 3, k)]
                                       - 2.0 * u[IDX(i, 4, k)]
                                       + u[IDX(i, 5, k)]
                                      ) * idy_sqrd;

                    }

                
            }

        if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && jy_e==je)) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        DyDyu[IDX(i, je - 2, k)] = (u[IDX(i, je - 3, k)]
                                            - 2.0 * u[IDX(i, je - 2, k)]
                                            + u[IDX(i, je - 1, k)]
                                           ) * idy_sqrd;

                        DyDyu[IDX(i, je - 1, k)] = (-u[IDX(i, je - 4, k)]
                                            + 4.0 * u[IDX(i, je - 3, k)]
                                            - 5.0 * u[IDX(i, je - 2, k)]
                                            + 2.0 * u[IDX(i, je - 1, k)]
                                           ) * idy_sqrd;

                    }

                

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }



        /*----------------------------------------------------------------------;
        *
        * compute second derivative in z direction
        *
        *----------------------------------------------------------------------*/


        __device__ void deriv42_zz(double * const  DzDzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idz_sqrd = 1.0 / (dz * dz);
            const double idz_sqrd_by_12 = idz_sqrd / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];


            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;
            const int n = nx * ny;

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);
                        DzDzu[pp] = (-u[pp - 2 * n] + 16.0 * u[pp - n] - 30.0 * u[pp]
                             + 16.0 * u[pp + n] - u[pp + 2 * n]) * idz_sqrd_by_12;
                    }



        if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && kz_b==kb) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        DzDzu[IDX(i, j, 3)] = (2.0 * u[IDX(i, j, 3)]
                                       - 5.0 * u[IDX(i, j, 4)]
                                       + 4.0 * u[IDX(i, j, 5)]
                                       - u[IDX(i, j, 6)]
                                      ) * idz_sqrd;

                        DzDzu[IDX(i, j, 4)] = (u[IDX(i, j, 3)]
                                       - 2.0 * u[IDX(i, j, 4)]
                                       + u[IDX(i, j, 5)]
                                      ) * idz_sqrd;

                    }

                
            }

        if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && kz_e==ke) ) {


                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        DzDzu[IDX(i, j, ke - 2)] = (u[IDX(i, j, ke - 3)]
                                            - 2.0 * u[IDX(i, j, ke - 2)]
                                            + u[IDX(i, j, ke - 1)]
                                           ) * idz_sqrd;

                        DzDzu[IDX(i, j, ke - 1)] = (-u[IDX(i, j, ke - 4)]
                                            + 4.0 * u[IDX(i, j, ke - 3)]
                                            - 5.0 * u[IDX(i, j, ke - 2)]
                                            + 2.0 * u[IDX(i, j, ke - 1)]
                                           ) * idz_sqrd;

                    }

                
            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute first advective derivative in x direction
 *
 *----------------------------------------------------------------------*/
        __device__    void deriv42adv_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const double * const betax, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idx = 1.0 / dx;
            const double idx_by_2 = 0.50 * idx;
            const double idx_by_12 = idx / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);


                        if (betax[pp] > 0.0) {
                            Dxu[pp] = (-3.0 * u[pp - 1]
                                    - 10.0 * u[pp]
                                    + 18.0 * u[pp + 1]
                                    - 6.0 * u[pp + 2]
                                    + u[pp + 3]
                                    ) * idx_by_12;
                        } else {
                            Dxu[pp] = (-u[pp - 3]
                                    + 6.0 * u[pp - 2]
                                    - 18.0 * u[pp - 1]
                                    + 10.0 * u[pp]
                                    + 3.0 * u[pp + 1]
                                    ) * idx_by_12;
                        }

                    }



        if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && ix_b==ib)  ) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        Dxu[IDX(3, j, k)] = (-3.0 * u[IDX(3, j, k)]
                                     + 4.0 * u[IDX(4, j, k)]
                                     - u[IDX(5, j, k)]
                                    ) * idx_by_2;

                        if (betax[IDX(4, j, k)] > 0.0) {
                            Dxu[IDX(4, j, k)] = (-3.0 * u[IDX(4, j, k)]
                                                + 4.0 * u[IDX(5, j, k)]
                                                - u[IDX(6, j, k)]
                                                ) * idx_by_2;
                        } else {
                            Dxu[IDX(4, j, k)] = (-u[IDX(3, j, k)]
                                                + u[IDX(5, j, k)]
                                                ) * idx_by_2;
                        }

                        if (betax[IDX(5, j, k)] > 0.0) {
                            Dxu[IDX(5, j, k)] = (-3.0 * u[IDX(4, j, k)]
                                                - 10.0 * u[IDX(5, j, k)]
                                                + 18.0 * u[IDX(6, j, k)]
                                                - 6.0 * u[IDX(7, j, k)]
                                                + u[IDX(8, j, k)]
                                                ) * idx_by_12;
                        } else {
                            Dxu[IDX(5, j, k)] = (u[IDX(3, j, k)]
                                                - 4.0 * u[IDX(4, j, k)]
                                                + 3.0 * u[IDX(5, j, k)]
                                                ) * idx_by_2;
                        }

                    }

                


            }

        if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && ix_e==ie) ) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        if (betax[IDX(ie - 3, j, k)] < 0.0) {
                            Dxu[IDX(ie - 3, j, k)] = (-3.0 * u[IDX(ie - 3, j, k)]
                                                      + 4.0 * u[IDX(ie - 2, j, k)]
                                                      - u[IDX(ie - 1, j, k)]
                                                     ) * idx_by_2;
                        } else {
                            Dxu[IDX(ie - 3, j, k)] = (-u[IDX(ie - 6, j, k)]
                                                      + 6.0 * u[IDX(ie - 5, j, k)]
                                                      - 18.0 * u[IDX(ie - 4, j, k)]
                                                      + 10.0 * u[IDX(ie - 3, j, k)]
                                                      + 3.0 * u[IDX(ie - 2, j, k)]
                                                     ) * idx_by_12;
                        }
        
                        if (betax[IDX(ie - 2, j, k)] > 0.0) {
                            Dxu[IDX(ie - 2, j, k)] = (-u[IDX(ie - 3, j, k)]
                                                      + u[IDX(ie - 1, j, k)]
                                                     ) * idx_by_2;
                        } else {
                            Dxu[IDX(ie - 2, j, k)] = (u[IDX(ie - 4, j, k)]
                                                      - 4.0 * u[IDX(ie - 3, j, k)]
                                                      + 3.0 * u[IDX(ie - 2, j, k)]
                                                     ) * idx_by_2;
                        }
        
                        Dxu[IDX(ie - 1, j, k)] = (u[IDX(ie - 3, j, k)]
                                                  - 4.0 * u[IDX(ie - 2, j, k)]
                                                  + 3.0 * u[IDX(ie - 1, j, k)]
                                                 ) * idx_by_2;

                    }


                

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute first advective derivative in y direction
 *
 *----------------------------------------------------------------------*/
        __device__  void deriv42adv_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const double * const betay, unsigned int pw, unsigned bflag) {

            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idy = 1.0 / dy;
            const double idy_by_2 = 0.50 * idy;
            const double idy_by_12 = idy / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);


                        if (betay[pp] > 0.0) {
                            Dyu[pp] = (-3.0 * u[pp - nx]
                                    - 10.0 * u[pp]
                                    + 18.0 * u[pp + nx]
                                    - 6.0 * u[pp + 2 * nx]
                                    + u[pp + 3 * nx]
                                    ) * idy_by_12;
                        } else {
                            Dyu[pp] = (-u[pp - 3 * nx]
                                    + 6.0 * u[pp - 2 * nx]
                                    - 18.0 * u[pp - nx]
                                    + 10.0 * u[pp]
                                    + 3.0 * u[pp + nx]
                                    ) * idy_by_12;
                        }

                    }







        if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && jy_b==jb) ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Dyu[IDX(i, 3, k)] = (-3.0 * u[IDX(i, 3, k)]
                                     + 4.0 * u[IDX(i, 4, k)]
                                     - u[IDX(i, 5, k)]
                                    ) * idy_by_2;

                        if (betay[IDX(i, 4, k)] > 0.0) {
                            Dyu[IDX(i, 4, k)] = (-3.0 * u[IDX(i, 4, k)]
                                                + 4.0 * u[IDX(i, 5, k)]
                                                - u[IDX(i, 6, k)]
                                                ) * idy_by_2;
                        } else {
                            Dyu[IDX(i, 4, k)] = (-u[IDX(i, 3, k)]
                                                + u[IDX(i, 5, k)]
                                                ) * idy_by_2;
                        }

                        if (betay[IDX(i, 5, k)] > 0.0) {
                            Dyu[IDX(i, 5, k)] = (-3.0 * u[IDX(i, 4, k)]
                                                - 10.0 * u[IDX(i, 5, k)]
                                                + 18.0 * u[IDX(i, 6, k)]
                                                - 6.0 * u[IDX(i, 7, k)]
                                                + u[IDX(i, 8, k)]
                                                ) * idy_by_12;
                        } else {
                            Dyu[IDX(i, 5, k)] = (u[IDX(i, 3, k)]
                                                - 4.0 * u[IDX(i, 4, k)]
                                                + 3.0 * u[IDX(i, 5, k)]
                                                ) * idy_by_2;
                        }

                    }

                

            }

        if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && jy_e==je)) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        if (betay[IDX(i, je - 3, k)] < 0.0) {
                            Dyu[IDX(i, je - 3, k)] = (-3.0 * u[IDX(i, je - 3, k)]
                                                      + 4.0 * u[IDX(i, je - 2, k)]
                                                      - u[IDX(i, je - 1, k)]
                                                     ) * idy_by_2;
                        } else {
                            Dyu[IDX(i, je - 3, k)] = (-u[IDX(i, je - 6, k)]
                                                      + 6.0 * u[IDX(i, je - 5, k)]
                                                      - 18.0 * u[IDX(i, je - 4, k)]
                                                      + 10.0 * u[IDX(i, je - 3, k)]
                                                      + 3.0 * u[IDX(i, je - 2, k)]
                                                     ) * idy_by_12;
                        }
        
                        if (betay[IDX(i, je - 2, k)] > 0.0) {
                            Dyu[IDX(i, je - 2, k)] = (-u[IDX(i, je - 3, k)]
                                                      + u[IDX(i, je - 1, k)]
                                                     ) * idy_by_2;
                        } else {
                            Dyu[IDX(i, je - 2, k)] = (u[IDX(i, je - 4, k)]
                                                      - 4.0 * u[IDX(i, je - 3, k)]
                                                      + 3.0 * u[IDX(i, je - 2, k)]
                                                     ) * idy_by_2;
                        }
        
                        Dyu[IDX(i, je - 1, k)] = (u[IDX(i, je - 3, k)]
                                                  - 4.0 * u[IDX(i, je - 2, k)]
                                                  + 3.0 * u[IDX(i, je - 1, k)]
                                                 ) * idy_by_2;

                    }

                


            }


#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }


/*----------------------------------------------------------------------;
 *
 * compute first advective derivative in z direction
 *
 *----------------------------------------------------------------------*/


        __device__  void deriv42adv_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const double * const betaz, unsigned int pw, unsigned bflag) {


            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            const double idz = 1.0 / dz;
            const double idz_by_2 = 0.50 * idz;
            const double idz_by_12 = idz / 12.0;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;

            const int n = nx * ny;

            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);


                        if (betaz[pp] > 0.0) {
                            Dzu[pp] = (-3.0 * u[pp - n]
                                    - 10.0 * u[pp]
                                    + 18.0 * u[pp + n]
                                    - 6.0 * u[pp + 2 * n]
                                    + u[pp + 3 * n]
                                    ) * idz_by_12;
                        } else {
                            Dzu[pp] = (-u[pp - 3 * n]
                                    + 6.0 * u[pp - 2 * n]
                                    - 18.0 * u[pp - n]
                                    + 10.0 * u[pp]
                                    + 3.0 * u[pp + n]
                                    ) * idz_by_12;
                        }


                    }





        if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && kz_b==kb) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        Dzu[IDX(i, j, 3)] = (-3.0 * u[IDX(i, j, 3)]
                                     + 4.0 * u[IDX(i, j, 4)]
                                     - u[IDX(i, j, 5)]
                                    ) * idz_by_2;

                        if (betaz[IDX(i, j, 4)] > 0.0) {
                            Dzu[IDX(i, j, 4)] = (-3.0 * u[IDX(i, j, 4)]
                                                + 4.0 * u[IDX(i, j, 5)]
                                                - u[IDX(i, j, 6)]
                                                ) * idz_by_2;
                        } else {
                            Dzu[IDX(i, j, 4)] = (-u[IDX(i, j, 3)]
                                                + u[IDX(i, j, 5)]
                                                ) * idz_by_2;
                        }

                        if (betaz[IDX(i, j, 5)] > 0.0) {
                            Dzu[IDX(i, j, 5)] = (-3.0 * u[IDX(i, j, 4)]
                                                - 10.0 * u[IDX(i, j, 5)]
                                                + 18.0 * u[IDX(i, j, 6)]
                                                - 6.0 * u[IDX(i, j, 7)]
                                                + u[IDX(i, j, 8)]
                                                ) * idz_by_12;
                        } else {
                            Dzu[IDX(i, j, 5)] = (u[IDX(i, j, 3)]
                                                - 4.0 * u[IDX(i, j, 4)]
                                                + 3.0 * u[IDX(i, j, 5)]
                                                ) * idz_by_2;
                        }

                    }

                

            }

        if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && kz_e==ke) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        if (betaz[IDX(i, j, ke - 3)] < 0.0) {
                            Dzu[IDX(i, j, ke - 3)] = (-3.0 * u[IDX(i, j, ke - 3)]
                                                      + 4.0 * u[IDX(i, j, ke - 2)]
                                                      - u[IDX(i, j, ke - 1)]
                                                     ) * idz_by_2;
                        } else {
                            Dzu[IDX(i, j, ke - 3)] = (-u[IDX(i, j, ke - 6)]
                                                      + 6.0 * u[IDX(i, j, ke - 5)]
                                                      - 18.0 * u[IDX(i, j, ke - 4)]
                                                      + 10.0 * u[IDX(i, j, ke - 3)]
                                                      + 3.0 * u[IDX(i, j, ke - 2)]
                                                     ) * idz_by_12;
                        }
        
                        if (betaz[IDX(i, j, ke - 2)] > 0.0) {
                            Dzu[IDX(i, j, ke - 2)] = (-u[IDX(i, j, ke - 3)]
                                                      + u[IDX(i, j, ke - 1)]
                                                     ) * idz_by_2;
                        } else {
                            Dzu[IDX(i, j, ke - 2)] = (u[IDX(i, j, ke - 4)]
                                                      - 4.0 * u[IDX(i, j, ke - 3)]
                                                      + 3.0 * u[IDX(i, j, ke - 2)]
                                                     ) * idz_by_2;
                        }
        
                        Dzu[IDX(i, j, ke - 1)] = (u[IDX(i, j, ke - 3)]
                                                  - 4.0 * u[IDX(i, j, ke - 2)]
                                                  + 3.0 * u[IDX(i, j, ke - 1)]
                                                 ) * idz_by_2;

                    }

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }

/*----------------------------------------------------------------------
 *
 * compute Kriess-Oliger derivative in x direction
 *
 *----------------------------------------------------------------------*/


        __device__  void ko_deriv42_x(double * const  Du, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {


            const unsigned int i_b=ijk_lm[2*0+0]+pw;
            const unsigned int i_e=ijk_lm[2*0+1]-pw;

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;


            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            double pre_factor_6_dx = -1.0 / 64.0 / dx;

            double smr3 = 59.0 / 48.0 * 64 * dx;
            double smr2 = 43.0 / 48.0 * 64 * dx;
            double smr1 = 49.0 / 48.0 * 64 * dx;
            double spr3 = smr3;
            double spr2 = smr2;
            double spr1 = smr1;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;


            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);

                        Du[pp] = pre_factor_6_dx *
                                 (
                                         -u[pp - 3]
                                         + 6.0 * u[pp - 2]
                                         - 15.0 * u[pp - 1]
                                         + 20.0 * u[pp]
                                         - 15.0 * u[pp + 1]
                                         + 6.0 * u[pp + 2]
                                         - u[pp + 3]
                                 );
                    }



            if(i_b==ib && ix_b==ib)
            {
                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {

                        Du[IDX(ix_b, j, k)] = pre_factor_6_dx *
                                              (
                                                      -u[IDX(ix_b + 4, j, k)]
                                                      + 6.0 * u[IDX(ix_b + 3, j, k)]
                                                      - 15.0 * u[IDX(ix_b + 2, j, k)]
                                                      + 20.0 * u[IDX(ix_b + 1, j, k)]
                                                      - 15.0 * u[IDX(ix_b, j, k)]
                                                      + 6.0 * u[IDX(ix_b - 1, j, k)]
                                                      - u[IDX(ix_b - 2, j, k)]
                                              );

                    }

            }


            if(i_e==ie && ix_e==ie)
            {
                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {

                        Du[IDX(ix_e - 1, j, k)] = pre_factor_6_dx *
                                                  (
                                                          -u[IDX(ix_e + 1, j, k)]
                                                          + 6.0 * u[IDX(ix_e, j, k)]
                                                          - 15.0 * u[IDX(ix_e - 1, j, k)]
                                                          + 20.0 * u[IDX(ix_e - 2, j, k)]
                                                          - 15.0 * u[IDX(ix_e - 3, j, k)]
                                                          + 6.0 * u[IDX(ix_e - 4, j, k)]
                                                          - u[IDX(ix_e - 5, j, k)]
                                                  );
                    }
            }


            if ((bflag & (1u << OCT_DIR_LEFT)) && (i_b==ib && ix_b==ib)) {


                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        Du[IDX(3, j, k)] = (u[IDX(6, j, k)]
                                    - 3.0 * u[IDX(5, j, k)]
                                    + 3.0 * u[IDX(4, j, k)]
                                    - u[IDX(3, j, k)]
                                   ) / smr3;

                        Du[IDX(4, j, k)] = (
                                                u[IDX(7, j, k)]
                                                - 6.0 * u[IDX(6, j, k)]
                                                + 12.0 * u[IDX(5, j, k)]
                                                - 10.0 * u[IDX(4, j, k)]
                                                + 3.0 * u[IDX(3, j, k)]
                                        ) / smr2;
                        Du[IDX(5, j, k)] = (
                                                u[IDX(8, j, k)]
                                                - 6.0 * u[IDX(7, j, k)]
                                                + 15.0 * u[IDX(6, j, k)]
                                                - 19.0 * u[IDX(5, j, k)]
                                                + 12.0 * u[IDX(4, j, k)]
                                                - 3.0 * u[IDX(3, j, k)]
                                        ) / smr1;
                            }


                

            }

            if ((bflag & (1u << OCT_DIR_RIGHT)) && (i_e==ie && ix_e==ie)) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int j=jy_b;j<jy_e;j++)
                    {
                        Du[IDX(ie - 3, j, k)] = (
                            u[IDX(ie - 6, j, k)]
                            - 6.0 * u[IDX(ie - 5, j, k)]
                            + 15.0 * u[IDX(ie - 4, j, k)]
                            - 19.0 * u[IDX(ie - 3, j, k)]
                            + 12.0 * u[IDX(ie - 2, j, k)]
                            - 3.0 * u[IDX(ie - 1, j, k)]
                    ) / spr1;

                    Du[IDX(ie - 2, j, k)] = (
                                                u[IDX(ie - 5, j, k)]
                                                - 6.0 * u[IDX(ie - 4, j, k)]
                                                + 12.0 * u[IDX(ie - 3, j, k)]
                                                - 10.0 * u[IDX(ie - 2, j, k)]
                                                + 3.0 * u[IDX(ie - 1, j, k)]
                                        ) / spr2;

                    Du[IDX(ie - 1, j, k)] = (
                                                u[IDX(ie - 4, j, k)]
                                                - 3.0 * u[IDX(ie - 3, j, k)]
                                                + 3.0 * u[IDX(ie - 2, j, k)]
                                                - u[IDX(ie - 1, j, k)]
                                        ) / spr3;
                    }

                

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }



/*----------------------------------------------------------------------
 *
 * compute Kriess-Oliger derivative in y direction
 *
 *----------------------------------------------------------------------*/

        __device__  void ko_deriv42_y(double * const  Du, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {



            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=ijk_lm[2*1+0]+pw;
            const unsigned int j_e=ijk_lm[2*1+1]-pw;

            const unsigned int k_b=max((int)ijk_lm[2*2+0],int(3));
            const unsigned int k_e=min((int)ijk_lm[2*2+1],sz[2]-3);

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            double pre_factor_6_dy = -1.0 / 64.0 / dy;

            double smr3 = 59.0 / 48.0 * 64 * dy;
            double smr2 = 43.0 / 48.0 * 64 * dy;
            double smr1 = 49.0 / 48.0 * 64 * dy;
            double spr3 = smr3;
            double spr2 = smr2;
            double spr1 = smr1;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;



            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);

                        Du[pp] = pre_factor_6_dy *
                                (
                                        -u[pp - 3 * nx]
                                        + 6.0 * u[pp - 2 * nx]
                                        - 15.0 * u[pp - nx]
                                        + 20.0 * u[pp]
                                        - 15.0 * u[pp + nx]
                                        + 6.0 * u[pp + 2 * nx]
                                        - u[pp + 3 * nx]
                                );

                    }


            if(j_b==jb && jy_b==jb)
            {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Du[IDX(i, jy_b, k)] = pre_factor_6_dy *
                                              (
                                                      -u[IDX(i, jy_b + 4, k)]
                                                      + 6.0 * u[IDX(i, jy_b + 3, k)]
                                                      - 15.0 * u[IDX(i, jy_b + 2, k)]
                                                      + 20.0 * u[IDX(i, jy_b + 1, k)]
                                                      - 15.0 * u[IDX(i, jy_b, k)]
                                                      + 6.0 * u[IDX(i, jy_b - 1, k)]
                                                      - u[IDX(i, jy_b - 2, k)]
                                              );

                    }

            }


            if(j_e==je && jy_e==je)
            {
                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        Du[IDX(i, jy_e - 1, k)] = pre_factor_6_dy *
                                                  (
                                                          -u[IDX(i, jy_e + 1, k)]
                                                          + 6.0 * u[IDX(i, jy_e, k)]
                                                          - 15.0 * u[IDX(i, jy_e - 1, k)]
                                                          + 20.0 * u[IDX(i, jy_e - 2, k)]
                                                          - 15.0 * u[IDX(i, jy_e - 3, k)]
                                                          + 6.0 * u[IDX(i, jy_e - 4, k)]
                                                          - u[IDX(i, jy_e - 5, k)]
                                                  );
                    }

            }


            if ((bflag & (1u << OCT_DIR_DOWN)) && (j_b==jb && jy_b==jb) ) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Du[IDX(i, 3, k)] = (u[IDX(i, 6, k)]
                        - 3.0 * u[IDX(i, 5, k)]
                        + 3.0 * u[IDX(i, 4, k)]
                        - u[IDX(i, 3, k)]
                       ) / smr3;

                       Du[IDX(i, 4, k)] = (
                                                u[IDX(i, 7, k)]
                                                - 6.0 * u[IDX(i, 6, k)]
                                                + 12.0 * u[IDX(i, 5, k)]
                                                - 10.0 * u[IDX(i, 4, k)]
                                                + 3.0 * u[IDX(i, 3, k)]
                                        ) / smr2;
                       Du[IDX(i, 5, k)] = (
                                                u[IDX(i, 8, k)]
                                                - 6.0 * u[IDX(i, 7, k)]
                                                + 15.0 * u[IDX(i, 6, k)]
                                                - 19.0 * u[IDX(i, 5, k)]
                                                + 12.0 * u[IDX(i, 4, k)]
                                                - 3.0 * u[IDX(i, 3, k)]
                                        ) / smr1;

                    }


                
            }

            if ((bflag & (1u << OCT_DIR_UP)) && (j_e==je && jy_e==je)) {

                for(unsigned int k=kz_b;k<kz_e;k++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Du[IDX(i, je - 3, k)] = (
                            u[IDX(i, je - 6, k)]
                            - 6.0 * u[IDX(i, je - 5, k)]
                            + 15.0 * u[IDX(i, je - 4, k)]
                            - 19.0 * u[IDX(i, je - 3, k)]
                            + 12.0 * u[IDX(i, je - 2, k)]
                            - 3.0 * u[IDX(i, je - 1, k)]
                    ) / spr1;

                    Du[IDX(i, je - 2, k)] = (
                                                u[IDX(i, je - 5, k)]
                                                - 6.0 * u[IDX(i, je - 4, k)]
                                                + 12.0 * u[IDX(i, je - 3, k)]
                                                - 10.0 * u[IDX(i, je - 2, k)]
                                                + 3.0 * u[IDX(i, je - 1, k)]
                                        ) / spr2;

                    Du[IDX(i, je - 1, k)] = (
                                                u[IDX(i, je - 4, k)]
                                                - 3.0 * u[IDX(i, je - 3, k)]
                                                + 3.0 * u[IDX(i, je - 2, k)]
                                                - u[IDX(i, je - 1, k)]
                                        ) / spr3;

                    }

                

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }



/*----------------------------------------------------------------------
 *
 * compute Kriess-Oliger derivative in z direction
 *
 *----------------------------------------------------------------------*/

        __device__  void ko_deriv42_z(double * const  Du, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag) {



            const unsigned int i_b=max((int)ijk_lm[2*0+0],(int)3);
            const unsigned int i_e=min((int)ijk_lm[2*0+1],sz[0]-3);

            const unsigned int j_b=max((int)ijk_lm[2*1+0],(int)3);
            const unsigned int j_e=min((int)ijk_lm[2*1+1],sz[1]-3);

            const unsigned int k_b=ijk_lm[2*2+0]+pw;
            const unsigned int k_e=ijk_lm[2*2+1]-pw;

            unsigned int l_x=i_e-i_b;
            unsigned int l_y=j_e-j_b;
            unsigned int l_z=k_e-k_b;

            if(threadIdx.x>=l_x || threadIdx.y >= l_y || threadIdx.z>=l_z) return;

            if(l_x<blockDim.x) l_x=blockDim.x;
            if(l_y<blockDim.y) l_y=blockDim.y;
            if(l_z<blockDim.z) l_z=blockDim.z;

            const unsigned int ix_b= (i_b + (threadIdx.x * l_x)/blockDim.x)-ijk_lm[0];
            const unsigned int ix_e= (i_b + ((threadIdx.x+1) * l_x)/blockDim.x)-ijk_lm[0];

            const unsigned int jy_b= (j_b + (threadIdx.y * l_y)/blockDim.y)-ijk_lm[2];
            const unsigned int jy_e= (j_b + ((threadIdx.y+1) * l_y)/blockDim.y)-ijk_lm[2];

            const unsigned int kz_b= (k_b + (threadIdx.z * (l_z))/blockDim.z)-ijk_lm[4];
            const unsigned int kz_e= (k_b + ((threadIdx.z+1) * (l_z))/blockDim.z)-ijk_lm[4];

            double pre_factor_6_dz = -1.0 / 64.0 / dz;

            double smr3 = 59.0 / 48.0 * 64 * dz;
            double smr2 = 43.0 / 48.0 * 64 * dz;
            double smr1 = 49.0 / 48.0 * 64 * dz;
            double spr3 = smr3;
            double spr2 = smr2;
            double spr1 = smr1;

            const int nx = tile_sz[0];
            const int ny = tile_sz[1];
            const int nz = tile_sz[2];

            const int ib = 3;
            const int jb = 3;
            const int kb = 3;
            const int ie = sz[0] - 3;
            const int je = sz[1] - 3;
            const int ke = sz[2] - 3;

            const int n = nx * ny;


            for(unsigned int k=kz_b;k<kz_e;k++)
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        const int pp = IDX(i, j, k);

                        Du[pp] = pre_factor_6_dz *
                                (
                                        -u[pp - 3 * n]
                                        + 6.0 * u[pp - 2 * n]
                                        - 15.0 * u[pp - n]
                                        + 20.0 * u[pp]
                                        - 15.0 * u[pp + n]
                                        + 6.0 * u[pp + 2 * n]
                                        - u[pp + 3 * n]
                                );
                    }

            if(k_b==kb && kz_b==kb)
            {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        Du[IDX(i, j, kz_b)] = pre_factor_6_dz *
                                              (
                                                      -u[IDX(i, j, kz_b + 4)]
                                                      + 6.0 * u[IDX(i, j, kz_b + 3)]
                                                      - 15.0 * u[IDX(i, j, kz_b + 2)]
                                                      + 20.0 * u[IDX(i, j, kz_b + 1)]
                                                      - 15.0 * u[IDX(i, j, kz_b)]
                                                      + 6.0 * u[IDX(i, j, kz_b - 1)]
                                                      - u[IDX(i, j, kz_b - 2)]
                                              );

                    }

            }


            if(k_e==ke && kz_e==ke)
            {
                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        Du[IDX(i, j, kz_e - 1)] = pre_factor_6_dz *
                                                  (
                                                          -u[IDX(i, j, kz_e + 1)]
                                                          + 6.0 * u[IDX(i, j, kz_e)]
                                                          - 15.0 * u[IDX(i, j, kz_e - 1)]
                                                          + 20.0 * u[IDX(i, j, kz_e - 2)]
                                                          - 15.0 * u[IDX(i, j, kz_e - 3)]
                                                          + 6.0 * u[IDX(i, j, kz_e - 4)]
                                                          - u[IDX(i, j, kz_e - 5)]
                                                  );

                    }

            }


            if ((bflag & (1u << OCT_DIR_BACK)) && (k_b==kb && kz_b==kb) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {
                        Du[IDX(i, j, 3)] = (u[IDX(i, j, 6)]
                        - 3.0 * u[IDX(i, j, 5)]
                        + 3.0 * u[IDX(i, j, 4)]
                        - u[IDX(i, j, 3)]
                       ) / smr3;

                        Du[IDX(i, j, 4)] = (
                                                u[IDX(i, j, 7)]
                                                - 6.0 * u[IDX(i, j, 6)]
                                                + 12.0 * u[IDX(i, j, 5)]
                                                - 10.0 * u[IDX(i, j, 4)]
                                                + 3.0 * u[IDX(i, j, 3)]
                                        ) / smr2;
                        Du[IDX(i, j, 5)] = (
                                                u[IDX(i, j, 8)]
                                                - 6.0 * u[IDX(i, j, 7)]
                                                + 15.0 * u[IDX(i, j, 6)]
                                                - 19.0 * u[IDX(i, j, 5)]
                                                + 12.0 * u[IDX(i, j, 4)]
                                                - 3.0 * u[IDX(i, j, 3)]
                                        ) / smr1;

                    }



            }

            if ((bflag & (1u << OCT_DIR_FRONT)) && (k_e==ke && kz_e==ke) ) {

                for(unsigned int j=jy_b;j<jy_e;j++)
                    for(unsigned int i=ix_b;i<ix_e;i++)
                    {

                        Du[IDX(i, j, ke - 3)] = (
                            u[IDX(i, j, ke - 6)]
                            - 6.0 * u[IDX(i, j, ke - 5)]
                            + 15.0 * u[IDX(i, j, ke - 4)]
                            - 19.0 * u[IDX(i, j, ke - 3)]
                            + 12.0 * u[IDX(i, j, ke - 2)]
                            - 3.0 * u[IDX(i, j, ke - 1)]
                    ) / spr1;

                    Du[IDX(i, j, ke - 2)] = (
                                                u[IDX(i, j, ke - 5)]
                                                - 6.0 * u[IDX(i, j, ke - 4)]
                                                + 12.0 * u[IDX(i, j, ke - 3)]
                                                - 10.0 * u[IDX(i, j, ke - 2)]
                                                + 3.0 * u[IDX(i, j, ke - 1)]
                                        ) / spr2;

                    Du[IDX(i, j, ke - 1)] = (
                                                u[IDX(i, j, ke - 4)]
                                                - 3.0 * u[IDX(i, j, ke - 3)]
                                                + 3.0 * u[IDX(i, j, ke - 2)]
                                                - u[IDX(i, j, ke - 1)]
                                        ) / spr3;

                    }

                

            }

#ifdef DEBUG_DERIVS_COMP
            if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
#endif


        }



} //end of namespace cuda




