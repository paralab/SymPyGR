#ifndef _DERVS_H
#define _DERVS_H

#include <cmath>
#include "dendro.h"


#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

// stores the padding width. 
extern unsigned int DERIV_PW;


void deriv42_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);
void deriv42_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);


void deriv42adv_z(double * const  Dzu, const double * const  u,const double dz, const unsigned int *sz, const double * const betaz, unsigned bflag);
void deriv42adv_y(double * const  Dyu, const double * const  u,const double dy, const unsigned int *sz, const double * const betay, unsigned bflag);
void deriv42adv_x(double * const  Dxu, const double * const  u,const double dx, const unsigned int *sz, const double * const betax, unsigned bflag);

void deriv42_zz(double * const  Du, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);
void deriv42_yy(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int *sz, unsigned bflag);


void ko_deriv42_z(double * const Du, const double * const u, const double dz, const unsigned *sz, unsigned bflag);
void ko_deriv42_y(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void ko_deriv42_x(double * const  Du, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);

void disstvb3_x(double * const  Du, const double * const  u, const double * const lam, const double dx, const unsigned int *sz, unsigned bflag);
void disstvb3_y(double * const  Du, const double * const  u, const double * const lam, const double dx, const unsigned int *sz, unsigned bflag);
void disstvb3_z(double * const  Du, const double * const  u, const double * const lam, const double dx, const unsigned int *sz, unsigned bflag);

void disstvb5_x(double * const  Du, const double * const  u, const double * const lam, const double dx, const unsigned int *sz, unsigned bflag);
void disstvb5_y(double * const  Du, const double * const  u, const double * const lam, const double dx, const unsigned int *sz, unsigned bflag);
void disstvb5_z(double * const  Du, const double * const  u, const double * const lam, const double dx, const unsigned int *sz, unsigned bflag);

/**@brief: copies unzipped padding to the computed derivatives.
 * (this is essential for the mixed derivatives.)
 * */
void cpy_unzip_padd(double * const  Du, const double * const  u,const unsigned int *sz, unsigned bflag);


// 6th order derivatives. 
void deriv64_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);
void deriv64_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv64_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);

void deriv64_zz(double * const  Du, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);
void deriv64_yy(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv64_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int *sz, unsigned bflag);

void deriv64adv_z(double * const  Dzu, const double * const  u,const double dz, const unsigned int *sz, const double * const betaz, unsigned bflag);
void deriv64adv_y(double * const  Dyu, const double * const  u,const double dy, const unsigned int *sz, const double * const betay, unsigned bflag);
void deriv64adv_x(double * const  Dxu, const double * const  u,const double dx, const unsigned int *sz, const double * const betax, unsigned bflag);


// these are the derivs that will be used in the BSSN CODE based on the FD derivative order. 



#ifdef BSSN_USE_4TH_ORDER_DERIVS
    #define deriv_x deriv42_x
    #define deriv_y deriv42_y
    #define deriv_z deriv42_z

    #define deriv_xx deriv42_xx
    #define deriv_yy deriv42_yy
    #define deriv_zz deriv42_zz

    #define adv_deriv_x deriv42adv_x
    #define adv_deriv_y deriv42adv_y
    #define adv_deriv_z deriv42adv_z

    #define ko_deriv_x ko_deriv42_x
    #define ko_deriv_y ko_deriv42_y
    #define ko_deriv_z ko_deriv42_z
#endif


#ifdef BSSN_USE_6TH_ORDER_DERIVS
    #define deriv_x deriv64_x
    #define deriv_y deriv64_y
    #define deriv_z deriv64_z

    #define deriv_xx deriv64_xx
    #define deriv_yy deriv64_yy
    #define deriv_zz deriv64_zz

    #define adv_deriv_x deriv64adv_x
    #define adv_deriv_y deriv64adv_y
    #define adv_deriv_z deriv64adv_z

    #define ko_deriv_x ko_deriv42_x
    #define ko_deriv_y ko_deriv42_y
    #define ko_deriv_z ko_deriv42_z
#endif




#endif
