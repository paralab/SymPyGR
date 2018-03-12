#ifndef RHS_H
#define RHS_H

#include <cmath>
#include <iostream>
#include <time.h>
#include "derivs.h"
#include "profile_param.h"
#include "derivs_cuda.h"




#define IDX(i,j,k) ( (i) + nx * ( (j) + ny * (k) ) )

#define deriv_x deriv42_x
#define deriv_y deriv42_yWithCuda
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

void bssnrhs(double * dev_var_in, const unsigned int unzip_dof,
             const unsigned int &offset,
             const double *ptmin, const double *ptmax, const unsigned int *sz,
             const unsigned int &bflag);


void bssn_bcs(double *f_rhs, const double *f,
              const double *dxf, const double *dyf, const double *dzf,
              const double *pmin, const double *pmax,
              const double f_falloff, const double f_asymptotic,
              const unsigned int *sz, const unsigned int &bflag);

void fake_initial_data(double x, double y, double z, double *u);

#endif
