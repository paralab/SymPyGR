#include <iostream>
#include <cmath>
#include <assert.h>
#include "hawaii.h"
#include "bssn.h"

double t0;
double tf;

void rhs_stage1(coll_point_t *point, const int gen);
void rhs_stage2(coll_point_t *point, const int gen);
void rhs_stage3(coll_point_t *point, const int gen);

/// ---------------------------------------------------------------------------
/// Solving the BSSN formulation of the Einstein equations.
/// ---------------------------------------------------------------------------

// get_local_dt  {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
double get_local_dt (const coll_point_t *point)
{
  double dx0 = L_dim[x_dir] / max_index[x_dir];
  double dy0 = L_dim[y_dir] / max_index[y_dir];
  double dz0 = L_dim[z_dir] / max_index[z_dir];

  int h = 1 << (JJ - point->level);
  double dx = dx0 * h;
  double dy = dy0 * h;
  double dz = dz0 * h;
  double delta = fmin(dx, dy);
  delta = fmin(delta, dz);
  double dt = pars.cfl * delta;
  return dt;
}
// }}}

// problem_init {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void problem_init(void) 
{
  t0 = 0.0; // Simulation starting time

  // FIXME set simulation stopping time with a parameter
  tf = 1;   // Simulation stopping time

  double ranges[n_variab];
  initial_ranges(ranges);

  // FIXME Check definition of epsilon
  for (int i=0;i<n_variab;i++) {
    eps[i] = pars.epsilon; // Error tolerance for grid refinement
    eps_scale[i] = ranges[1] * eps[i]; // Error tolerance for grid refinement
  }

  assert(n_dim == 3);
  L_dim[x_dir] = pars.xmax - pars.xmin; // Problem physical dimension
  L_dim[y_dir] = pars.ymax - pars.ymin; // Problem physical dimension
  L_dim[z_dir] = pars.zmax - pars.zmin; // Problem physical dimension

  validate_parameters();
  time_stamp = 0; // Initialize global time stamp
}
// }}}

// {{{
/*----------------------------------------------------------------------
 *
 * initial_ranges:
 *                 This subroutine makes a guess for the initial range
 *                 of the BSSN variables.
 *
 *----------------------------------------------------------------------*/
void initial_ranges(double *ranges) {

  for (int i = 0; i < n_variab; ++i) {
    ranges[i] = 1.0;
  }
}
// }}}

// initial_condition {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void initial_condition(const coord_t *coord, double *u) 
{
  if (pars.id_type == 0) {
    init_data_puncture(coord, u);
  }
  else if (pars.id_type == 1) {
    init_data_fake(coord, u);
  }
  else if (pars.id_type == 2) {
    init_data_kerrschild(coord, u);
  }
  //else if (pars.id_type == 3) {
    //init_teukosky_waves(coord, u);
  //}
  else {
    printf("Unknown initial data type = %d\n",pars.id_type);
  }
}
// }}}

// init_data_fake {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void init_data_fake(const coord_t *coord, double *u) 
{
 /* This is fake data just for testing. */
  
  double x1 = coord->pos[x_dir];
  double y1 = coord->pos[y_dir];
  double z1 = coord->pos[z_dir];

  
  u[U_ALPHA]  =  1.0 - 0.25*sin(x1);
  u[U_SHIFTX] =  0.0;
  u[U_SHIFTY] =  0.0;
  u[U_SHIFTZ] =  0.0;

  u[U_GAMTX]  =  5.0*cos(x1)/(10.0*sin(x1)+26.0-1.0*pow(cos(x1),2));
  u[U_GAMTY]  = -5.0*sin(y1)/(25.0+10.0*cos(y1)+pow(cos(y1),2));
  u[U_GAMTZ]  =  0.0;

  u[U_GBX]    =  0.0;
  u[U_GBY]    =  0.0;
  u[U_GBZ]    =  0.0;

  u[U_GTXX]   =  1.0+0.2*sin(x1);
  u[U_GTXY]   =  0.0;
  u[U_GTXZ]   =  0.0;
  u[U_GTYY]   =  1.0+0.2*cos(y1);
  u[U_GTYZ]   =  0.0;
  u[U_GTZZ]   =  1.0/(u[U_GTXX]*u[U_GTYY]);

  double Atd[3][3];
  Atd[0][0] = exp(-4.0*cos(x1)*sin(y1))
            * (cos(x1)-0.3333333333e0*exp(4.0*cos(x1)*sin(y1))
            * (1.0+0.2*sin(x1))*(5.0*exp(-4.0*cos(x1)*sin(y1))
            / (5.0+sin(x1))*cos(x1)+5.0*exp(-4.0*cos(x1)*sin(y1))
            / (5.0+cos(y1))*cos(y1)+0.4e-1*(25.0+5.0*cos(y1) + 5.0*sin(x1)
            + sin(x1)*cos(y1))*exp(-4.0*cos(x1)*sin(y1))*cos(z1)));
  Atd[0][1] = 0.0;
  Atd[0][2] = 0.0;
  Atd[1][0] = 0.0;
  Atd[1][1] = exp(-0.4e1*cos(x1)*sin(y1))
            * (cos(y1)-0.3333333333e0*exp(4.0*cos(x1)*sin(y1))
            * (1.0+0.2*cos(y1))*(5.0*exp(-4.e0*cos(x1)*sin(y1))
            / (5.e0+sin(x1))*cos(x1)+5.e0*exp(-4.e0*cos(x1)*sin(y1))
            /(5.e0+cos(y1))*cos(y1)+0.4e-1*(25.e0+5.e0*cos(y1)+5.e0*sin(x1)
            +sin(x1)*cos(y1))*exp(-4.e0*cos(x1)*sin(y1))*cos(z1)));
  Atd[1][2] = 0.0;
  Atd[2][0] = 0.0;
  Atd[2][1] = 0.0;
  Atd[2][2] = exp(-0.4e1*cos(x1)*sin(y1))
            * (cos(z1)-0.3333333333e0*exp(4.0*cos(x1)*sin(y1))
            / (1.0+0.2*sin(x1))/(1.0+0.2*cos(y1))
            * (5.e0*exp(-4.0*cos(x1)*sin(y1))/(5.0+sin(x1))*cos(x1)
            + 5.e0*exp(-4.0*cos(x1)*sin(y1))/(5.0+cos(y1))*cos(y1)
            + 0.4e-1*(25.e0+5.e0*cos(y1)+5.e0*sin(x1)+sin(x1)*cos(y1))
            * exp(-4.e0*cos(x1)*sin(y1))*cos(z1)));

  u[U_ATXX] = Atd[0][0];
  u[U_ATXY] = Atd[0][1];
  u[U_ATXZ] = Atd[0][2];
  u[U_ATYY] = Atd[1][1];
  u[U_ATYZ] = Atd[1][2];
  u[U_ATZZ] = Atd[2][2];

  u[U_CHI] = exp(-4.0*cos(x1)*sin(y1));
  u[U_TRK] = 5.0*exp(-4.0*cos(x1)*sin(y1))/(5.0+sin(x1))*cos(x1)
           + 5.0*exp(-4.0*cos(x1)*sin(y1))/(5.0+cos(y1))*cos(y1)
           + 0.4e-1*(25.0+5.0*cos(y1)+5.0*sin(x1)
           + sin(x1)*cos(y1))*exp(-4.0*cos(x1)*sin(y1))*cos(z1);

}
// }}}

// init_data_kerrschild (doesn't exist yet...) {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void init_data_kerrschild(const coord_t *coord, double *u) 
{
  printf("Kerr-Schild data not ready yet.\n");
  exit(2);
}
// }}}

// init_data_puncture {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void init_data_puncture(const coord_t *coord, double *u) 
{
  double epijk[3][3][3];
  int i,j,k;
  for (k=0;k<3;k++) {
    for (j=0;j<3;j++) {
      for (i=0;i<3;i++) {
        epijk[k][j][i] = 0.0;
      }
    }
  }
  epijk[0][1][2] = 1.0;epijk[1][2][0] = 1.0;epijk[2][0][1] = 1.0;
  epijk[0][2][1] = -1.0;epijk[2][1][0] = -1.0;epijk[1][0][2] = -1.0;

  double deltaij[3][3];
  for (j=0;j<3;j++) {
    for (i=0;i<3;i++) {
      deltaij[j][i] = 0.0;
    }
  }
  deltaij[0][0] = 1.0;deltaij[1][1] = 1.0;deltaij[2][2] = 1.0;

  double mass1 = pars.mass1;
  double bh1x = pars.bh1x;
  double bh1y =  pars.bh1y;
  double bh1z =  pars.bh1z;
  double vp1[3];
  vp1[0] = pars.vx1;
  vp1[1] = pars.vy1;
  vp1[2] = pars.vz1;
  double vp1tot = sqrt( vp1[0]*vp1[0] + vp1[1]*vp1[1] + vp1[2]*vp1[2] );
  double spin1 = pars.spin1;
  double spin1_th = pars.spin1_th;
  double spin1_phi = pars.spin1_phi;
  double vs1[3];
  vs1[0] = spin1*sin(spin1_th)*cos(spin1_phi);
  vs1[1] = spin1*sin(spin1_th)*sin(spin1_phi);
  vs1[2] = spin1*cos(spin1_th);

  // bh 2
  double mass2 = pars.mass2;
  double bh2x =  pars.bh2x;
  double bh2y =  pars.bh2y;
  double bh2z =  pars.bh2z;
  double vp2[3];
  vp2[0] = pars.vx2;
  vp2[1] = pars.vy2;
  vp2[2] = pars.vz2;
  double vp2tot = sqrt( vp2[0]*vp2[0] + vp2[1]*vp2[1] + vp2[2]*vp2[2] );
  double spin2 = pars.spin2;
  double spin2_th = pars.spin2_th;
  double spin2_phi = pars.spin2_phi;
  double vs2[3];
  vs2[0] = spin2*sin(spin2_th)*cos(spin2_phi);
  vs2[1] = spin2*sin(spin2_th)*sin(spin2_phi);
  vs2[2] = spin2*cos(spin2_th);

  double xx = coord->pos[x_dir];
  double yy = coord->pos[y_dir];
  double zz = coord->pos[z_dir];

  double x1,y1,z1,rv1;
  double x2,y2,z2,rv2;
  double vn1[3],vn2[3];
  double vpsibl;
  double v_u_corr,amp_capj,amp_capr,l_r,u0_j,u2_j,mu_j,p2_mu_j,v_u_j1;
  double v1,v2,v3,v4,vt1,vt2;
  int i1,i2,i3,i4;
  double amp_capp,u0_p,u2_p,mu_p,p2_mu_p;
  double v_u_p1,v_u_c1,v_u_j2,v_u_p2;
  double v_u_c2,vpsibl_u,vpsibl_u2;

  //allocating the location of BH1

  x1 = xx - bh1x;
  y1 = yy - bh1y;
  z1 = zz - bh1z;

  //locating as a radial form
  rv1 = sqrt(x1*x1 + y1*y1 + z1*z1);
  vn1[0] = x1/rv1;
  vn1[1] = y1/rv1;
  vn1[2] = z1/rv1;

  //same as BH2
  x2 = xx - bh2x;
  y2 = yy - bh2y;
  z2 = zz - bh2z;

  rv2 = sqrt(x2*x2 + y2*y2 + z2*z2);
  vn2[0] = x2/rv2;
  vn2[1] = y2/rv2;
  vn2[2] = z2/rv2;

  //Initial data is related with the paper: http://arxiv.org/abs/0711.1165
  //Brill-Lindquist conformal factor
  vpsibl = 1.0 + mass1/(2.0*rv1);
  vpsibl = vpsibl + mass2/(2.0*rv2);

  v_u_corr = 0.0;
  // bh 1
  //For spinning puncture
  if ( fabs(spin1) > 1.e-6 ) {
      amp_capj = 4.0*spin1/(mass1*mass1);
      amp_capr = 2.0*rv1/mass1;
      l_r = 1.0/(1.0+amp_capr);
      u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
      u2_j = -pow(l_r,5)/20.0;
      mu_j = vn1[0]*vs1[0];
      mu_j = mu_j + vn1[1]*vs1[1];
      mu_j = (mu_j + vn1[2]*vs1[2])/fabs(spin1);
      p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
      v_u_j1 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
      v_u_corr = v_u_corr + v_u_j1;
  }
  //For boosting puncture
  if (vp1tot > 1.e-6) {
      amp_capp = 2.0*vp1tot/mass1;
      amp_capr = 2.0*rv1/mass1;
      l_r = 1.0/(1.0 + amp_capr);
      u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
      u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
      u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
      u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
      u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
      u2_p = (u2_p)/(80.0*amp_capr);
      mu_p =        vn1[0]*vp1[0]/vp1tot;
      mu_p = mu_p + vn1[1]*vp1[1]/vp1tot;
      mu_p = mu_p + vn1[2]*vp1[2]/vp1tot;
      p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
      v_u_p1 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
      v_u_corr = v_u_corr + v_u_p1;
  }
  //For spinning boosted pucture
  if ( vp1tot > 1.e-6 && fabs(spin1) > 1.e-6 ) {
        v1 =      (vp1[1]*vs1[2]-vp1[2]*vs1[1])*vn1[0];
        v1 = v1 + (vp1[2]*vs1[0]-vp1[0]*vs1[2])*vn1[1];
        v1 = v1 + (vp1[0]*vs1[1]-vp1[1]*vs1[0])*vn1[2];
        v1 = v1*(16.0/pow(mass1,4))*rv1;
        
        amp_capr = 2.0*rv1/mass1;
        l_r = 1.0/(1.0 + amp_capr);
        
        v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);
        
        v_u_c1 = (v1*v2*pow(l_r,5))/80.0;
        v_u_corr = v_u_corr + v_u_c1;
  }
  // bh 2 same puncture as bh 1
  if ( fabs(spin2) > 1.e-6 ) {
        amp_capj = 4.0*spin2/(mass2*mass2);
        amp_capr = 2.0*rv2/mass2;
        l_r = 1.0/(1.0+amp_capr);
        u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
        u2_j = -pow(l_r,5)/20.0;
        mu_j = vn2[0]*vs2[0];
        mu_j = mu_j + vn2[1]*vs2[1];
        mu_j = (mu_j + vn2[2]*vs2[2])/fabs(spin2);
        p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
        v_u_j2 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
        v_u_corr = v_u_corr + v_u_j2;
  }
    
  if (vp2tot > 1.e-6) {
        amp_capp = 2.0*vp2tot/mass2;
        amp_capr = 2.0*rv2/mass2;
        l_r = 1.0/(1.0 + amp_capr);
        u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
        u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
        u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
        u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
        u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
        u2_p = (u2_p)/(80.0*amp_capr);
        mu_p =        vn2[0]*vp2[0]/vp2tot;
        mu_p = mu_p + vn2[1]*vp2[1]/vp2tot;
        mu_p = mu_p + vn2[2]*vp2[2]/vp2tot;
        p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
        v_u_p2 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
        v_u_corr = v_u_corr + v_u_p2;
  }
    
  if ( vp2tot > 1.e-6 && fabs(spin2) > 1.e-6 ) {
        v1 =      (vp2[1]*vs2[2]-vp2[2]*vs2[1])*vn2[0];
        v1 = v1 + (vp2[2]*vs2[0]-vp2[0]*vs2[2])*vn2[1];
        v1 = v1 + (vp2[0]*vs2[1]-vp2[1]*vs2[0])*vn2[2];
        v1 = v1*(16.0/pow(mass2,4))*rv2;
        
        amp_capr = 2.0*rv2/mass2;
        l_r = 1.0/(1.0 + amp_capr);
        
        v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);
        
        v_u_c2 = (v1*v2*pow(l_r,5))/80.0;
        v_u_corr = v_u_corr + v_u_c2;
  }
 
  // vpsibl_u will be used for the conformal factor,
  vpsibl_u  = vpsibl + v_u_corr;
  // vpsibl_u2 is for the Aij terms...
  // ! since the corrections are first order...
  // ! adding half of the correction seems to give the best results...
  // ! update - do a fit for spin = 0.6...
  vpsibl_u2 = vpsibl + v_u_corr;

  u[U_ALPHA] = 1.0/(vpsibl_u*vpsibl_u);
  v2 = 1.0/pow(vpsibl_u,4);
  u[U_CHI] = v2;
  u[U_TRK] = 0.0;

  u[U_SHIFTX] = 0.0;
  u[U_SHIFTY] = 0.0;
  u[U_SHIFTZ] = 0.0;

  u[U_GAMTX] = 0.0;
  u[U_GAMTY] = 0.0;
  u[U_GAMTZ] = 0.0;

  u[U_GBX] = 0.0;
  u[U_GBY] = 0.0;
  u[U_GBZ] = 0.0;

  u[U_GTXX] = 1.0;
  u[U_GTXY] = 0.0;
  u[U_GTXZ] = 0.0;
  u[U_GTYY] = 1.0;
  u[U_GTYZ] = 0.0;
  u[U_GTZZ] = 1.0;

  for (i1=0;i1<3;i1++) {
    for (i2=0;i2<3;i2++) {
      // first BH
      v2 = 0.0;
      for (i3=0;i3<3;i3++) {
        for (i4=0;i4<3;i4++) {
          vt1 = epijk[i1][i3][i4]*vs1[i3]*vn1[i4]*vn1[i2];
          vt2 = epijk[i2][i3][i4]*vs1[i3]*vn1[i4]*vn1[i1];
          v2 = v2 + vt1 + vt2;
        }
      }

      v3 = vp1[i1]*vn1[i2] + vp1[i2]*vn1[i1];
      vt1 = 0.0;
      for (i3=0;i3<3;i3++) {
        vt1 = vt1 + vp1[i3]*vn1[i3];
      }
      vt1 = vt1*(vn1[i1]*vn1[i2] - deltaij[i1][i2]);
      v3 = v3 + vt1;

      v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv1,3));
      v4 = v1*(v2+(rv1/2.0)*v3);

      // second BH
      v2 = 0.0;
      for (i3=0;i3<3;i3++) {
        for (i4=0;i4<3;i4++) {
          vt1 = epijk[i1][i3][i4]*vs2[i3]*vn2[i4]*vn2[i2];
          vt2 = epijk[i2][i3][i4]*vs2[i3]*vn2[i4]*vn2[i1];
          v2 = v2 + vt1 + vt2;
        }
      }

      v3 = vp2[i1]*vn2[i2] + vp2[i2]*vn2[i1];
      vt1 = 0.0;
      for (i3=0;i3<3;i3++) {
        vt1 = vt1 + vp2[i3]*vn2[i3];
      }
      vt1 = vt1*(vn2[i1]*vn2[i2] - deltaij[i1][i2]);
      v3 = v3 + vt1;

      v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv2,3));
      v4 = v4 + v1*(v2+(rv2/2.0)*v3);

      if ( i1 == 0 && i2 == 0 ) {
        u[U_ATXX] = v4;
      } else if ( i1 == 0 && i2 == 1 ) {
        u[U_ATXY] = v4;
      } else if ( i1 == 0 && i2 == 2 ) {
        u[U_ATXZ] = v4;
      } else if ( i1 == 1 && i2 == 1 ) {
        u[U_ATYY] = v4;
      } else if ( i1 == 1 && i2 == 2 ) {
        u[U_ATYZ] = v4;
      } else if ( i1 == 2 && i2 == 2 ) {
        u[U_ATZZ] = v4;
      }

    }
  }

}
// }}}

// compute_rhs {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void compute_rhs(const double t, const int gen) {

 /* compute finite difference derivatives everywhere */
  for (int i = 0; i < npts_in_array; i++) {
    coll_point_t *point = &coll_points->array[i];
    rhs_stage1(point, gen);
  }
  for (int i = 0; i < HASH_TBL_SIZE; ++i) {
    hash_entry_t *ptr = coll_points->hash_table[i];
    while (ptr != NULL) {
      if (ptr->point->status[CURRENT_STATUS] > nonessential) {
        rhs_stage1(ptr->point, gen);
      }
      ptr = ptr->next;
    }
  }

 /* call the RHS */
  for (int i = 0; i < npts_in_array; i++) {
    coll_point_t *point = &coll_points->array[i];
    rhs_stage2(point, gen);
  }

  for (int i = 0; i < HASH_TBL_SIZE; ++i) {
    hash_entry_t *ptr = coll_points->hash_table[i];
    while (ptr != NULL) {
      if (ptr->point->status[CURRENT_STATUS] > nonessential) {
        rhs_stage2(ptr->point, gen);
      }
      ptr = ptr->next;
    }
  }

 /* compute KO dissipation */
  for (int i = 0; i < npts_in_array; i++) {
    coll_point_t *point = &coll_points->array[i];
    rhs_stage3(point, gen);
  }
  for (int i = 0; i < HASH_TBL_SIZE; ++i) {
    hash_entry_t *ptr = coll_points->hash_table[i];
    while (ptr != NULL) {
      if (ptr->point->status[CURRENT_STATUS] > nonessential) {
        rhs_stage3(ptr->point, gen);
      }
      ptr = ptr->next;
    }
  }
}
// }}}

// rhs_stage1 {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void rhs_stage1(coll_point_t *point, const int gen) {
    if (gen > 0)
      check_derivative_stencil(point, gen);

    // FIXME We aren't using a mask currently.
    int mask[n_variab + n_aux] = {0};
    for (int ivar = 0; ivar < n_variab; ivar++)
      mask[ivar] = 1;

   /* Now calling Finite Difference routine for derivatives */
    for (dir_t dir = x_dir; dir < n_dim; dir++)
      compute_fd_derivative(point, gen, dir, mask);

}
// }}}

// rhs_stage2 {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void rhs_stage2(coll_point_t *point, const int gen) {
    for (int ivar = 0; ivar < n_variab; ivar++) {
      point->rhs[gen][0] = 0.0;
    }

   /* call RHS routine */
    bssn_rhs(point, gen);
}
// }}}

// rhs_stage3 {{{
/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void rhs_stage3(coll_point_t *point, const int gen) {
   /* FIXME We aren't using a mask currently. This will determine which fields
    *       get K-O dissipation. 
    */
    int mask[n_variab + n_aux] = {0};
    for (int ivar = 0; ivar < n_variab; ivar++)
      mask[ivar] = 1;

   /* Compute K-O dissipation derivatives in each direction. 
    * K-O derivatives are stored in the standard location for derivatives.
    */
    for (dir_t dir = x_dir; dir < n_dim; dir++)
      compute_KOdiss(point, gen, dir, mask);

    apply_KOdiss(point, gen, mask);

}
// }}}
