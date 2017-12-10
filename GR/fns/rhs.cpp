#include <iostream>
#include <cmath>
#include <cstring>
#include <assert.h>
#include "hawaii.h"

#include <stdint.h>
#include <time.h>

// bssn_rhs {{{
/*---------------------------------------------------------------------------*
 *
 *
 *
 *---------------------------------------------------------------------------*/
void bssn_rhs(coll_point_t *point, const int gen)
{

  double rho_ADM = 0.0;
  double Jtd_ADM[3] = {0.0, 0.0, 0.0};
  double Tu[4][4] = { {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, 
                      {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0} };
  double pTtd_ADM[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
  double tr_pT = 0.0;

  double Alpha, chi, trK, inv_chi;
  double Betau[3];
  double gtd[3][3], gtu[3][3], Atd[3][3], Atu[3][3], Atud[3][3];
  double Bu[3];
  double Rtd[3][3], Rpd_1[3][3];
  double Ct[3][3][3], Ctd[3][3][3], Psi1[3][3], Psi1TF[3][3];
  double CalGamt[3];

  double third_trPsi1, idetgtd, detgtd, div_Beta;
  double dd_Alpha[3][3];
  double dd_Betau[3][3][3], d_div_Beta[3];
  double dd_chi[3][3];
  double dd_gtd[3][3][3][3];

  double gtd_rhs[3][3], Atd_rhs[3][3];
  double chi_rhs, trK_rhs;
  double Gamt_rhs[3];
  double Alpha_rhs, Betau_rhs[3], Bu_rhs[3];

  const double four_pi_G  = 4.0*acos(-1.0);
  const double eight_pi_G  = 8.0*acos(-1.0);
  const double fourth = 0.25;
  const double third = 1.0/3.0;
  const double half = 0.5;
  const double twothirds = 2.0/3.0;
  const double threefourths = 0.75;
  const double threehalves = 1.5;
  const double two = 2.0;

  const int lambda_1 = pars.lambda_1;
  const int lambda_2 = pars.lambda_2;
  const int lambda_3 = pars.lambda_3;
  const int lambda_4 = pars.lambda_4;
  const int lambda_f0 = pars.lambda_f0;
  const int lambda_f1 = pars.lambda_f1;
  const double trK0 = pars.trK0;
  const double feta = pars.eta;
  
  double MapleGenVar1, MapleGenVar2;


  double *rhs = point->rhs[gen];


#if 0
  printf(" time_stamp=%d, gen=%d\n",point->time_stamp,gen);
  printf(" idx = (%d, %d, %d), (x,y,z)=( %g, %g, %g)\n",
           point->index.idx[0], point->index.idx[1], point->index.idx[2],
           point->coords.pos[0],point->coords.pos[1],point->coords.pos[2]);
#endif
 

  index_t id;
  for (int i = x_dir; i < n_dim; i++) {
    id.idx[i] = point->index.idx[i];
  }

  double dx[n_dim];
  int h[n_dim];

  for (int dir = x_dir; dir < n_dim; dir++) {
    int closest_level = get_closest_level(point, dir);
    assert(closest_level  >= 1 );
    dx[dir] = L_dim[dir] / ns[dir] / (1 << closest_level);
    h[dir] = 1 << (JJ - closest_level);
  }

 /* 
  * Build the stencil. This is the original choice
  */
  // FIXME: This whole section needs to be rewritten.
  double *u[3][7], *du[3][7], *pos[3][7];
  coll_point_t *stcl[3][7];

  int nnull = 0;
  for (int dir = 0; dir < n_dim; dir++) {
    for (int qq = x_dir; qq < n_dim; qq++) {
      id.idx[qq] = point->index.idx[qq];
    }
    for (int ll = 0; ll < 7; ll++) {
      if ( ll == 3 ) {
        u[dir][ll] = point->u[gen];
        du[dir][ll] = point->du[gen];
        pos[dir][ll] = point->coords.pos;
        stcl[dir][ll] = point;
      } else {
        coll_point_t *cpm=NULL;
        id.idx[dir] = point->index.idx[dir] + (ll-3)*h[dir];
        if (check_index(&id)) {
          cpm = get_coll_point(&id);
          //Stage 1 should have assured that all points exist, so...
          if (cpm == NULL) {
            printf("...trap me...\n");
            printf("     RHS point = (%d, %d, %d)\n",point->index.idx[0],point->index.idx[1],point->index.idx[2]);
            printf("     dir = %d, id.idx = (%d, %d, %d)\n",dir,id.idx[0],id.idx[1],id.idx[2]);
            printf("     level = %d, max_level = %d, ll=%d\n",point->level, max_level,ll);
            printf("     status= %d",point->status[0]);
          }
          assert(cpm != NULL);
  
          u[dir][ll] = cpm->u[gen];
          du[dir][ll] = cpm->du[gen];
          stcl[dir][ll] = cpm;
          //TODO, this should actually be the closest position to the center;
          // this will be corrected below
          pos[dir][ll] = cpm->coords.pos;
        } else {
          u[dir][ll] = NULL;
          du[dir][ll] = NULL;
          pos[dir][ll] = NULL;
          stcl[dir][ll] = NULL;
          nnull++;
        }
      }
    }
  }

  int bi[3] = {0, 0, 0};
  int npts[3] = {7, 7, 7}; 
  int indx[3] = {3, 3, 3};

  if (nnull > 0) {
   /* Apply boundary conditions and return */
    for (int d = 0; d < 3; d++) {
     /* check to see if neighbors don't exist. If so, then we are on boundary */
      if (u[d][2] == NULL || u[d][4]== NULL) {
        sommerfeld_boundary(point, gen);
        return;
      }

      if (u[d][1] == NULL) {
        bi[d]   = 2;
        npts[d] = 5;
        indx[d] = 1;
      }
      else if (u[d][0] == NULL) {
        bi[d]   = 1;
        npts[d] = 6; 
        indx[d] = 2;
      }
      else if (u[d][5] == NULL) {
        bi[d]   = 0;
        npts[d] = 5;
        indx[d] = 3;
      }
      else if (u[d][6] == NULL) {
        bi[d]   = 0;
        npts[d] = 6;
        indx[d] = 3;
      }
    } 
}
if (nnull > 0){
    //* check to see if neighbors don't exist. If so, then we are on boundary *
     for(int dd = 0; dd < 3; dd++) {
      if (du[dd][2] == NULL || du[dd][4]== NULL) {
        sommerfeld_boundary(point, gen);
        return;
      }
/*
      if (du[dd][1] == NULL) {
        bi[dd]   = 2;
        npts[dd] = 5;
        indx[dd] = 1;
      }
      else if (du[dd][0] == NULL) {
        bi[dd]   = 1;
        npts[dd] = 6; 
        indx[dd] = 2;
      }
      else if (du[dd][5] == NULL) {
        bi[dd]   = 0;
        npts[dd] = 5;
        indx[dd] = 3;
      }
      else if (du[dd][6] == NULL) {
        bi[dd]   = 0;
        npts[dd] = 6;
        indx[dd] = 3;
      } */
    } 
  }

  double *upt = point->u[gen];
  double *dupt = point->du[gen];

  Alpha = upt[U_ALPHA];
  Betau[0] = upt[U_SHIFTX];
  Betau[1] = upt[U_SHIFTY];
  Betau[2] = upt[U_SHIFTZ];
  Bu[0] = upt[U_GBX];
  Bu[1] = upt[U_GBY];
  Bu[2] = upt[U_GBZ];

  //Gamt[0] = upt[U_GAMTX];
  //Gamt[1] = upt[U_GAMTY];
  //Gamt[2] = upt[U_GAMTZ];

  chi = upt[U_CHI];
  trK = upt[U_TRK];
  gtd[0][0] = upt[U_GTXX];
  gtd[0][1] = upt[U_GTXY];
  gtd[0][2] = upt[U_GTXZ];
  gtd[1][0] = upt[U_GTXY];
  gtd[1][1] = upt[U_GTYY];
  gtd[1][2] = upt[U_GTYZ];
  gtd[2][0] = upt[U_GTXZ];
  gtd[2][1] = upt[U_GTYZ];
  gtd[2][2] = upt[U_GTZZ];

  Atd[0][0] = upt[U_ATXX];
  Atd[0][1] = upt[U_ATXY];
  Atd[0][2] = upt[U_ATXZ];
  Atd[1][0] = upt[U_ATXY];
  Atd[1][1] = upt[U_ATYY];
  Atd[1][2] = upt[U_ATYZ];
  Atd[2][0] = upt[U_ATXZ];
  Atd[2][1] = upt[U_ATYZ];
  Atd[2][2] = upt[U_ATZZ];

 /* define derivatives first */
 /* metric first derivatives */
  double d_chi[3];
  double d_gtd[3][3][3];

  d_chi[0] = dupt[U_CHI*n_dim + 0];
  d_chi[1] = dupt[U_CHI*n_dim + 1];
  d_chi[2] = dupt[U_CHI*n_dim + 2];

  d_gtd[0][0][0] = dupt[U_GTXX*n_dim + 0];
  d_gtd[1][0][0] = dupt[U_GTXX*n_dim + 1];
  d_gtd[2][0][0] = dupt[U_GTXX*n_dim + 2];

  d_gtd[0][0][1] = dupt[U_GTXY*n_dim + 0];
  d_gtd[1][0][1] = dupt[U_GTXY*n_dim + 1];
  d_gtd[2][0][1] = dupt[U_GTXY*n_dim + 2];

  d_gtd[0][0][2] = dupt[U_GTXZ*n_dim + 0];
  d_gtd[1][0][2] = dupt[U_GTXZ*n_dim + 1];
  d_gtd[2][0][2] = dupt[U_GTXZ*n_dim + 2];

  d_gtd[0][1][0] = d_gtd[0][0][1];
  d_gtd[1][1][0] = d_gtd[1][0][1];
  d_gtd[2][1][0] = d_gtd[2][0][1];

  d_gtd[0][1][1] = dupt[U_GTYY*n_dim + 0];
  d_gtd[1][1][1] = dupt[U_GTYY*n_dim + 1];
  d_gtd[2][1][1] = dupt[U_GTYY*n_dim + 2];

  d_gtd[0][1][2] = dupt[U_GTYZ*n_dim + 0];
  d_gtd[1][1][2] = dupt[U_GTYZ*n_dim + 1];
  d_gtd[2][1][2] = dupt[U_GTYZ*n_dim + 2];

  d_gtd[0][2][0] = d_gtd[0][0][2];
  d_gtd[1][2][0] = d_gtd[1][0][2];
  d_gtd[2][2][0] = d_gtd[2][0][2];

  d_gtd[0][2][1] = d_gtd[0][1][2];
  d_gtd[1][2][1] = d_gtd[1][1][2];
  d_gtd[2][2][1] = d_gtd[2][1][2];

  d_gtd[0][2][2] = dupt[U_GTZZ*n_dim + 0];
  d_gtd[1][2][2] = dupt[U_GTZZ*n_dim + 1];
  d_gtd[2][2][2] = dupt[U_GTZZ*n_dim + 2];

 /* trK first derivatives */
  double d_trK[3];
  d_trK[0] = dupt[U_TRK*n_dim + 0];
  d_trK[1] = dupt[U_TRK*n_dim + 1];
  d_trK[2] = dupt[U_TRK*n_dim + 2];

 /* Gamt first derivatives */
  double d_Gamt[3][3];
  d_Gamt[0][0] = dupt[U_GAMTX*n_dim + 0];
  d_Gamt[1][0] = dupt[U_GAMTX*n_dim + 1];
  d_Gamt[2][0] = dupt[U_GAMTX*n_dim + 2];

  d_Gamt[0][1] = dupt[U_GAMTY*n_dim + 0];
  d_Gamt[1][1] = dupt[U_GAMTY*n_dim + 1];
  d_Gamt[2][1] = dupt[U_GAMTY*n_dim + 2];

  d_Gamt[0][2] = dupt[U_GAMTZ*n_dim + 0];
  d_Gamt[1][2] = dupt[U_GAMTZ*n_dim + 1];
  d_Gamt[2][2] = dupt[U_GAMTZ*n_dim + 2];


 /* Gauge first derivatives */
  double d_Alpha[3];
  d_Alpha[0] = dupt[U_ALPHA*n_dim + 0];
  d_Alpha[1] = dupt[U_ALPHA*n_dim + 1];
  d_Alpha[2] = dupt[U_ALPHA*n_dim + 2];

  double d_Betau[3][3];
  d_Betau[0][0] = dupt[U_SHIFTX*n_dim + 0];
  d_Betau[1][0] = dupt[U_SHIFTX*n_dim + 1];
  d_Betau[2][0] = dupt[U_SHIFTX*n_dim + 2];

  d_Betau[0][1] = dupt[U_SHIFTY*n_dim + 0];
  d_Betau[1][1] = dupt[U_SHIFTY*n_dim + 1];
  d_Betau[2][1] = dupt[U_SHIFTY*n_dim + 2];

  d_Betau[0][2] = dupt[U_SHIFTZ*n_dim + 0];
  d_Betau[1][2] = dupt[U_SHIFTZ*n_dim + 1];
  d_Betau[2][2] = dupt[U_SHIFTZ*n_dim + 2];

 /*
  *  2d arrays are 
          ( a00, a01, a02 )
          ( a10, a11, a12 )
          ( a20, a21, a22 )
     But for symmetric 2D tensors, we only store the unique components on the
     grid. The grid functions are arranged in a 1D array
          ( a00 , a01, a02, a11, a12, a22 )

     To map the 2D index notation to the sparse storage 1d storage, we use
     func_map. Given 2 indices, func_map returns the location of the component
     in the 1d array.

  *  Alternative:  2*i + j -i*j/4
  */

  double adv_d_gtd[3][3][3];
  double adv_d_Atd[3][3][3];
  double adv_d_Bu[3][3];
  double adv_d_Betau[3][3];
  double adv_d_Gamt[3][3];
  double adv_d_Alpha[3];
  double adv_d_chi[3];
  double adv_d_trK[3];
  double u1d[7];

  for (int d = 0; d < 3; d++) {
    for (int m = 0; m < 3; m++) {
     /* Gamt */
      for (int p = 0; p < npts[d]; p++) {
        u1d[p] = u[d][p+bi[d]][U_GAMTX+m];
      }
      adv_d_Gamt[d][m] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);

     /* Beta */
      for (int p = 0; p < npts[d]; p++) {
        u1d[p] = u[d][p+bi[d]][U_SHIFTX+m];
      }
      adv_d_Betau[d][m] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);

     /* Gauge B */
      for (int p = 0; p < npts[d]; p++) {
        u1d[p] = u[d][p+bi[d]][U_GBX+m];
      }
      adv_d_Bu[d][m] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);
    }

   /* chi */
    for (int p = 0; p < npts[d]; p++) {
      u1d[p] = u[d][p+bi[d]][U_CHI];
    }
    adv_d_chi[d] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);

   /* trK */
    for (int p = 0; p < npts[d]; p++) {
      u1d[p] = u[d][p+bi[d]][U_TRK];
    }
    adv_d_trK[d] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);

   /* Alpha */
    for (int p = 0; p < npts[d]; p++) {
      u1d[p] = u[d][p+bi[d]][U_ALPHA];
    }
    adv_d_Alpha[d] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);
 
    for (int j = 0; j < 3; j++) {
      for (int k = j; k < 3; k++) {
       /* gtd */
        for (int p = 0; p < npts[d]; p++) {
          u1d[p] = u[d][p+bi[d]][U_GTXX + 2*j+k-j*k/4];
        }
        adv_d_gtd[d][j][k] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);

       /* Atd */
        for (int p = 0; p < npts[d]; p++) {
          u1d[p] = u[d][p+bi[d]][U_ATXX + 2*j+k-j*k/4];
        }
        adv_d_Atd[d][j][k] = deriv42adv_x(u1d, dx[d], Betau[d], indx[d], npts[d]);
      }
    }
  }



  // second derivatives -- xx, yy, and zz second derivatives
  double xd_chi[3];

  for (int d = x_dir; d < 3; d++) {

    for (int i = 0; i < 3; i++) {  
      for (int j = i; j < 3; j++) {  
        for (int p = 0; p < npts[d]; p++) {
          u1d[p] = u[d][p+bi[d]][U_GTXX + 2*i + j - i*j/4];
        }
        dd_gtd[d][d][i][j] = deriv42_xx(u1d, dx[d], indx[d], npts[d]);
        dd_gtd[d][d][j][i] = dd_gtd[d][d][i][j];
      }
    }

    for (int m = 0; m < 3; m++) {
      for (int p = 0; p < npts[d]; p++) {
        u1d[p] = u[d][p+bi[d]][U_SHIFTX+m];
      }
      dd_Betau[d][d][m] = deriv42_xx(u1d, dx[d], indx[d], npts[d]);
    }

    for (int p = 0; p < npts[d]; p++) {
      u1d[p] = u[d][p+bi[d]][U_CHI];
    }
    dd_chi[d][d] = deriv42_xx(u1d, dx[d], indx[d], npts[d]);
    xd_chi[d] = deriv42_x(u1d, dx[d], indx[d], npts[d]);

    for (int p = 0; p < npts[d]; p++) {
      u1d[p] = u[d][p+bi[d]][U_ALPHA];
    }
    dd_Alpha[d][d] = deriv42_xx(u1d, dx[d], indx[d], npts[d]);

/*
    if ( point->coords.pos[0] > -7.21 &&  point->coords.pos[0] < -7.19 &&
       point->coords.pos[1] > -7.21 &&  point->coords.pos[1] < -7.19 &&
       point->coords.pos[2] > -7.21 &&  point->coords.pos[2] < -7.19) {

      printf("d = %d",d);
    }
*/
 
    double rdchi = xd_chi[d] - d_chi[d];
    //if (fabs(rdchi) > 1.0e-5) {
     // printf("  >>> DCHI d=%d:  rdchi=%g, indx=%d, index=%d, d_chi=%g, xd_chi=%g, \n",d,rdchi,indx[d], point->index.idx[d], d_chi[d], xd_chi[d]);
    //}
  }

  // mixed second derivatives -- xy, xz, and yz second derivatives
  for (int m = 0; m < 2; m++) {
    for (int p = m+1; p < 3; p++) {
      int k;
      for (int i = 0; i < 3; i++) {
        for (int j = i; j < 3; j++) {
          k = (U_GTXX + 2*i + j - i*j/4)*n_dim + m;

          for (int q = 0; q < npts[p]; q++) {
            u1d[q] = du[p][q+bi[p]][k];
          }
          dd_gtd[p][m][i][j] = deriv42_x(u1d, dx[p], indx[p], npts[p]);
          dd_gtd[m][p][i][j] = dd_gtd[p][m][i][j];
          dd_gtd[m][p][j][i] = dd_gtd[m][p][i][j];
          dd_gtd[p][m][j][i] = dd_gtd[m][p][i][j];
        }
      }

      for (int i = 0; i < 3; i++) {
        k = (U_SHIFTX + i)*n_dim + m;
        for (int q = 0; q < npts[p]; q++) {
          u1d[q] = du[p][q+bi[p]][k];
        }
        dd_Betau[p][m][i] = deriv42_x(u1d, dx[p], indx[p], npts[p]); 
        dd_Betau[m][p][i] = dd_Betau[p][m][i];
      }

      k = U_ALPHA * n_dim + m;
      for (int q = 0; q < npts[p]; q++) {
        u1d[q] = du[p][q+bi[p]][k];
      }
      dd_Alpha[p][m] = deriv42_x(u1d, dx[p], indx[p], npts[p]);
      dd_Alpha[m][p] = dd_Alpha[p][m];

      k = U_CHI * n_dim + m;
      for (int q = 0; q < npts[p]; q++) {
        u1d[q] = du[p][q+bi[p]][k];
      }
      dd_chi[p][m] = deriv42_x(u1d, dx[p], indx[p], npts[p]);
      dd_chi[m][p] = dd_chi[p][m];

    }
  }


#include "include/bssn_vars.h"

#include "include/bssn_rhs.h"

#if 0
  if ( point->coords.pos[0] > -7.21 &&  point->coords.pos[0] < -7.19 &&
       point->coords.pos[1] > -7.21 &&  point->coords.pos[1] < -7.19 &&
       point->coords.pos[2] > -7.21 &&  point->coords.pos[2] < -7.19) {

    printf(" x=%g, y=%g, z=%g\n",point->coords.pos[0],point->coords.pos[1],point->coords.pos[2]); 
    printf("gtu[0][0]     = %g\n",gtu[0][0]);
    printf("gtu[1][1]     = %g\n",gtu[1][1]);
    printf("d_Gamt[0][0]  = %g\n",d_Gamt[0][0]);
    printf("d_Gamt[1][0]  = %g\n",d_Gamt[1][0]);
    printf("d_Gamt[1][1]  = %g\n",d_Gamt[1][1]);
    printf("d_Gamt[2][2]  = %g\n",d_Gamt[2][2]);
    printf("Atu[0][0]     = %g\n",Atu[0][0]);
    printf("Atu[1][1]     = %g\n",Atu[1][1]);
    printf("Atu[2][2]     = %g\n",Atu[2][2]);
    printf("Atd_rhs[0][0] = %g\n",Atd_rhs[0][0]);
    printf("Atd_rhs[0][1] = %g\n",Atd_rhs[0][1]);
    printf("Atd_rhs[1][1] = %g\n",Atd_rhs[1][1]);
    printf("Atd_rhs[2][2] = %g\n",Atd_rhs[2][2]);
    printf("Gamt_rhs[0]   = %g\n",Gamt_rhs[0]);
    printf("Gamt_rhs[1]   = %g\n",Gamt_rhs[1]);
    printf("Gamt_rhs[2]   = %g\n",Gamt_rhs[2]);
    exit(1);
  }
#endif

#include "include/bssn_adv_rhs.h"

  rhs[U_ALPHA] = Alpha_rhs;
  rhs[U_SHIFTX] = Betau_rhs[0];
  rhs[U_SHIFTY] = Betau_rhs[1];
  rhs[U_SHIFTZ] = Betau_rhs[2];
  rhs[U_GBX] = Bu_rhs[0];
  rhs[U_GBY] = Bu_rhs[1];
  rhs[U_GBZ] = Bu_rhs[2];

  rhs[U_CHI] = chi_rhs;
  rhs[U_TRK] = trK_rhs;
  rhs[U_GTXX] = gtd_rhs[0][0];
  rhs[U_GTXY] = gtd_rhs[0][1];
  rhs[U_GTXZ] = gtd_rhs[0][2];
  rhs[U_GTYY] = gtd_rhs[1][1];
  rhs[U_GTYZ] = gtd_rhs[1][2];
  rhs[U_GTZZ] = gtd_rhs[2][2];

  rhs[U_ATXX] = Atd_rhs[0][0];
  rhs[U_ATXY] = Atd_rhs[0][1];
  rhs[U_ATXZ] = Atd_rhs[0][2];
  rhs[U_ATYY] = Atd_rhs[1][1];
  rhs[U_ATYZ] = Atd_rhs[1][2];
  rhs[U_ATZZ] = Atd_rhs[2][2];

  rhs[U_GAMTX] = Gamt_rhs[0];
  rhs[U_GAMTY] = Gamt_rhs[1];
  rhs[U_GAMTZ] = Gamt_rhs[2];

#if 0
  if ( point->coords.pos[0] > -7.21 &&  point->coords.pos[0] < -7.19 &&
       point->coords.pos[1] > -7.21 &&  point->coords.pos[1] < -7.19 &&
       point->coords.pos[2] > -7.21 &&  point->coords.pos[2] < -7.19) {

  printf(" x=%g, y=%g, z=%g\n",point->coords.pos[0],point->coords.pos[1],point->coords.pos[2]); 
    printf("gtu[0][0]     = %g\n",gtu[0][0]);
    printf("gtu[1][1]     = %g\n",gtu[1][1]);
    printf("d_Gamt[0][0]  = %g\n",d_Gamt[0][0]);
    printf("d_Gamt[1][0]  = %g\n",d_Gamt[1][0]);
    printf("d_Gamt[1][1]  = %g\n",d_Gamt[1][1]);
    printf("d_Gamt[2][2]  = %g\n",d_Gamt[2][2]);
    printf("Atu[0][0]     = %g\n",Atu[0][0]);
    printf("Atu[1][1]     = %g\n",Atu[1][1]);
    printf("Atu[2][2]     = %g\n",Atu[2][2]);
    printf("Atd_rhs[0][0] = %g\n",Atd_rhs[0][0]);
    printf("Atd_rhs[0][1] = %g\n",Atd_rhs[0][1]);
    printf("Atd_rhs[1][1] = %g\n",Atd_rhs[1][1]);
    printf("Atd_rhs[2][2] = %g\n",Atd_rhs[2][2]);
    printf("Gamt_rhs[0]   = %g\n",Gamt_rhs[0]);
    printf("Gamt_rhs[1]   = %g\n",Gamt_rhs[1]);
    printf("Gamt_rhs[2]   = %g\n",Gamt_rhs[2]);
    exit(1);
  }
#endif
}
// }}}

// sommerfeld_boundary {{{
/*---------------------------------------------------------------------------*
 *
 *
 *
 *---------------------------------------------------------------------------*/
void sommerfeld_boundary(coll_point_t *point, const int gen)
{


  const double x = point->coords.pos[0];
  const double y = point->coords.pos[1];
  const double z = point->coords.pos[2];
  const double inv_r = 1.0 / sqrt(x*x + y*y + z*z);

  double *rhs = point->rhs[gen];
  double *u = point->u[gen];
  double *du = point->du[gen];

  const double f_falloff[24] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 
                                 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                                 2.0, 2.0, 2.0, 2.0, 2.0, 2.0 };
  const double f_asymptotic[24] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
                                    0.0, 0.0, 0.0, 
                                    0.0, 0.0, 0.0, 
                                    1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };

  for (int i = 0; i < n_variab; i++) {
    rhs[i] = - inv_r * (
                         x * du[i*n_dim + 0]
                       + y * du[i*n_dim + 1]
                       + z * du[i*n_dim + 2]
                       + f_falloff[i]
                           * (   u[i]
                               - f_asymptotic[i] )
                             );
    
  }

#if 0
  /* Some debugging stuff */
  if (y < -11.99 && y > -12.01 && x < -4.3 && x > -4.33 
                               && z > 3.82 && z < 3.89) {
    printf(" x = %g\n",x);
    printf(" y = %g\n",x);
    printf(" z = %g\n",x);
    //printf(" x = %g, y = %g, z = %g\n",x,y,z);
  }
#endif

}
// }}}

// enforce_bssn_constraints {{{
/*---------------------------------------------------------------------------*
 *
 *
 *
 *---------------------------------------------------------------------------*/
void enforce_bssn_constraints(double *upt)
{

  double gtd[3][3], Atd[3][3], gtu[3][3];
  double chi = upt[U_CHI];
  const double one_third = 1.0/3.0;

  gtd[0][0] = upt[U_GTXX];
  gtd[0][1] = upt[U_GTXY];
  gtd[0][2] = upt[U_GTXZ];
  gtd[1][0] = upt[U_GTXY];
  gtd[1][1] = upt[U_GTYY];
  gtd[1][2] = upt[U_GTYZ];
  gtd[2][0] = upt[U_GTXZ];
  gtd[2][1] = upt[U_GTYZ];
  gtd[2][2] = upt[U_GTZZ];

  Atd[0][0] = upt[U_ATXX];
  Atd[0][1] = upt[U_ATXY];
  Atd[0][2] = upt[U_ATXZ];
  Atd[1][0] = upt[U_ATXY];
  Atd[1][1] = upt[U_ATYY];
  Atd[1][2] = upt[U_ATYZ];
  Atd[2][0] = upt[U_ATXZ];
  Atd[2][1] = upt[U_ATYZ];
  Atd[2][2] = upt[U_ATZZ];


 /* Require gtd to have unit determinant */
  double det_gtd =   gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
            - gtd[0][1]*gtd[0][1]*gtd[2][2]
            + 2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
            - gtd[0][2]*gtd[0][2]*gtd[1][1];

  double det_gtd_to_third = pow(det_gtd,one_third);
  if (det_gtd_to_third < 0.0) {   /* FIXME is this needed? */
    det_gtd_to_third = 1.0;
  }
  double det_gtd_to_neg_third = 1.0 / det_gtd_to_third;

  gtd[0][0] = det_gtd_to_neg_third * gtd[0][0];
  gtd[0][1] = det_gtd_to_neg_third * gtd[0][1];
  gtd[0][2] = det_gtd_to_neg_third * gtd[0][2];
  gtd[1][1] = det_gtd_to_neg_third * gtd[1][1];
  gtd[1][2] = det_gtd_to_neg_third * gtd[1][2];
  gtd[2][2] = det_gtd_to_neg_third * gtd[2][2];

  det_gtd =   gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
            - gtd[0][1]*gtd[0][1]*gtd[2][2]
            + 2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
            - gtd[0][2]*gtd[0][2]*gtd[1][1];



  double detgt_m1 = det_gtd - 1.0;

  if (fabs(detgt_m1) > 1.0e-6) {
    printf("enforce_bssn_constraint: det(gtd) != 1. det=%g\n",det_gtd);
    printf("      gtd(1,1)=%g\n",gtd[0][0]);
    printf("      gtd(1,2)=%g\n",gtd[0][1]);
    printf("      gtd(1,3)=%g\n",gtd[0][2]);
    printf("      gtd(2,2)=%g\n",gtd[1][1]);
    printf("      gtd(2,3)=%g\n",gtd[1][2]);
    printf("      gtd(3,3)=%g\n",gtd[2][2]);
  }

  double idet_gtd = 1.0/det_gtd;
  gtu[0][0] = idet_gtd*(gtd[1][1]*gtd[2][2]-gtd[1][2]*gtd[1][2]);
  gtu[0][1] = idet_gtd*(-gtd[0][1]*gtd[2][2]+gtd[0][2]*gtd[1][2]);
  gtu[0][2] = idet_gtd*(gtd[0][1]*gtd[1][2]-gtd[0][2]*gtd[1][1]);
  gtu[1][0] = gtu[0][1];
  gtu[1][1] = idet_gtd*(gtd[0][0]*gtd[2][2]-gtd[0][2]*gtd[0][2]);
  gtu[1][2] = idet_gtd*(-gtd[0][0]*gtd[1][2]+gtd[0][1]*gtd[0][2]);
  gtu[2][0] = gtu[0][2];
  gtu[2][1] = gtu[1][2];
  gtu[2][2] = idet_gtd*(gtd[0][0]*gtd[1][1]-gtd[0][1]*gtd[0][1]);

 /* Require Atd to be traceless. */
  double trace_Atd =    Atd[0][0]*gtu[0][0]
                      + Atd[1][1]*gtu[1][1]
                      + Atd[2][2]*gtu[2][2]
                      + 2.0 * (   Atd[0][1]*gtu[0][1]
                                + Atd[0][2]*gtu[0][2]
                                + Atd[1][2]*gtu[1][2]  );

  double neg_one_third_trace_Atd = - one_third * trace_Atd;

  Atd[0][0] = Atd[0][0] + neg_one_third_trace_Atd * gtd[0][0];
  Atd[0][1] = Atd[0][1] + neg_one_third_trace_Atd * gtd[0][1];
  Atd[0][2] = Atd[0][2] + neg_one_third_trace_Atd * gtd[0][2];
  Atd[1][1] = Atd[1][1] + neg_one_third_trace_Atd * gtd[1][1];
  Atd[1][2] = Atd[1][2] + neg_one_third_trace_Atd * gtd[1][2];
  Atd[2][2] = Atd[2][2] + neg_one_third_trace_Atd * gtd[2][2];

  double tr_A =    Atd[0][0]*gtu[0][0]
                 + Atd[1][1]*gtu[1][1]
                 + Atd[2][2]*gtu[2][2]
                 + 2.0 * (   Atd[0][1]*gtu[0][1]
                           + Atd[0][2]*gtu[0][2]
                           + Atd[1][2]*gtu[1][2]  );


  //assert(tr_A < 1.0e-8);

  if (fabs(tr_A) > 1.0e-6) {
    printf("enforce_bssn_constraint: tr_A != 0. tr_A=%g\n",tr_A);
    printf("      Atd(1,1)=%g\n",Atd[0][0]);
    printf("      Atd(1,2)=%g\n",Atd[0][1]);
    printf("      Atd(1,3)=%g\n",Atd[0][2]);
    printf("      Atd(2,2)=%g\n",Atd[1][1]);
    printf("      Atd(2,3)=%g\n",Atd[1][2]);
    printf("      Atd(3,3)=%g\n",Atd[2][2]);
  }


  upt[U_ATXX] = Atd[0][0];
  upt[U_ATXY] = Atd[0][1];
  upt[U_ATXZ] = Atd[0][2];
  upt[U_ATYY] = Atd[1][1];
  upt[U_ATYZ] = Atd[1][2];
  upt[U_ATZZ] = Atd[2][2];

  upt[U_GTXX] = gtd[0][0];
  upt[U_GTXY] = gtd[0][1];
  upt[U_GTXZ] = gtd[0][2];
  upt[U_GTYY] = gtd[1][1];
  upt[U_GTYZ] = gtd[1][2];
  upt[U_GTZZ] = gtd[2][2];

  if ( chi < pars.chi_floor ) {
   /* FIXME This needs to be fixed when we add a fluid to the code. */
   /* ! First rescale the densitized fluid variables.
      ! The include file bssn_puncture_fluid_rescale.inc
      ! must be provided in the BSSN_*MHD project.

      ! Chi must be positive to do the rescaling of fluid variables.
      if ( chi <= 0.0) {
        chi = pars.chi_floor;
      }
      else {
        // ok... go ahead and rescale the fluid variables.
      }


    */

   /* now place the floor on chi */
    upt[U_CHI] = pars.chi_floor;
  }

  if ( upt[U_ALPHA] < pars.chi_floor ) {
    /* now place the floor on the lapse */
    upt[U_ALPHA] = pars.chi_floor;
  }

}
// }}}

