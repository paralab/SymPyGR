#ifndef _HAVE_TESTCODE_H
#define _HAVE_TESTCODE_H

#define AVX_SIMD_LENGTH 4


#define restrict __restrict__

void define_coordinates(double *x, double xmin, double xmax, double &h, int n);
void initial_data(double *u[], double *xi[], int shp[]);
void set_zero(double *u, int shp[]);
double L2norm(double * const u, const int shp[]);
void cal_norms(double df[], double *f[], int shp[], int bounds);
void cal_diffs(double *df[], double *f1[], double *f2[], int shp[], int bounds);

double grad(int dir, double *u, double h, int idx[], int shp[]);
double grad2(int dir1, int dir2, double *u, double h, int idx[], int shp[]);
double agrad(int dir1, double beta, double *u, double h, int idx[], int shp[]);

void compute_derivatives(const double* u, const int *shp, double h, const double *beta, double** dd);

void sdfoutput(char *fname, int *shp, double *coords, double *data);
void printdiff(double **dtu, double **dth, int *shp);

void transpose_xy(double * const restrict fout,
                  const double * const restrict fin,
                  const int nx, const int ny, const int nz);
void transpose_xz(double * const restrict fout,
                  const double * const restrict fin,
                  const int nx, const int ny, const int nz);
void deriv42_x(double * const Dxu, const double * const u, const double dx,
               const int nx, const int ny, const int nz);
void deriv42_y(double * const Dyu, const double * const u, const double dy,
               const int nx, const int ny, const int nz);
void deriv42_z(double * const Dzu, const double * const u, const double dz,
               const int nx, const int ny, const int nz);
void deriv42adv_x(double * const Dxu, const double * const u, const double dx,
                  const int nx, const int ny, const int nz,
                  const double * const betax);
void deriv42adv_y(double * const restrict Dyu,
                  const double * const restrict u,
                  const double dy,
                  const int nx, const int ny, const int nz,
                  const double * const betay);
void deriv42adv_z(double * const restrict Dzu,
                  const double * const restrict u,
                  const double dz,
                  const int nx, const int ny, const int nz,
                  const double * const betaz);
void deriv42_xx(double * const DxDxu, const double * const u,
                const double dx, const int nx, const int ny, const int nz);
void deriv42_yy(double * const restrict Du,
                const double * const restrict u,
                const double dy, const int nx, const int ny, const int nz);
void deriv42_zz(double * const restrict Du,
                const double * const restrict u,
                const double dz, const int nx, const int ny, const int nz);

int apply_stencil_x(const double* in, double* out, int *sz, double h);
int apply_stencil_y(const double* in, double* out, int *sz, double h);
int apply_stencil_z(const double* in, double* out, int *sz, double h);
int avx_apply_stencil_x(const double* in, double* out, int *sz, double h);
int avx_apply_stencil_y(const double* in, double* out, int *sz, double h);
int avx_apply_stencil_z(const double* in, double* out, int *sz, double h);
int avx_apply_stencil_d_ddx(const double* in, double* const out,
                          double *const ddout, int *sz, double h);
int avx_apply_stencil_d_ddy(const double* in, double* const out,
                          double *const ddout, int *sz, double h);
int avx_apply_stencil_d_ddz(const double* in, double* const out,
                           double* const ddout, int *sz, double h);
int avx_apply_adv_stencil_x(const double* in, double* out, double* beta,
                            int *sz, double h);
int avx_apply_adv_stencil_y(const double* in, double* out, double* beta,
                            int *sz, double h);
int avx_apply_adv_stencil_z(const double* in, double* out, double* beta,
                            int *sz, double h);



void rhs(double *a_rhs_3D, double *b_rhs0_3D, double *b_rhs1_3D,
         double *b_rhs2_3D, double *chi_rhs_3D, double *K_rhs_3D,
         double *gt_rhs00_3D, double *gt_rhs01_3D, double *gt_rhs02_3D,
         double *gt_rhs11_3D, double *gt_rhs12_3D, double *gt_rhs22_3D,
         double *At_rhs00_3D, double *At_rhs01_3D, double *At_rhs02_3D,
         double *At_rhs11_3D, double *At_rhs12_3D, double *At_rhs22_3D,
         double *Gt_rhs0_3D, double *Gt_rhs1_3D, double *Gt_rhs2_3D,
         double *B_rhs0_3D, double *B_rhs1_3D, double *B_rhs2_3D,
         double *alpha_3D, double *beta0_3D, double *beta1_3D,
         double *beta2_3D, double *chi_3D, double *K_3D,
         double *gt0_3D, double *gt1_3D, double *gt2_3D,
         double *gt3_3D, double *gt4_3D, double *gt5_3D,
         double *At0_3D, double *At1_3D, double *At2_3D,
         double *At3_3D, double *At4_3D, double *At5_3D,
         double *Gt0_3D, double *Gt1_3D, double *Gt2_3D,
         double *B0_3D, double *B1_3D, double *B2_3D,
         int *lambda, int *lambda_f,
         double eta_const, int eta_damping, int eta_damping_exp, double eta_R0,
         double trK0, double chi_floor,
         double *xi[], int shp[]);

void rhs_opt(double *a_rhs_3D, double *b_rhs0_3D, double *b_rhs1_3D,
         double *b_rhs2_3D, double *chi_rhs_3D, double *K_rhs_3D,
         double *gt_rhs00_3D, double *gt_rhs01_3D, double *gt_rhs02_3D,
         double *gt_rhs11_3D, double *gt_rhs12_3D, double *gt_rhs22_3D,
         double *At_rhs00_3D, double *At_rhs01_3D, double *At_rhs02_3D,
         double *At_rhs11_3D, double *At_rhs12_3D, double *At_rhs22_3D,
         double *Gt_rhs0_3D, double *Gt_rhs1_3D, double *Gt_rhs2_3D,
         double *B_rhs0_3D, double *B_rhs1_3D, double *B_rhs2_3D,
         double *alpha_3D, double *beta0_3D, double *beta1_3D,
         double *beta2_3D, double *chi_3D, double *K_3D,
         double *gt0_3D, double *gt1_3D, double *gt2_3D,
         double *gt3_3D, double *gt4_3D, double *gt5_3D,
         double *At0_3D, double *At1_3D, double *At2_3D,
         double *At3_3D, double *At4_3D, double *At5_3D,
         double *Gt0_3D, double *Gt1_3D, double *Gt2_3D,
         double *B0_3D, double *B1_3D, double *B2_3D,
         int *lambda, int *lambda_f,
         double eta_const, int eta_damping, int eta_damping_exp, double eta_R0,
         double trK0, double chi_floor,
         double *xi[], int shp[]);

void rhs_vec(double *a_rhs_3D, double *b_rhs0_3D, double *b_rhs1_3D,
         double *b_rhs2_3D, double *chi_rhs_3D, double *K_rhs_3D,
         double *gt_rhs00_3D, double *gt_rhs01_3D, double *gt_rhs02_3D,
         double *gt_rhs11_3D, double *gt_rhs12_3D, double *gt_rhs22_3D,
         double *At_rhs00_3D, double *At_rhs01_3D, double *At_rhs02_3D,
         double *At_rhs11_3D, double *At_rhs12_3D, double *At_rhs22_3D,
         double *Gt_rhs0_3D, double *Gt_rhs1_3D, double *Gt_rhs2_3D,
         double *B_rhs0_3D, double *B_rhs1_3D, double *B_rhs2_3D,
         double *alpha_3D, double *beta0_3D, double *beta1_3D,
         double *beta2_3D, double *chi_3D, double *K_3D,
         double *gt0_3D, double *gt1_3D, double *gt2_3D,
         double *gt3_3D, double *gt4_3D, double *gt5_3D,
         double *At0_3D, double *At1_3D, double *At2_3D,
         double *At3_3D, double *At4_3D, double *At5_3D,
         double *Gt0_3D, double *Gt1_3D, double *Gt2_3D,
         double *B0_3D, double *B1_3D, double *B2_3D,
         int *lambda, int *lambda_f,
         double eta_const, int eta_damping, int eta_damping_exp, double eta_R0,
         double trK0, double chi_floor,
         double *xi[], int shp[]);

void hadrhs( double *  dt_alpha, double *  dt_shiftx,
             double *  dt_shifty, double *  dt_shiftz,
             double *  dt_chi, double *  dt_trK,
             double *  dt_gtxx, double *  dt_gtxy,
             double *  dt_gtxz, double *  dt_gtyy,
             double *  dt_gtyz, double *  dt_gtzz,
             double *  dt_Atxx, double *  dt_Atxy,
             double *  dt_Atxz, double *  dt_Atyy,
             double *  dt_Atyz, double *  dt_Atzz,
             double *  dt_Gamtx, double *  dt_Gamty,
             double *  dt_Gamtz, double *  dt_gbx,
             double *  dt_gby, double *  dt_gbz,
             double *  alpha,  double *  shiftx,
             double *  shifty,  double *  shiftz,
             double *  chi,  double *  trK,
             double *  gtxx,  double *  gtxy,
             double *  gtxz,  double *  gtyy,
             double *  gtyz,  double *  gtzz,
             double *  Atxx,  double *  Atxy,
             double *  Atxz,  double *  Atyy,
             double *  Atyz,  double *  Atzz,
             double *  Gamtx,  double *  Gamty,
             double *  Gamtz,  double *  gbx,
             double *  gby,  double *  gbz,
             int *  lambda,  int *  lambda_f,
             double eta_const,  int eta_damping,
             int eta_damping_exp,  double eta_R0,
             double trK0,  double chi_floor,
             double *  xi[],  int shp[] );


void hadrhs_transpose( double *  dt_alpha, double *  dt_shiftx,
             double *  dt_shifty, double *  dt_shiftz,
             double *  dt_chi, double *  dt_trK,
             double *  dt_gtxx, double *  dt_gtxy,
             double *  dt_gtxz, double *  dt_gtyy,
             double *  dt_gtyz, double *  dt_gtzz,
             double *  dt_Atxx, double *  dt_Atxy,
             double *  dt_Atxz, double *  dt_Atyy,
             double *  dt_Atyz, double *  dt_Atzz,
             double *  dt_Gamtx, double *  dt_Gamty,
             double *  dt_Gamtz, double *  dt_gbx,
             double *  dt_gby, double *  dt_gbz,
              double *  alpha,  double *  shiftx,
              double *  shifty,  double *  shiftz,
              double *  chi,  double *  trK,
              double *  gtxx,  double *  gtxy,
              double *  gtxz,  double *  gtyy,
              double *  gtyz,  double *  gtzz,
              double *  Atxx,  double *  Atxy,
              double *  Atxz,  double *  Atyy,
              double *  Atyz,  double *  Atzz,
              double *  Gamtx,  double *  Gamty,
              double *  Gamtz,  double *  gbx,
              double *  gby,  double *  gbz,
              int *  lambda,  int *  lambda_f,
              double eta_,  int eta_damping,
              int eta_damping_exp,  double eta_R0,
              double trK0,  double chi_floor,
              double *  xi[],  int shp[] );


extern "C"
void cal_bssn_rhs_(
       double *dt_Alpha3D, double *dt_shiftx, double *dt_shifty,
       double *dt_shiftz, double *dt_chi3D,
       double *dt_trK3D, double *dt_gtxx, double *dt_gtxy, double *dt_gtxz,
       double *dt_gtyy, double *dt_gtyz, double *dt_gtzz,
       double *dt_Atxx, double *dt_Atxy, double *dt_Atxz,
       double *dt_Atyy, double *dt_Atyz, double *dt_Atzz,
       double *dt_Gamtx, double *dt_Gamty, double *dt_Gamtz,
       double *dt_gbx, double *dt_gby, double *dt_gbz,
       double *Alpha3D, double *shiftx, double *shifty, double *shiftz,
       double *chi3D, double *trK3D,
       double *gtxx, double *gtxy, double *gtxz,
       double *gtyy, double *gtyz, double *gtzz,
       double *Atxx, double *Atxy, double *Atxz, double *Atyy, double *Atyz,
       double *Atzz, double *Gamtx, double *Gamty, double *Gamtz,
       double *gbx, double *gby, double *gbz,
       int *lambda, int *lambda_f,
       double *eta, int *eta_damping, int *eta_damping_exp, double *eta_R0,
       double *trk0, double *chi_floor,
       double *x1d, double *y1d, double *z1d, int *nx, int *ny, int *nz );

extern "C"
void cal_had_bssn_rhs_(
       double * dt_Alpha3D,
       double * dt_shiftx, double * dt_shifty, double * dt_shiftz,
       double * dt_chi3D, double * dt_trK3D,
       double * dt_gtxx, double * dt_gtxy, double * dt_gtxz,
       double * dt_gtyy, double * dt_gtyz, double * dt_gtzz,
       double * dt_Atxx, double * dt_Atxy, double * dt_Atxz,
       double * dt_Atyy, double * dt_Atyz, double * dt_Atzz,
       double * dt_Gamtx, double * dt_Gamty, double * dt_Gamtz,
       double * dt_gbx, double * dt_gby, double * dt_gbz,
       double * Alpha3D, double * shiftx, double * shifty, double * shiftz,
       double * chi3D, double * trK3D,
       double * gtxx, double * gtxy, double * gtxz,
       double * gtyy, double * gtyz, double * gtzz,
       double * Atxx, double * Atxy, double * Atxz,
       double * Atyy, double * Atyz, double * Atzz,
       double * Gamtx, double * Gamty, double * Gamtz,
       double * gbx, double * gby, double * gbz,
       double *dx_alpha,double *dy_alpha,double *dz_alpha,double *dx_shiftx,
       double *dy_shiftx,double *dz_shiftx,double *dx_shifty,double *dy_shifty,
       double *dz_shifty,double *dx_shiftz,double *dy_shiftz,double *dz_shiftz,
       double *dx_gbx,double *dy_gbx,double *dz_gbx,double *dx_gby,
       double *dy_gby,double *dz_gby,double *dx_gbz,double *dy_gbz,
       double *dz_gbz,double *dx_chi,double *dy_chi,double *dz_chi,
       double *dx_Gamtx,double *dy_Gamtx,double *dz_Gamtx,double *dx_Gamty,
       double *dy_Gamty,double *dz_Gamty,double *dx_Gamtz,double *dy_Gamtz,
       double *dz_Gamtz,double *dx_trK,double *dy_trK,double *dz_trK,
       double *dx_gtxx,double *dy_gtxx,double *dz_gtxx,double *dx_gtxy,
       double *dy_gtxy,double *dz_gtxy,double *dx_gtxz,double *dy_gtxz,
       double *dz_gtxz,double *dx_gtyy,double *dy_gtyy,double *dz_gtyy,
       double *dx_gtyz,double *dy_gtyz,double *dz_gtyz,double *dx_gtzz,
       double *dy_gtzz,double *dz_gtzz,double *dx_Atxx,double *dy_Atxx,
       double *dz_Atxx,double *dx_Atxy,double *dy_Atxy,double *dz_Atxy,
       double *dx_Atxz,double *dy_Atxz,double *dz_Atxz,double *dx_Atyy,
       double *dy_Atyy,double *dz_Atyy,double *dx_Atyz,double *dy_Atyz,
       double *dz_Atyz,double *dx_Atzz,double *dy_Atzz,double *dz_Atzz,
       double *dxx_gtxx,double *dxy_gtxx,double *dxz_gtxx,double *dyy_gtxx,
       double *dyz_gtxx,double *dzz_gtxx,double *dxx_gtxy,double *dxy_gtxy,
       double *dxz_gtxy,double *dyy_gtxy,double *dyz_gtxy,double *dzz_gtxy,
       double *dxx_gtxz,double *dxy_gtxz,double *dxz_gtxz,double *dyy_gtxz,
       double *dyz_gtxz,double *dzz_gtxz,double *dxx_gtyy,double *dxy_gtyy,
       double *dxz_gtyy,double *dyy_gtyy,double *dyz_gtyy,double *dzz_gtyy,
       double *dxx_gtyz,double *dxy_gtyz,double *dxz_gtyz,double *dyy_gtyz,
       double *dyz_gtyz,double *dzz_gtyz,double *dxx_gtzz,double *dxy_gtzz,
       double *dxz_gtzz,double *dyy_gtzz,double *dyz_gtzz,double *dzz_gtzz,
       double *dxx_chi,double *dxy_chi,double *dxz_chi,double *dyy_chi,
       double *dyz_chi,double *dzz_chi,double *dxx_alpha,double *dxy_alpha,
       double *dxz_alpha,double *dyy_alpha,double *dyz_alpha,double *dzz_alpha,
       double *dxx_shiftx,double *dxy_shiftx,double *dxz_shiftx,
       double *dyy_shiftx,
       double *dyz_shiftx,double *dzz_shiftx,double *dxx_shifty,
       double *dxy_shifty,
       double *dxz_shifty,double *dyy_shifty,double *dyz_shifty,
       double *dzz_shifty,
       double *dxx_shiftz,double *dxy_shiftz,double *dxz_shiftz,
       double *dyy_shiftz,
       double *dyz_shiftz,double *dzz_shiftz,
       int *lambda, int *lambda_f,
       double *eta, int *eta_damping, int *eta_damping_exp, double *R_0,
       double *trk0, double *chi_floor,
       double *x1d, double *y1d, double *z1d, int *nx, int *ny, int *nz);

extern "C"
void cal_had_bssn_rhs_lie_(
       double * dt_Alpha3D,
       double * dt_shiftx, double * dt_shifty, double * dt_shiftz,
       double * dt_chi3D, double * dt_trK3D,
       double * dt_gtxx, double * dt_gtxy, double * dt_gtxz,
       double * dt_gtyy, double * dt_gtyz,
       double * dt_gtzz, double * dt_Atxx, double * dt_Atxy,
       double * dt_Atxz, double * dt_Atyy, double * dt_Atyz,
       double * dt_Atzz, double * dt_Gamtx, double * dt_Gamty,
       double * dt_Gamtz, double * dt_gbx, double * dt_gby, double * dt_gbz,
       double *shiftx, double *shifty, double *shiftz,
       double * adx_gtxx,double * ady_gtxx,double * adz_gtxx,double * adx_gtxy,
       double * ady_gtxy,double * adz_gtxy,double * adx_gtxz,double * ady_gtxz,
       double * adz_gtxz,double * adx_gtyy,double * ady_gtyy,double * adz_gtyy,
       double * adx_gtyz,double * ady_gtyz,double * adz_gtyz,double * adx_gtzz,
       double * ady_gtzz,double * adz_gtzz,double * adx_Atxx,double * ady_Atxx,
       double * adz_Atxx,double * adx_Atxy,double * ady_Atxy,double * adz_Atxy,
       double * adx_Atxz,double * ady_Atxz,double * adz_Atxz,double * adx_Atyy,
       double * ady_Atyy,double * adz_Atyy,double * adx_Atyz,double * ady_Atyz,
       double * adz_Atyz,double * adx_Atzz,double * ady_Atzz,double * adz_Atzz,
       double * adx_alpha, double * ady_alpha,double * adz_alpha,
       double * adx_shiftx, double * ady_shiftx,double * adz_shiftx,
       double * adx_shifty,double * ady_shifty, double * adz_shifty,
       double * adx_shiftz,double * ady_shiftz,double * adz_shiftz,
       double * adx_chi,double * ady_chi,double * adz_chi,
       double * adx_Gamtx, double * ady_Gamtx,double * adz_Gamtx,
       double * adx_Gamty,double * ady_Gamty, double * adz_Gamty,
       double * adx_Gamtz,double * ady_Gamtz,double * adz_Gamtz,
       double * adx_trK,double * ady_trK,double * adz_trK,double * adx_gbx,
       double * ady_gbx,double * adz_gbx,double * adx_gby,double * ady_gby,
       double * adz_gby,double * adx_gbz,double * ady_gbz,double * adz_gbz,
       int *lambda, int *lambda_f,
       double *x1d, double *y1d, double *z1d,
       int *nx, int *ny, int *nz);

#endif
