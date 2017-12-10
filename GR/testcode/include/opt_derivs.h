  std::cout << "calling avx_apply_stencil_ddx alpha"; 
  avx_apply_stencil_d_ddx(alpha_3D, dx_alpha, dxx_alpha, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx beta0"; 
  avx_apply_stencil_d_ddx(beta0_3D, dx_beta0, dxx_beta0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx beta1"; 
  avx_apply_stencil_d_ddx(beta1_3D, dx_beta1, dxx_beta1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx beta2"; 
  avx_apply_stencil_d_ddx(beta2_3D, dx_beta2, dxx_beta2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x B0"; 
  avx_apply_stencil_x(B0_3D, dx_B0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_B0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x B1"; 
  avx_apply_stencil_x(B1_3D, dx_B1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_B1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x B2"; 
  avx_apply_stencil_x(B2_3D, dx_B2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_B2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx chi"; 
  avx_apply_stencil_d_ddx(chi_3D, dx_chi, dxx_chi, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x Gt0"; 
  avx_apply_stencil_x(Gt0_3D, dx_Gt0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_Gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x Gt1"; 
  avx_apply_stencil_x(Gt1_3D, dx_Gt1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_Gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x Gt2"; 
  avx_apply_stencil_x(Gt2_3D, dx_Gt2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_Gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x K"; 
  avx_apply_stencil_x(K_3D, dx_K, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_K[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx gt0"; 
  avx_apply_stencil_d_ddx(gt0_3D, dx_gt0, dxx_gt0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx gt1"; 
  avx_apply_stencil_d_ddx(gt1_3D, dx_gt1, dxx_gt1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx gt2"; 
  avx_apply_stencil_d_ddx(gt2_3D, dx_gt2, dxx_gt2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx gt3"; 
  avx_apply_stencil_d_ddx(gt3_3D, dx_gt3, dxx_gt3, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx gt4"; 
  avx_apply_stencil_d_ddx(gt4_3D, dx_gt4, dxx_gt4, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddx gt5"; 
  avx_apply_stencil_d_ddx(gt5_3D, dx_gt5, dxx_gt5, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dxx_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x At0"; 
  avx_apply_stencil_x(At0_3D, dx_At0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x At1"; 
  avx_apply_stencil_x(At1_3D, dx_At1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x At2"; 
  avx_apply_stencil_x(At2_3D, dx_At2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x At3"; 
  avx_apply_stencil_x(At3_3D, dx_At3, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x At4"; 
  avx_apply_stencil_x(At4_3D, dx_At4, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_x At5"; 
  avx_apply_stencil_x(At5_3D, dx_At5, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy alpha"; 
  avx_apply_stencil_d_ddy(alpha_3D, dy_alpha, dyy_alpha, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy beta0"; 
  avx_apply_stencil_d_ddy(beta0_3D, dy_beta0, dyy_beta0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy beta1"; 
  avx_apply_stencil_d_ddy(beta1_3D, dy_beta1, dyy_beta1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy beta2"; 
  avx_apply_stencil_d_ddy(beta2_3D, dy_beta2, dyy_beta2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y B0"; 
  avx_apply_stencil_y(B0_3D, dy_B0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_B0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y B1"; 
  avx_apply_stencil_y(B1_3D, dy_B1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_B1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y B2"; 
  avx_apply_stencil_y(B2_3D, dy_B2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_B2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy chi"; 
  avx_apply_stencil_d_ddy(chi_3D, dy_chi, dyy_chi, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y Gt0"; 
  avx_apply_stencil_y(Gt0_3D, dy_Gt0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_Gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y Gt1"; 
  avx_apply_stencil_y(Gt1_3D, dy_Gt1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_Gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y Gt2"; 
  avx_apply_stencil_y(Gt2_3D, dy_Gt2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_Gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y K"; 
  avx_apply_stencil_y(K_3D, dy_K, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_K[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy gt0"; 
  avx_apply_stencil_d_ddy(gt0_3D, dy_gt0, dyy_gt0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy gt1"; 
  avx_apply_stencil_d_ddy(gt1_3D, dy_gt1, dyy_gt1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy gt2"; 
  avx_apply_stencil_d_ddy(gt2_3D, dy_gt2, dyy_gt2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy gt3"; 
  avx_apply_stencil_d_ddy(gt3_3D, dy_gt3, dyy_gt3, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy gt4"; 
  avx_apply_stencil_d_ddy(gt4_3D, dy_gt4, dyy_gt4, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddy gt5"; 
  avx_apply_stencil_d_ddy(gt5_3D, dy_gt5, dyy_gt5, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dyy_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y At0"; 
  avx_apply_stencil_y(At0_3D, dy_At0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y At1"; 
  avx_apply_stencil_y(At1_3D, dy_At1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y At2"; 
  avx_apply_stencil_y(At2_3D, dy_At2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y At3"; 
  avx_apply_stencil_y(At3_3D, dy_At3, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y At4"; 
  avx_apply_stencil_y(At4_3D, dy_At4, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_y At5"; 
  avx_apply_stencil_y(At5_3D, dy_At5, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz alpha"; 
  avx_apply_stencil_d_ddz(alpha_3D, dz_alpha, dzz_alpha, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz beta0"; 
  avx_apply_stencil_d_ddz(beta0_3D, dz_beta0, dzz_beta0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz beta1"; 
  avx_apply_stencil_d_ddz(beta1_3D, dz_beta1, dzz_beta1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz beta2"; 
  avx_apply_stencil_d_ddz(beta2_3D, dz_beta2, dzz_beta2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z B0"; 
  avx_apply_stencil_z(B0_3D, dz_B0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_B0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z B1"; 
  avx_apply_stencil_z(B1_3D, dz_B1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_B1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z B2"; 
  avx_apply_stencil_z(B2_3D, dz_B2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_B2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz chi"; 
  avx_apply_stencil_d_ddz(chi_3D, dz_chi, dzz_chi, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z Gt0"; 
  avx_apply_stencil_z(Gt0_3D, dz_Gt0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_Gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z Gt1"; 
  avx_apply_stencil_z(Gt1_3D, dz_Gt1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_Gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z Gt2"; 
  avx_apply_stencil_z(Gt2_3D, dz_Gt2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_Gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z K"; 
  avx_apply_stencil_z(K_3D, dz_K, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_K[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz gt0"; 
  avx_apply_stencil_d_ddz(gt0_3D, dz_gt0, dzz_gt0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz gt1"; 
  avx_apply_stencil_d_ddz(gt1_3D, dz_gt1, dzz_gt1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz gt2"; 
  avx_apply_stencil_d_ddz(gt2_3D, dz_gt2, dzz_gt2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz gt3"; 
  avx_apply_stencil_d_ddz(gt3_3D, dz_gt3, dzz_gt3, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz gt4"; 
  avx_apply_stencil_d_ddz(gt4_3D, dz_gt4, dzz_gt4, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_ddz gt5"; 
  avx_apply_stencil_d_ddz(gt5_3D, dz_gt5, dzz_gt5, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dzz_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z At0"; 
  avx_apply_stencil_z(At0_3D, dz_At0, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z At1"; 
  avx_apply_stencil_z(At1_3D, dz_At1, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z At2"; 
  avx_apply_stencil_z(At2_3D, dz_At2, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z At3"; 
  avx_apply_stencil_z(At3_3D, dz_At3, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z At4"; 
  avx_apply_stencil_z(At4_3D, dz_At4, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_stencil_z At5"; 
  avx_apply_stencil_z(At5_3D, dz_At5, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "...calling mixed second: dxy_gt0" << endl;
  avx_apply_stencil_y(dx_gt0, dxy_gt0, sz, h);
  std::cout << "...calling mixed second: dxz_gt0" << endl;
  avx_apply_stencil_z(dx_gt0, dxz_gt0, sz, h);
  std::cout << "...calling mixed second: dyz_gt0" << endl;
  avx_apply_stencil_z(dy_gt0, dyz_gt0, sz, h);
  std::cout << "...calling mixed second: dxy_gt1" << endl;
  avx_apply_stencil_y(dx_gt1, dxy_gt1, sz, h);
  std::cout << "...calling mixed second: dxz_gt1" << endl;
  avx_apply_stencil_z(dx_gt1, dxz_gt1, sz, h);
  std::cout << "...calling mixed second: dyz_gt1" << endl;
  avx_apply_stencil_z(dy_gt1, dyz_gt1, sz, h);
  std::cout << "...calling mixed second: dxy_gt2" << endl;
  avx_apply_stencil_y(dx_gt2, dxy_gt2, sz, h);
  std::cout << "...calling mixed second: dxz_gt2" << endl;
  avx_apply_stencil_z(dx_gt2, dxz_gt2, sz, h);
  std::cout << "...calling mixed second: dyz_gt2" << endl;
  avx_apply_stencil_z(dy_gt2, dyz_gt2, sz, h);
  std::cout << "...calling mixed second: dxy_gt3" << endl;
  avx_apply_stencil_y(dx_gt3, dxy_gt3, sz, h);
  std::cout << "...calling mixed second: dxz_gt3" << endl;
  avx_apply_stencil_z(dx_gt3, dxz_gt3, sz, h);
  std::cout << "...calling mixed second: dyz_gt3" << endl;
  avx_apply_stencil_z(dy_gt3, dyz_gt3, sz, h);
  std::cout << "...calling mixed second: dxy_gt4" << endl;
  avx_apply_stencil_y(dx_gt4, dxy_gt4, sz, h);
  std::cout << "...calling mixed second: dxz_gt4" << endl;
  avx_apply_stencil_z(dx_gt4, dxz_gt4, sz, h);
  std::cout << "...calling mixed second: dyz_gt4" << endl;
  avx_apply_stencil_z(dy_gt4, dyz_gt4, sz, h);
  std::cout << "...calling mixed second: dxy_gt5" << endl;
  avx_apply_stencil_y(dx_gt5, dxy_gt5, sz, h);
  std::cout << "...calling mixed second: dxz_gt5" << endl;
  avx_apply_stencil_z(dx_gt5, dxz_gt5, sz, h);
  std::cout << "...calling mixed second: dyz_gt5" << endl;
  avx_apply_stencil_z(dy_gt5, dyz_gt5, sz, h);
  std::cout << "...calling mixed second: dxy_chi" << endl;
  avx_apply_stencil_y(dx_chi, dxy_chi, sz, h);
  std::cout << "...calling mixed second: dxz_chi" << endl;
  avx_apply_stencil_z(dx_chi, dxz_chi, sz, h);
  std::cout << "...calling mixed second: dyz_chi" << endl;
  avx_apply_stencil_z(dy_chi, dyz_chi, sz, h);
  std::cout << "...calling mixed second: dxy_alpha" << endl;
  avx_apply_stencil_y(dx_alpha, dxy_alpha, sz, h);
  std::cout << "...calling mixed second: dxz_alpha" << endl;
  avx_apply_stencil_z(dx_alpha, dxz_alpha, sz, h);
  std::cout << "...calling mixed second: dyz_alpha" << endl;
  avx_apply_stencil_z(dy_alpha, dyz_alpha, sz, h);
  std::cout << "...calling mixed second: dxy_beta0" << endl;
  avx_apply_stencil_y(dx_beta0, dxy_beta0, sz, h);
  std::cout << "...calling mixed second: dxz_beta0" << endl;
  avx_apply_stencil_z(dx_beta0, dxz_beta0, sz, h);
  std::cout << "...calling mixed second: dyz_beta0" << endl;
  avx_apply_stencil_z(dy_beta0, dyz_beta0, sz, h);
  std::cout << "...calling mixed second: dxy_beta1" << endl;
  avx_apply_stencil_y(dx_beta1, dxy_beta1, sz, h);
  std::cout << "...calling mixed second: dxz_beta1" << endl;
  avx_apply_stencil_z(dx_beta1, dxz_beta1, sz, h);
  std::cout << "...calling mixed second: dyz_beta1" << endl;
  avx_apply_stencil_z(dy_beta1, dyz_beta1, sz, h);
  std::cout << "...calling mixed second: dxy_beta2" << endl;
  avx_apply_stencil_y(dx_beta2, dxy_beta2, sz, h);
  std::cout << "...calling mixed second: dxz_beta2" << endl;
  avx_apply_stencil_z(dx_beta2, dxz_beta2, sz, h);
  std::cout << "...calling mixed second: dyz_beta2" << endl;
  avx_apply_stencil_z(dy_beta2, dyz_beta2, sz, h);
