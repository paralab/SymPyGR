  std::cout << "calling avx_apply_adv_stencil_x gt0"; 
  avx_apply_adv_stencil_x(gt0_3D, dx_gt0, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x gt1"; 
  avx_apply_adv_stencil_x(gt1_3D, dx_gt1, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x gt2"; 
  avx_apply_adv_stencil_x(gt2_3D, dx_gt2, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x gt3"; 
  avx_apply_adv_stencil_x(gt3_3D, dx_gt3, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x gt4"; 
  avx_apply_adv_stencil_x(gt4_3D, dx_gt4, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x gt5"; 
  avx_apply_adv_stencil_x(gt5_3D, dx_gt5, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x At0"; 
  avx_apply_adv_stencil_x(At0_3D, dx_At0, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x At1"; 
  avx_apply_adv_stencil_x(At1_3D, dx_At1, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x At2"; 
  avx_apply_adv_stencil_x(At2_3D, dx_At2, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x At3"; 
  avx_apply_adv_stencil_x(At3_3D, dx_At3, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x At4"; 
  avx_apply_adv_stencil_x(At4_3D, dx_At4, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x At5"; 
  avx_apply_adv_stencil_x(At5_3D, dx_At5, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_At5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x alpha"; 
  avx_apply_adv_stencil_x(alpha_3D, dx_alpha, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x beta0"; 
  avx_apply_adv_stencil_x(beta0_3D, dx_beta0, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x beta1"; 
  avx_apply_adv_stencil_x(beta1_3D, dx_beta1, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x beta2"; 
  avx_apply_adv_stencil_x(beta2_3D, dx_beta2, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x chi"; 
  avx_apply_adv_stencil_x(chi_3D, dx_chi, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x Gt0"; 
  avx_apply_adv_stencil_x(Gt0_3D, dx_Gt0, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_Gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x Gt1"; 
  avx_apply_adv_stencil_x(Gt1_3D, dx_Gt1, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_Gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x Gt2"; 
  avx_apply_adv_stencil_x(Gt2_3D, dx_Gt2, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_Gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x K"; 
  avx_apply_adv_stencil_x(K_3D, dx_K, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_K[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x B0"; 
  avx_apply_adv_stencil_x(B0_3D, dx_B0, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_B0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x B1"; 
  avx_apply_adv_stencil_x(B1_3D, dx_B1, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_B1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_x B2"; 
  avx_apply_adv_stencil_x(B2_3D, dx_B2, beta0_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dx_B2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y gt0"; 
  avx_apply_adv_stencil_y(gt0_3D, dy_gt0, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y gt1"; 
  avx_apply_adv_stencil_y(gt1_3D, dy_gt1, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y gt2"; 
  avx_apply_adv_stencil_y(gt2_3D, dy_gt2, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y gt3"; 
  avx_apply_adv_stencil_y(gt3_3D, dy_gt3, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y gt4"; 
  avx_apply_adv_stencil_y(gt4_3D, dy_gt4, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y gt5"; 
  avx_apply_adv_stencil_y(gt5_3D, dy_gt5, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y At0"; 
  avx_apply_adv_stencil_y(At0_3D, dy_At0, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y At1"; 
  avx_apply_adv_stencil_y(At1_3D, dy_At1, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y At2"; 
  avx_apply_adv_stencil_y(At2_3D, dy_At2, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y At3"; 
  avx_apply_adv_stencil_y(At3_3D, dy_At3, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y At4"; 
  avx_apply_adv_stencil_y(At4_3D, dy_At4, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y At5"; 
  avx_apply_adv_stencil_y(At5_3D, dy_At5, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_At5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y alpha"; 
  avx_apply_adv_stencil_y(alpha_3D, dy_alpha, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y beta0"; 
  avx_apply_adv_stencil_y(beta0_3D, dy_beta0, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y beta1"; 
  avx_apply_adv_stencil_y(beta1_3D, dy_beta1, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y beta2"; 
  avx_apply_adv_stencil_y(beta2_3D, dy_beta2, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y chi"; 
  avx_apply_adv_stencil_y(chi_3D, dy_chi, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y Gt0"; 
  avx_apply_adv_stencil_y(Gt0_3D, dy_Gt0, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_Gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y Gt1"; 
  avx_apply_adv_stencil_y(Gt1_3D, dy_Gt1, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_Gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y Gt2"; 
  avx_apply_adv_stencil_y(Gt2_3D, dy_Gt2, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_Gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y K"; 
  avx_apply_adv_stencil_y(K_3D, dy_K, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_K[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y B0"; 
  avx_apply_adv_stencil_y(B0_3D, dy_B0, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_B0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y B1"; 
  avx_apply_adv_stencil_y(B1_3D, dy_B1, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_B1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_y B2"; 
  avx_apply_adv_stencil_y(B2_3D, dy_B2, beta1_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dy_B2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z gt0"; 
  avx_apply_adv_stencil_z(gt0_3D, dz_gt0, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z gt1"; 
  avx_apply_adv_stencil_z(gt1_3D, dz_gt1, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z gt2"; 
  avx_apply_adv_stencil_z(gt2_3D, dz_gt2, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z gt3"; 
  avx_apply_adv_stencil_z(gt3_3D, dz_gt3, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z gt4"; 
  avx_apply_adv_stencil_z(gt4_3D, dz_gt4, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z gt5"; 
  avx_apply_adv_stencil_z(gt5_3D, dz_gt5, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_gt5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z At0"; 
  avx_apply_adv_stencil_z(At0_3D, dz_At0, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z At1"; 
  avx_apply_adv_stencil_z(At1_3D, dz_At1, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z At2"; 
  avx_apply_adv_stencil_z(At2_3D, dz_At2, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z At3"; 
  avx_apply_adv_stencil_z(At3_3D, dz_At3, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At3[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z At4"; 
  avx_apply_adv_stencil_z(At4_3D, dz_At4, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At4[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z At5"; 
  avx_apply_adv_stencil_z(At5_3D, dz_At5, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_At5[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z alpha"; 
  avx_apply_adv_stencil_z(alpha_3D, dz_alpha, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_alpha[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z beta0"; 
  avx_apply_adv_stencil_z(beta0_3D, dz_beta0, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_beta0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z beta1"; 
  avx_apply_adv_stencil_z(beta1_3D, dz_beta1, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_beta1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z beta2"; 
  avx_apply_adv_stencil_z(beta2_3D, dz_beta2, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_beta2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z chi"; 
  avx_apply_adv_stencil_z(chi_3D, dz_chi, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_chi[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z Gt0"; 
  avx_apply_adv_stencil_z(Gt0_3D, dz_Gt0, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_Gt0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z Gt1"; 
  avx_apply_adv_stencil_z(Gt1_3D, dz_Gt1, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_Gt1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z Gt2"; 
  avx_apply_adv_stencil_z(Gt2_3D, dz_Gt2, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_Gt2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z K"; 
  avx_apply_adv_stencil_z(K_3D, dz_K, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_K[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z B0"; 
  avx_apply_adv_stencil_z(B0_3D, dz_B0, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_B0[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z B1"; 
  avx_apply_adv_stencil_z(B1_3D, dz_B1, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_B1[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
  std::cout << "calling avx_apply_adv_stencil_z B2"; 
  avx_apply_adv_stencil_z(B2_3D, dz_B2, beta2_3D, sz, h);
{ bool isn=false;
 for (int iii=0; iii<n; ++iii) 
	isn |= std::isnan(dz_B2[iii]);
 std::cout << "  isNAN: " << isn << std::endl;}
