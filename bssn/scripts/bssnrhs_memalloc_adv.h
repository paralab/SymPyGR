#ifdef BSSN_USE_ADVECTIVE_DERIVS
  double *agrad_0_gt0 = __mem_pool->allocate(n);
  double *agrad_1_gt0 = __mem_pool->allocate(n);
  double *agrad_2_gt0 = __mem_pool->allocate(n);
  double *agrad_0_gt1 = __mem_pool->allocate(n);
  double *agrad_1_gt1 = __mem_pool->allocate(n);
  double *agrad_2_gt1 = __mem_pool->allocate(n);
  double *agrad_0_gt2 = __mem_pool->allocate(n);
  double *agrad_1_gt2 = __mem_pool->allocate(n);
  double *agrad_2_gt2 = __mem_pool->allocate(n);
  double *agrad_0_gt3 = __mem_pool->allocate(n);
  double *agrad_1_gt3 = __mem_pool->allocate(n);
  double *agrad_2_gt3 = __mem_pool->allocate(n);
  double *agrad_0_gt4 = __mem_pool->allocate(n);
  double *agrad_1_gt4 = __mem_pool->allocate(n);
  double *agrad_2_gt4 = __mem_pool->allocate(n);
  double *agrad_0_gt5 = __mem_pool->allocate(n);
  double *agrad_1_gt5 = __mem_pool->allocate(n);
  double *agrad_2_gt5 = __mem_pool->allocate(n);
  double *agrad_0_At0 = __mem_pool->allocate(n);
  double *agrad_1_At0 = __mem_pool->allocate(n);
  double *agrad_2_At0 = __mem_pool->allocate(n);
  double *agrad_0_At1 = __mem_pool->allocate(n);
  double *agrad_1_At1 = __mem_pool->allocate(n);
  double *agrad_2_At1 = __mem_pool->allocate(n);
  double *agrad_0_At2 = __mem_pool->allocate(n);
  double *agrad_1_At2 = __mem_pool->allocate(n);
  double *agrad_2_At2 = __mem_pool->allocate(n);
  double *agrad_0_At3 = __mem_pool->allocate(n);
  double *agrad_1_At3 = __mem_pool->allocate(n);
  double *agrad_2_At3 = __mem_pool->allocate(n);
  double *agrad_0_At4 = __mem_pool->allocate(n);
  double *agrad_1_At4 = __mem_pool->allocate(n);
  double *agrad_2_At4 = __mem_pool->allocate(n);
  double *agrad_0_At5 = __mem_pool->allocate(n);
  double *agrad_1_At5 = __mem_pool->allocate(n);
  double *agrad_2_At5 = __mem_pool->allocate(n);
  double *agrad_0_alpha = __mem_pool->allocate(n);
  double *agrad_1_alpha = __mem_pool->allocate(n);
  double *agrad_2_alpha = __mem_pool->allocate(n);
  double *agrad_0_beta0 = __mem_pool->allocate(n);
  double *agrad_1_beta0 = __mem_pool->allocate(n);
  double *agrad_2_beta0 = __mem_pool->allocate(n);
  double *agrad_0_beta1 = __mem_pool->allocate(n);
  double *agrad_1_beta1 = __mem_pool->allocate(n);
  double *agrad_2_beta1 = __mem_pool->allocate(n);
  double *agrad_0_beta2 = __mem_pool->allocate(n);
  double *agrad_1_beta2 = __mem_pool->allocate(n);
  double *agrad_2_beta2 = __mem_pool->allocate(n);
  double *agrad_0_chi = __mem_pool->allocate(n);
  double *agrad_1_chi = __mem_pool->allocate(n);
  double *agrad_2_chi = __mem_pool->allocate(n);
  double *agrad_0_Gt0 = __mem_pool->allocate(n);
  double *agrad_1_Gt0 = __mem_pool->allocate(n);
  double *agrad_2_Gt0 = __mem_pool->allocate(n);
  double *agrad_0_Gt1 = __mem_pool->allocate(n);
  double *agrad_1_Gt1 = __mem_pool->allocate(n);
  double *agrad_2_Gt1 = __mem_pool->allocate(n);
  double *agrad_0_Gt2 = __mem_pool->allocate(n);
  double *agrad_1_Gt2 = __mem_pool->allocate(n);
  double *agrad_2_Gt2 = __mem_pool->allocate(n);
  double *agrad_0_K = __mem_pool->allocate(n);
  double *agrad_1_K = __mem_pool->allocate(n);
  double *agrad_2_K = __mem_pool->allocate(n);
  double *agrad_0_B0 = __mem_pool->allocate(n);
  double *agrad_1_B0 = __mem_pool->allocate(n);
  double *agrad_2_B0 = __mem_pool->allocate(n);
  double *agrad_0_B1 = __mem_pool->allocate(n);
  double *agrad_1_B1 = __mem_pool->allocate(n);
  double *agrad_2_B1 = __mem_pool->allocate(n);
  double *agrad_0_B2 = __mem_pool->allocate(n);
  double *agrad_1_B2 = __mem_pool->allocate(n);
  double *agrad_2_B2 = __mem_pool->allocate(n);
#else
  double *agrad_0_gt0 = grad_0_gt0;
  double *agrad_1_gt0 = grad_1_gt0;
  double *agrad_2_gt0 = grad_2_gt0;
  double *agrad_0_gt1 = grad_0_gt1;
  double *agrad_1_gt1 = grad_1_gt1;
  double *agrad_2_gt1 = grad_2_gt1;
  double *agrad_0_gt2 = grad_0_gt2;
  double *agrad_1_gt2 = grad_1_gt2;
  double *agrad_2_gt2 = grad_2_gt2;
  double *agrad_0_gt3 = grad_0_gt3;
  double *agrad_1_gt3 = grad_1_gt3;
  double *agrad_2_gt3 = grad_2_gt3;
  double *agrad_0_gt4 = grad_0_gt4;
  double *agrad_1_gt4 = grad_1_gt4;
  double *agrad_2_gt4 = grad_2_gt4;
  double *agrad_0_gt5 = grad_0_gt5;
  double *agrad_1_gt5 = grad_1_gt5;
  double *agrad_2_gt5 = grad_2_gt5;
  double *agrad_0_At0 = grad_0_At0;
  double *agrad_1_At0 = grad_1_At0;
  double *agrad_2_At0 = grad_2_At0;
  double *agrad_0_At1 = grad_0_At1;
  double *agrad_1_At1 = grad_1_At1;
  double *agrad_2_At1 = grad_2_At1;
  double *agrad_0_At2 = grad_0_At2;
  double *agrad_1_At2 = grad_1_At2;
  double *agrad_2_At2 = grad_2_At2;
  double *agrad_0_At3 = grad_0_At3;
  double *agrad_1_At3 = grad_1_At3;
  double *agrad_2_At3 = grad_2_At3;
  double *agrad_0_At4 = grad_0_At4;
  double *agrad_1_At4 = grad_1_At4;
  double *agrad_2_At4 = grad_2_At4;
  double *agrad_0_At5 = grad_0_At5;
  double *agrad_1_At5 = grad_1_At5;
  double *agrad_2_At5 = grad_2_At5;
  double *agrad_0_alpha = grad_0_alpha;
  double *agrad_1_alpha = grad_1_alpha;
  double *agrad_2_alpha = grad_2_alpha;
  double *agrad_0_beta0 = grad_0_beta0;
  double *agrad_1_beta0 = grad_1_beta0;
  double *agrad_2_beta0 = grad_2_beta0;
  double *agrad_0_beta1 = grad_0_beta1;
  double *agrad_1_beta1 = grad_1_beta1;
  double *agrad_2_beta1 = grad_2_beta1;
  double *agrad_0_beta2 = grad_0_beta2;
  double *agrad_1_beta2 = grad_1_beta2;
  double *agrad_2_beta2 = grad_2_beta2;
  double *agrad_0_chi = grad_0_chi;
  double *agrad_1_chi = grad_1_chi;
  double *agrad_2_chi = grad_2_chi;
  double *agrad_0_Gt0 = grad_0_Gt0;
  double *agrad_1_Gt0 = grad_1_Gt0;
  double *agrad_2_Gt0 = grad_2_Gt0;
  double *agrad_0_Gt1 = grad_0_Gt1;
  double *agrad_1_Gt1 = grad_1_Gt1;
  double *agrad_2_Gt1 = grad_2_Gt1;
  double *agrad_0_Gt2 = grad_0_Gt2;
  double *agrad_1_Gt2 = grad_1_Gt2;
  double *agrad_2_Gt2 = grad_2_Gt2;
  double *agrad_0_K = grad_0_K;
  double *agrad_1_K = grad_1_K;
  double *agrad_2_K = grad_2_K;
  double *agrad_0_B0 = grad_0_B0;
  double *agrad_1_B0 = grad_1_B0;
  double *agrad_2_B0 = grad_2_B0;
  double *agrad_0_B1 = grad_0_B1;
  double *agrad_1_B1 = grad_1_B1;
  double *agrad_2_B1 = grad_2_B1;
  double *agrad_0_B2 = grad_0_B2;
  double *agrad_1_B2 = grad_1_B2;
  double *agrad_2_B2 = grad_2_B2;
#endif
