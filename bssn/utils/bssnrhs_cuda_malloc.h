cudaStatus = cudaMalloc((void **) &grad_0_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_B0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_B0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_B0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_B0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_B0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_B0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_B1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_B1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_B1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_B1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_B1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_B1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_B2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_B2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_B2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_B2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_B2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_B2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_Gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_Gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_Gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_Gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_Gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_Gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_Gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_Gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_Gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_Gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_Gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_Gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_Gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_Gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_Gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_Gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_Gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_Gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_K, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_K cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_K, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_K cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_K, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_K cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_At0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_At0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_At0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_At0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_At0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_At0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_At1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_At1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_At1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_At1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_At1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_At1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_At2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_At2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_At2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_At2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_At2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_At2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_At3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_At3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_At3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_At3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_At3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_At3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_At4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_At4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_At4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_At4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_At4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_At4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_0_At5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_0_At5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_1_At5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_1_At5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad_2_At5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad_2_At5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_gt0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_gt0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_gt1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_gt1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_gt2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_gt2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_gt3, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_gt3 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_gt4, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_gt4 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_gt5, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_gt5 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_chi, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_chi cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_alpha, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_alpha cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_beta0, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_beta0 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_beta1, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_beta1 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_0_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_0_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_1_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_1_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_0_2_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_0_2_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_1_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_1_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_1_2_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_1_2_beta2 cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &grad2_2_2_beta2, size);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "grad2_2_2_beta2 cudaMalloc failed!\n"); }

