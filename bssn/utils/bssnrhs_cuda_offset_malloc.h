int * dev_alphaInt;
cudaStatus = cudaMalloc((void **) &dev_alphaInt, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "alphaInt cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_alphaInt, &alphaInt, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "alphaInt cudaMemcpy failed!\n"); return;}

int * dev_chiInt;
cudaStatus = cudaMalloc((void **) &dev_chiInt, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "chiInt cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_chiInt, &chiInt, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "chiInt cudaMemcpy failed!\n"); return;}

int * dev_KInt;
cudaStatus = cudaMalloc((void **) &dev_KInt, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "KInt cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_KInt, &KInt, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "KInt cudaMemcpy failed!\n"); return;}

int * dev_gt0Int;
cudaStatus = cudaMalloc((void **) &dev_gt0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt0Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_gt0Int, &gt0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt0Int cudaMemcpy failed!\n"); return;}

int * dev_gt1Int;
cudaStatus = cudaMalloc((void **) &dev_gt1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt1Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_gt1Int, &gt1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt1Int cudaMemcpy failed!\n"); return;}

int * dev_gt2Int;
cudaStatus = cudaMalloc((void **) &dev_gt2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt2Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_gt2Int, &gt2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt2Int cudaMemcpy failed!\n"); return;}

int * dev_gt3Int;
cudaStatus = cudaMalloc((void **) &dev_gt3Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt3Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_gt3Int, &gt3Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt3Int cudaMemcpy failed!\n"); return;}

int * dev_gt4Int;
cudaStatus = cudaMalloc((void **) &dev_gt4Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt4Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_gt4Int, &gt4Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt4Int cudaMemcpy failed!\n"); return;}

int * dev_gt5Int;
cudaStatus = cudaMalloc((void **) &dev_gt5Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt5Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_gt5Int, &gt5Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt5Int cudaMemcpy failed!\n"); return;}

int * dev_beta0Int;
cudaStatus = cudaMalloc((void **) &dev_beta0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta0Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_beta0Int, &beta0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta0Int cudaMemcpy failed!\n"); return;}

int * dev_beta1Int;
cudaStatus = cudaMalloc((void **) &dev_beta1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta1Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_beta1Int, &beta1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta1Int cudaMemcpy failed!\n"); return;}

int * dev_beta2Int;
cudaStatus = cudaMalloc((void **) &dev_beta2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta2Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_beta2Int, &beta2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta2Int cudaMemcpy failed!\n"); return;}

int * dev_At0Int;
cudaStatus = cudaMalloc((void **) &dev_At0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At0Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_At0Int, &At0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At0Int cudaMemcpy failed!\n"); return;}

int * dev_At1Int;
cudaStatus = cudaMalloc((void **) &dev_At1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At1Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_At1Int, &At1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At1Int cudaMemcpy failed!\n"); return;}

int * dev_At2Int;
cudaStatus = cudaMalloc((void **) &dev_At2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At2Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_At2Int, &At2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At2Int cudaMemcpy failed!\n"); return;}

int * dev_At3Int;
cudaStatus = cudaMalloc((void **) &dev_At3Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At3Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_At3Int, &At3Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At3Int cudaMemcpy failed!\n"); return;}

int * dev_At4Int;
cudaStatus = cudaMalloc((void **) &dev_At4Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At4Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_At4Int, &At4Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At4Int cudaMemcpy failed!\n"); return;}

int * dev_At5Int;
cudaStatus = cudaMalloc((void **) &dev_At5Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At5Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_At5Int, &At5Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At5Int cudaMemcpy failed!\n"); return;}

int * dev_Gt0Int;
cudaStatus = cudaMalloc((void **) &dev_Gt0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt0Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_Gt0Int, &Gt0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt0Int cudaMemcpy failed!\n"); return;}

int * dev_Gt1Int;
cudaStatus = cudaMalloc((void **) &dev_Gt1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt1Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_Gt1Int, &Gt1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt1Int cudaMemcpy failed!\n"); return;}

int * dev_Gt2Int;
cudaStatus = cudaMalloc((void **) &dev_Gt2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt2Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_Gt2Int, &Gt2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt2Int cudaMemcpy failed!\n"); return;}

int * dev_B0Int;
cudaStatus = cudaMalloc((void **) &dev_B0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B0Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_B0Int, &B0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B0Int cudaMemcpy failed!\n"); return;}

int * dev_B1Int;
cudaStatus = cudaMalloc((void **) &dev_B1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B1Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_B1Int, &B1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B1Int cudaMemcpy failed!\n"); return;}

int * dev_B2Int;
cudaStatus = cudaMalloc((void **) &dev_B2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B2Int cudaMalloc failed!\n"); return;}
cudaStatus = cudaMemcpyAsync(dev_B2Int, &B2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B2Int cudaMemcpy failed!\n"); return;}

int * dev_bflag;
cudaStatus = cudaMalloc((void **) &dev_bflag, sizeof(int));
if (cudaStatus != cudaSuccess)
    fprintf(stderr, "bflag cudaMalloc failed!\n");
cudaStatus = cudaMemcpyAsync(dev_bflag, &bflag, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) 
    fprintf(stderr, "bflag cudaMemcpy failed!\n"); 

