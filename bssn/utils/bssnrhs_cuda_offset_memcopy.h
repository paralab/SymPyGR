cudaStatus = cudaMemcpyAsync(dev_alphaInt, &alphaInt, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "alphaInt cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_chiInt, &chiInt, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "chiInt cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_KInt, &KInt, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "KInt cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_gt0Int, &gt0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt0Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_gt1Int, &gt1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt1Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_gt2Int, &gt2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt2Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_gt3Int, &gt3Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt3Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_gt4Int, &gt4Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt4Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_gt5Int, &gt5Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt5Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_beta0Int, &beta0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta0Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_beta1Int, &beta1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta1Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_beta2Int, &beta2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta2Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_At0Int, &At0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At0Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_At1Int, &At1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At1Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_At2Int, &At2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At2Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_At3Int, &At3Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At3Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_At4Int, &At4Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At4Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_At5Int, &At5Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At5Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_Gt0Int, &Gt0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt0Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_Gt1Int, &Gt1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt1Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_Gt2Int, &Gt2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt2Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_B0Int, &B0Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B0Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_B1Int, &B1Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B1Int cudaMemcpyAsync failed!\n"); return;}

cudaStatus = cudaMemcpyAsync(dev_B2Int, &B2Int, sizeof(int), cudaMemcpyHostToDevice, stream);
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B2Int cudaMemcpyAsync failed!\n"); return;}

