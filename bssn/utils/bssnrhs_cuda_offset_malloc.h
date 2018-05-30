cudaStatus = cudaMalloc((void **) &dev_alphaInt, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "alphaInt cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_chiInt, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "chiInt cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_KInt, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "KInt cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_gt0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt0Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_gt1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt1Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_gt2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt2Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_gt3Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt3Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_gt4Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt4Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_gt5Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "gt5Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_beta0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta0Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_beta1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta1Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_beta2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "beta2Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_At0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At0Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_At1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At1Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_At2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At2Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_At3Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At3Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_At4Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At4Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_At5Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "At5Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_Gt0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt0Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_Gt1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt1Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_Gt2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "Gt2Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_B0Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B0Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_B1Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B1Int cudaMalloc failed!\n"); }

cudaStatus = cudaMalloc((void **) &dev_B2Int, sizeof(int));
if (cudaStatus != cudaSuccess) {fprintf(stderr, "B2Int cudaMalloc failed!\n"); }

