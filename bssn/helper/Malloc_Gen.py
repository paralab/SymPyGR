def bssnrhs_cudo_malloc_gen():
    # This function will generate cuda memory allocations for bssnrhs and its deallocation script
    "double * name;"
    "cudaMalloc((void **) &dev_dy, sizeof(double));"
    "cudaFree(&dev_sz);"

    """ -------------------sample allocation--------------------
    int * dev_u_offset;
    cudaStatus = cudaMalloc((void **) &dev_u_offset, sizeof(int));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "dy cudaMalloc failed!\n"); return;}
    """

    output_file1 = open("t_bssnrhs_cuda_malloc.h", 'w')
    output_file2 = open("t_bssnrhs_cuda_mdealloc.h", 'w')
    
    with open("test_bssnrhs_memalloc.h", "r") as f:
        output_file1.write("int size = n * sizeof(double);\n")
        for line in f:
            line = line.strip().split()
            if len(line)>1:
                output_memalloc = "double * %s; cudaStatus = cudaMalloc((void **) &%s, size);\n"%(line[1][1:], line[1][1:])
                output_memalloc_error_check = 'if (cudaStatus != cudaSuccess) {fprintf(stderr, "%s cudaMalloc failed!\\n"); return;}\n\n'%(line[1][1:])
                output_memdealloc = "cudaFree(&%s);\n"%(line[1][1:])

                output_file1.write(output_memalloc)
                output_file1.write(output_memalloc_error_check)
                output_file2.write(output_memdealloc)
    output_file1.close()
    output_file2.close()

def bssnrhs_derivs_gen():
    # This function will convert existing deriv calls into cuda supported version of method calls
    "deriv_y(grad_1_alpha, alpha, hy, sz, bflag);" #original
    "deriv_y(grad_1_B0, dev_var_in, B0Int, hy, sz, bflag);"
    "deriv_y(grad_1_alpha, dev_var_in, dev_alphaInt, dev_dy_hy, dev_sz, bflag, sz);"

    readyToUse = ["alphaInt", "chiInt", "KInt", "gt0Int", "gt1Int", "gt2Int", "gt3Int", "gt4Int", "gt5Int", "beta0Int", "beta1Int", "beta2Int",
    "At0Int", "At1Int", "At2Int", "At3Int", "At4Int", "At5Int", "Gt0Int", "Gt1Int", "Gt2Int", "B0Int", "B1Int", "B2Int"]

    output_file1 = open("t_bssnrhs_cuda_derivs.h", 'w')

    with open("test_bssnrhs_derivs.h", "r") as f:
        for line in f:
            line1 = line.strip().split("(")
            firstPara = line1[1].split(",")[0].strip() 
            line2 = line.strip().split(",")
            if (line2[1].strip()+"Int") in readyToUse:
                output_method_call = "%s(%s, dev_var_in, dev_%sInt, dev_dy_%s, dev_sz, bflag, sz);\n"%(line1[0], firstPara, line2[1].strip(), line2[2].strip())
            else:
                output_method_call = "%s(%s, %s, dev_zero, dev_dy_%s, dev_sz, bflag, sz);\n"%(line1[0], firstPara, line2[1].strip(), line2[2].strip())
            output_file1.write(output_method_call)
    output_file1.close()

def allocate_memory_for_offset_ints():
    """ -------------------sample allocation--------------------
    int * dev_u_offset;
    cudaStatus = cudaMalloc((void **) &dev_u_offset, sizeof(int));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "dy cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_u_offset, &u_offset, sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "dy cudaMemcpy failed!\n"); return;}
    """

    "cudaFree(&grad_2_B1);"

    readyToUse = ["alphaInt", "chiInt", "KInt", "gt0Int", "gt1Int", "gt2Int", "gt3Int", "gt4Int", "gt5Int", "beta0Int", "beta1Int", "beta2Int",
    "At0Int", "At1Int", "At2Int", "At3Int", "At4Int", "At5Int", "Gt0Int", "Gt1Int", "Gt2Int", "B0Int", "B1Int", "B2Int"]

    f = open("offset_ints", "w")
    f2 = open("offset_ints_dealloc", "w")

    for i in readyToUse:
        variableName = "dev_"+i
        line1 = "int * %s;\n"%(variableName)
        line2 = "cudaStatus = cudaMalloc((void **) &%s, sizeof(int));\n"%(variableName)
        line3 = 'if (cudaStatus != cudaSuccess) {fprintf(stderr, "%s cudaMalloc failed!\\n"); return;}\n'%(i)
        line4 = "cudaStatus = cudaMemcpy(%s, &%s, sizeof(int), cudaMemcpyHostToDevice);\n"%(variableName, i)
        line5 = 'if (cudaStatus != cudaSuccess) {fprintf(stderr, "%s cudaMemcpy failed!\\n"); return;}\n\n'%(i)
        f.write(line1)
        f.write(line2)
        f.write(line3)
        f.write(line4)
        f.write(line5)

        line6 = "cudaFree(&%s);\n"%(variableName)
        f2.write(line6)

    f.close()
    f2.close()


def main():
    allocate_memory_for_offset_ints()
    bssnrhs_cudo_malloc_gen()
    bssnrhs_derivs_gen()

main()