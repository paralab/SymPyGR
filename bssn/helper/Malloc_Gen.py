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
        # output_file1.write("int size = n * sizeof(double);\n")
        for line in f:
            line = line.strip().split()
            if len(line)>1:
                # output_memalloc = "double * %s;\n"%(line[1][1:])
                # output_memalloc = "double * %s; cudaStatus = cudaMalloc((void **) &%s, size);\n"%(line[1][1:], line[1][1:])
                # output_memalloc = 'CHECK_ERROR(cudaMalloc((void **) &%s, size), "%s");\n'%(line[1][1:], line[1][1:])
                # output_memalloc = "cudaStatus = cudaMalloc((void **) &%s, size);\n"%(line[1][1:])
                # output_memalloc_error_check = 'if (cudaStatus != cudaSuccess) {fprintf(stderr, "%s cudaMalloc failed!\\n"); return;}\n\n'%(line[1][1:])
                output_memdealloc = "cudaStatus = cudaFree(%s); if (cudaStatus != cudaSuccess) {fprintf(stderr, \"%s cudafree failed!\\n\");}\n"%(line[1][1:], line[1][1:])
                # output_memdealloc = 'CHECK_ERROR(cudaFree(%s), "%s cudafree");\n'%(line[1][1:], line[1][1:])

                # output_file1.write(output_memalloc)
                # output_file1.write(output_memalloc_error_check)
                # output_file2.write(output_memdealloc)
    output_file1.close()
    output_file2.close()

def bssnrhs_derivs_gen():
    # This function will convert existing deriv calls into cuda supported version of method calls
    "deriv_y(grad_1_alpha, alpha, hy, sz, bflag);" #original
    "deriv_y(grad_1_B0, dev_var_in, B0Int, hy, sz, bflag);"
    "deriv_y(grad_1_alpha, dev_var_in, dev_alphaInt, dev_dy_hy, dev_sz, bflag, sz);"
    # deriv_x(grad_0_alpha, dev_var_in, alphaInt, hx, bflag, sz, stream);

    readyToUse = ["alphaInt", "chiInt", "KInt", "gt0Int", "gt1Int", "gt2Int", "gt3Int", "gt4Int", "gt5Int", "beta0Int", "beta1Int", "beta2Int",
    "At0Int", "At1Int", "At2Int", "At3Int", "At4Int", "At5Int", "Gt0Int", "Gt1Int", "Gt2Int", "B0Int", "B1Int", "B2Int"]

    output_file1 = open("t_bssnrhs_cuda_derivs.h", 'w')

    with open("test_bssnrhs_derivs.h", "r") as f:
        for line in f:
            line1 = line.strip().split("(")
            firstPara = line1[1].split(",")[0].strip() 
            line2 = line.strip().split(",")
            if (line2[1].strip()+"Int") in readyToUse:
                output_method_call = "%s(%s, dev_var_in, %sInt, %s, bflag, sz, stream);\n"%(line1[0], firstPara, line2[1].strip(), line2[2].strip())
            else:
                output_method_call = "%s(%s, %s, 0, %s, bflag, sz, stream);\n"%(line1[0], firstPara, line2[1].strip(), line2[2].strip())
            output_file1.write(output_method_call)
    output_file1.close()

def allocate_memory_for_offset_ints():
    """ -------------------sample allocation--------------------
    int * dev_u_offset;
    cudaStatus = cudaMalloc((void **) &dev_u_offset, sizeof(int));
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "dy cudaMalloc failed!\n"); return;}
    cudaStatus = cudaMemcpy(dev_u_offset, &u_offset, sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "dy cudaMemcpy failed!\n"); return;}

    cudaStatus = cudaMemcpyAsync(dev_KInt, &KInt, sizeof(int), cudaMemcpyHostToDevice, stream);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "KInt cudaMemcpy failed!\n"); return;}
    """

    "cudaFree(grad_2_B1);"

    readyToUse = ["alphaInt", "chiInt", "KInt", "gt0Int", "gt1Int", "gt2Int", "gt3Int", "gt4Int", "gt5Int", "beta0Int", "beta1Int", "beta2Int",
    "At0Int", "At1Int", "At2Int", "At3Int", "At4Int", "At5Int", "Gt0Int", "Gt1Int", "Gt2Int", "B0Int", "B1Int", "B2Int"]

    f = open("offset_ints", "w")
    f2 = open("offset_ints_dealloc", "w")

    for i in readyToUse:
        variableName = "dev_"+i
        # line1 = "int * %s;\n"%(variableName)
        # line2 = 'CHECK_ERROR(cudaMalloc((void **) &%s, sizeof(int)), "%s");\n'%(variableName, variableName)
        # line2 = "cudaStatus = cudaMalloc((void **) &%s, sizeof(int));\n"%(variableName)
        # line3 = 'if (cudaStatus != cudaSuccess) {fprintf(stderr, "%s cudaMalloc failed!\\n"); return;}\n\n'%(i)
        # line4 = "cudaStatus = cudaMemcpyAsync(%s, &%s, sizeof(int), cudaMemcpyHostToDevice, stream);\n"%(variableName, i)
        line4 = 'CHECK_ERROR(cudaMemcpyAsync(%s, &%s, sizeof(int), cudaMemcpyHostToDevice, stream), "%s cudaMemcpyHostToDevice");\n'%(variableName, i, variableName)
        # line5 = 'if (cudaStatus != cudaSuccess) {fprintf(stderr, "%s cudaMemcpyAsync failed!\\n"); return;}\n\n'%(i)
        # f.write(line1)
        # f.write(line2)
        # f.write(line3)
        f.write(line4)
        # f.write(line5)

        # line6 = "cudaStatus = cudaFree(%s);\nif (cudaStatus != cudaSuccess) {fprintf(stderr, \"%s cudafree failed!\\n\");}\n\n"%(variableName, variableName)
        # line6 = 'CHECK_ERROR(cudaFree(%s), "%s cudafree");\n'%(variableName, variableName)
        # f2.write(line6)

    f.close()
    f2.close()

def bssnrhs_adv_derivs_gen():
    """
    adv_deriv_x(agrad_0_gt0, gt0, hx, sz, beta0, bflag); // old

    adv_deriv_x(agrad_0_gt0, dev_var_in, dev_gt0Int, dev_dy_hx, dev_sz, dev_beta0Int, dev_bflag, sz); // new

    adv_deriv_x(agrad_0_gt0, dev_var_in, gt0Int, hx, beta0Int, bflag, sz, stream); 
    """

    output_file1 = open("t_bssnrhs_derivs_adv.h", 'w')

    with open("test_bssnrhs_derivs_adv.h", "r") as f:
        for line in f:
            line1 = line.strip().split("(")
            method_name = line1[0]

            line1 = line1[1].split(",")
            para1 = line1[0].strip() # agrad_0_gt0
            para2 = line1[1].strip() # gt0
            para3 = line1[2].strip() # hx
            para4 = line1[3].strip() # sz
            para5 = line1[4].strip() # beta0
            para6 = line1[5].split(")")[0] # bflag

            # print(method_name, para1, para2, para3, para4, para5, para6)

            output_method_call = "%s(%s, dev_var_in, %sInt, %s, %sInt, bflag, sz, stream);\n"%(method_name, para1, para2, para3, para5)

            # print(output_method_call)

            output_file1.write(output_method_call)

    output_file1.close()



def bssnrhs_cudo_malloc_adv_gen():
    """
        // old
         double *agrad_0_B0 = (double *) malloc(bytes);
         double *agrad_1_B0 = (double *) malloc(bytes);


        // new
        double *agrad_2_gt2; cudaStatus = cudaMalloc((void**)&agrad_2_gt2, size);
        if (cudaStatus != cudaSuccess) {fprintf(stderr, "agrad_2_gt2 cudaMalloc failed!\n"); return;}
    """
    output_file1 = open("t_bssnrhs_cuda_malloc_adv.h", 'w')
    output_file2 = open("t_bssnrhs_cuda_mdealloc_adv.h", 'w')
    
    with open("test_bssnrhs_memalloc_adv.h", "r") as f:
        # output_file1.write("int size = n * sizeof(double);\n")
        for line in f:
            line = line.strip().split()
            if len(line)>1:
                variable_name = line[1][1:].strip()
                # output_memalloc = "double * %s;\n"%(variable_name)
                # output_memalloc = "double * %s; cudaStatus = cudaMalloc((void **) &%s, size);\n"%(variable_name, variable_name)
                # output_memalloc = 'CHECK_ERROR(cudaMalloc((void **) &%s, size), "%s");\n'%(variable_name, variable_name)
                # output_memalloc = "cudaStatus = cudaMalloc((void **) &%s, size);\n"%(variable_name)
                # output_memalloc_error_check = 'if (cudaStatus != cudaSuccess) {fprintf(stderr, "%s cudaMalloc failed!\\n"); return;}\n\n'%(variable_name)
                # output_memdealloc = "cudaStatus = cudaFree(%s); if (cudaStatus != cudaSuccess) {fprintf(stderr, \"%s cudafree failed!\\n\");}\n"%(variable_name, variable_name)
                output_memdealloc = 'CHECK_ERROR(cudaFree(%s), "%s cudafree");\n'%(variable_name, variable_name)
                # output_file1.write(output_memalloc)
                # output_file1.write(output_memalloc_error_check)
                output_file2.write(output_memdealloc)
    output_file1.close()
    output_file2.close()

def bssnrhs_ko_derivs_gen():
    """
      ko_deriv_x(grad_0_gt0, gt0, hx, sz, bflag); // old

      ko_deriv_x(grad_0_gt0, dev_var_in, dev_gt0Int, dev_dy_hx, dev_sz, dev_bflag, sz); // new

      ko_deriv_x(grad_0_gt0, dev_var_in, gt0Int, hx, bflag, sz, stream);
    """

    output_file1 = open("t_bssnrhs_ko_derivs.h", 'w')

    with open("test_bssnrhs_ko_derivs.h", "r") as f:
        for line in f:
            line1 = line.strip().split("(")
            method_name = line1[0]

            line1 = line1[1].split(",")
            para1 = line1[0].strip() # grad_0_gt0
            para2 = line1[1].strip() # gt0
            para3 = line1[2].strip() # hx
            para4 = line1[3].strip() # sz
            para5 = line1[4].split(")")[0].strip() # bflag

            # print(method_name, para1, para2, para3, para4, para5)

            output_method_call = "%s(%s, dev_var_in, %sInt, %s, %s, sz, stream);\n"%(method_name, para1, para2, para3, para5)

            # print(output_method_call)

            output_file1.write(output_method_call)

    output_file1.close()

def main():
    # allocate_memory_for_offset_ints()
    # bssnrhs_cudo_malloc_gen()
    bssnrhs_derivs_gen()
    # bssnrhs_cudo_malloc_adv_gen()
    bssnrhs_adv_derivs_gen()
    bssnrhs_ko_derivs_gen()

main()
