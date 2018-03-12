def bssnrhs_cudo_malloc_gen():
    # This function will generate cuda memory allocations for bssnrhs and its deallocation script
    "double * name;"
    "cudaMalloc((void **) &dev_dy, sizeof(double));"
    "cudaFree(&dev_sz);"

    output_file1 = open("bssnrhs_cuda_malloc.h", 'w')
    output_file2 = open("bssnrhs_cuda_mdealloc.h", 'w')
    
    with open("test_bssnrhs_memalloc.h", "r") as f:
        output_file1.write("int size = n * sizeof(double);\n")
        for line in f:
            line = line.strip().split()
            if len(line)>1:
                output_memalloc = "double * %s; cudaMalloc((void **) &%s, size);\n"%(line[1][1:], line[1][1:])
                output_memdealloc = "cudaFree(&%s);\n"%(line[1][1:])

                output_file1.write(output_memalloc)
                output_file2.write(output_memdealloc)

def bssnrhs_derivs_gen():
    # This function will convert existing deriv calls into cuda supported version of method calls
    "deriv_y(dev_var_in, alphaInt, hy, sz, bflag);"

    readyToUse = ["alphaInt", "chiInt", "KInt", "gt0Int", "gt1Int", "gt2Int", "gt3Int", "gt4Int", "gt5Int", "beta0Int", "beta1Int", "beta2Int",
    "At0Int", "At1Int", "At2Int", "At3Int", "At4Int", "At5Int", "Gt0Int", "Gt1Int", "Gt2Int", "B0Int", "B1Int", "B2Int"]

    output_file1 = open("test_bssnrhs_cuda_derivs.h", 'w')

    with open("bssnrhs_derivs.h", "r") as f:
        for line in f:
            line1 = line.strip().split("(")
            line2 = line.strip().split(",")
            if (line2[1].strip()+"Int") in readyToUse:
                output_method_call = "%s(dev_var_in, %sInt, %s, sz, bflag);\n"%(line1[0], line2[1].strip(), line2[2].strip())
            else:
                output_method_call = "%s(%s, 0, %s, sz, bflag);\n"%(line1[0], line2[1].strip(), line2[2].strip())
            output_file1.write(output_method_call)


def main():
    bssnrhs_cudo_malloc_gen()
    bssnrhs_derivs_gen()

main()