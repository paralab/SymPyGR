#----------------------------------------------------------------------------------------------------------------

# import re

# line = "DENDRO_28 + DENDRO_30 + DENDRO_32 - DENDRO_34 - DENDRO_36"

# dendro_var_list = list(set(re.findall("DENDRO_[0-9]+", line)))
# dendro_var_list = sorted(dendro_var_list, key=lambda x: int(x[7:]))

# with open("tempory.txt", "w") as data:
#     for dend_var in dendro_var_list:
#         sample_gpu = 'printf("%.20f -- {}\n", {});'.format(dend_var, dend_var) # printf("GPU - {}=%f\n", {});
#         byte_gpu = str(str.encode(sample_gpu))

#         data.write("%s\n"%(byte_gpu[2:-1]))




#----------------------------------------------------------------------------------------------------------------

# with open("tempory.txt", "w") as data:
#     for i in range(0, 992):
#         sample_gpu = 'printf("%.20f -- DENDRO_{}\n", DENDRO_{});'.format(i, i) # printf("GPU - {}=%f\n", {});
#         byte_gpu = str(str.encode(sample_gpu))
#         data.write("%s\n"%(byte_gpu[2:-1]))

#----------------------------------------------------------------------------------------------------------------

def bssnrhs_cudo_malloc_gen():
    # This function will generate cuda memory allocations for bssnrhs and its deallocation script
    "double *grad_1_alpha = (double *) malloc(bytes);"
    """test_file_write::appendToFile("output_cpu.txt", grad_1_K, n);"""

    """
    cudaStatus = cudaMemcpy(host_array_cpu, grad_1_K, size, cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {fprintf(stderr, "TEST: host_array_cpu cudaMemcpy from GPU to CPU failed!\n"); return;}
    test_file_write::writeToFile("output_cuda.txt", host_array_cpu, n);
    """

    output_file1 = open("tempory.txt", 'w')
    with open("test_bssnrhs_memalloc_adv.h", "r") as f:
        for line in f:
            line = line.strip().split()
            if len(line)>1:
                variable = line[1][1:]
                line_to_write = 'test_file_write::appendToFile("output_cpu.txt", %s, n);\n'%(variable) # CPU for derivs
                # line_to_write = 'cudaStatus = cudaMemcpy(host_array_cpu, %s, size, cudaMemcpyDeviceToHost);\nif (cudaStatus != cudaSuccess) {fprintf(stderr, "TEST: host_array_cpu cudaMemcpy from GPU %s to CPU failed!"); return;}\ntest_file_write::appendToFile("output_cuda.txt", host_array_cpu, n);\n'%(variable, variable)

                # print(line_to_write)

                output_file1.write(line_to_write)

    output_file1.close()

bssnrhs_cudo_malloc_gen()