#----------------------------------------------------------------------------------------------------------------

import re

line = "DENDRO_28 + DENDRO_30 + DENDRO_32 - DENDRO_34 - DENDRO_36"

dendro_var_list = list(set(re.findall("DENDRO_[0-9]+", line)))
dendro_var_list = sorted(dendro_var_list, key=lambda x: int(x[7:]))

with open("tempory.txt", "w") as data:
    for dend_var in dendro_var_list:
        sample_gpu = 'printf("%.20f -- {}\n", {});'.format(dend_var, dend_var) # printf("GPU - {}=%f\n", {});
        byte_gpu = str(str.encode(sample_gpu))

        data.write("%s\n"%(byte_gpu[2:-1]))




#----------------------------------------------------------------------------------------------------------------

# with open("tempory.txt", "w") as data:
#     for i in range(0, 992):
#         sample_gpu = 'printf("%.20f -- DENDRO_{}\n", DENDRO_{});'.format(i, i) # printf("GPU - {}=%f\n", {});
#         byte_gpu = str(str.encode(sample_gpu))
#         data.write("%s\n"%(byte_gpu[2:-1]))

#----------------------------------------------------------------------------------------------------------------
