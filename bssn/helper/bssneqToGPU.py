replaceKeyWords = {
    'alpha[pp]' : 'dev_var_in[*dev_alphaInt+pp]',
    'chi[pp]' : 'dev_var_in[*dev_chiInt+pp]',
    'K[pp]' : 'dev_var_in[*dev_KInt+pp]',
    'gt0[pp]' : 'dev_var_in[*dev_gt0Int+pp]',
    'gt1[pp]' : 'dev_var_in[*dev_gt1Int+pp]',
    'gt2[pp]' : 'dev_var_in[*dev_gt2Int+pp]',
    'gt3[pp]' : 'dev_var_in[*dev_gt3Int+pp]',
    'gt4[pp]' : 'dev_var_in[*dev_gt4Int+pp]',
    'gt5[pp]' : 'dev_var_in[*dev_gt5Int+pp]',
    'beta0[pp]' : 'dev_var_in[*dev_beta0Int+pp]',
    'beta1[pp]' : 'dev_var_in[*dev_beta1Int+pp]',
    'beta2[pp]' : 'dev_var_in[*dev_beta2Int+pp]',
    'At0[pp]' : 'dev_var_in[*dev_At0Int+pp]',
    'At1[pp]' : 'dev_var_in[*dev_At1Int+pp]',
    'At2[pp]' : 'dev_var_in[*dev_At2Int+pp]',
    'At3[pp]' : 'dev_var_in[*dev_At3Int+pp]',
    'At4[pp]' : 'dev_var_in[*dev_At4Int+pp]',
    'At5[pp]' : 'dev_var_in[*dev_At5Int+pp]',
    'Gt0[pp]' : 'dev_var_in[*dev_Gt0Int+pp]',
    'Gt1[pp]' : 'dev_var_in[*dev_Gt1Int+pp]',
    'Gt2[pp]' : 'dev_var_in[*dev_Gt2Int+pp]',
    'B0[pp]' : 'dev_var_in[*dev_B0Int+pp]',
    'B1[pp]' : 'dev_var_in[*dev_B1Int+pp]',
    'B2[pp]' : 'dev_var_in[*dev_B2Int+pp]',

    'a_rhs[pp]' : 'dev_var_out[*dev_alphaInt+pp]',
    'chi_rhs[pp]' : 'dev_var_out[*dev_chiInt+pp]',
    'K_rhs[pp]' : 'dev_var_out[*dev_KInt+pp]',
    'gt_rhs00[pp]' : 'dev_var_out[*dev_gt0Int+pp]',
    'gt_rhs01[pp]' : 'dev_var_out[*dev_gt1Int+pp]',
    'gt_rhs02[pp]' : 'dev_var_out[*dev_gt2Int+pp]',
    'gt_rhs11[pp]' : 'dev_var_out[*dev_gt3Int+pp]',
    'gt_rhs12[pp]' : 'dev_var_out[*dev_gt4Int+pp]',
    'gt_rhs22[pp]' : 'dev_var_out[*dev_gt5Int+pp]',
    'b_rhs0[pp]' : 'dev_var_out[*dev_beta0Int+pp]',
    'b_rhs1[pp]' : 'dev_var_out[*dev_beta1Int+pp]',
    'b_rhs2[pp]' : 'dev_var_out[*dev_beta2Int+pp]',
    'At_rhs00[pp]' : 'dev_var_out[*dev_At0Int+pp]',
    'At_rhs01[pp]' : 'dev_var_out[*dev_At1Int+pp]',
    'At_rhs02[pp]' : 'dev_var_out[*dev_At2Int+pp]',
    'At_rhs11[pp]' : 'dev_var_out[*dev_At3Int+pp]',
    'At_rhs12[pp]' : 'dev_var_out[*dev_At4Int+pp]',
    'At_rhs22[pp]' : 'dev_var_out[*dev_At5Int+pp]',
    'Gt_rhs0[pp]' : 'dev_var_out[*dev_Gt0Int+pp]',
    'Gt_rhs1[pp]' : 'dev_var_out[*dev_Gt1Int+pp]',
    'Gt_rhs2[pp]' : 'dev_var_out[*dev_Gt2Int+pp]',
    'B_rhs0[pp]' : 'dev_var_out[*dev_B0Int+pp]',
    'B_rhs1[pp]' : 'dev_var_out[*dev_B1Int+pp]',
    'B_rhs2[pp]' : 'dev_var_out[*dev_B2Int+pp]'
}

output_file = open("t_cuda_bssneqs.h", 'w')

with open('test_bssneqs.cpp', 'r') as file:
    for line in file:
        line = line.strip()
        newLine = ""
        for key in replaceKeyWords.keys():
            if line.count(key)>0:
                lineList = line.split(key)
                if len(lineList)==1:
                    newLine = line
                print(lineList)
            else:
                newLine = line
        break



        output_file.write(newLine)
        output_file.write("\n")



output_file.close()

"dev_var_out[*dev_alphaInt+pp]"
