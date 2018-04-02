
#include "test_param.h"
#include <fstream>
using namespace std;
if (test == 1)
{
  int n = sz[0] * sz[1] * sz[2];
  ofstream myfile("ouput_cuda.txt");
  if ((myfile.is_open()))
  {
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_alpha_cpu + i);
    }
    /*
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_B0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_B0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_B0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_B1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_B1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_B1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_B2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_B2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_B2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_Gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_Gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_Gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_Gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_Gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_Gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_Gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_Gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_Gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_K_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_K_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_K_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_At0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_At0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_At0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_At1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_At1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_At1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_At2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_At2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_At2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_At3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_At3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_At3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_At4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_At4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_At4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_At5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_1_At5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_2_At5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_gt0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_gt1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_gt2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_gt3_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_gt4_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_gt5_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_chi_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_alpha_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_beta0_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_beta1_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_0_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_1_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_0_2_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_1_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_1_2_beta2_cpu + i);
    }
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad2_2_2_beta2_cpu + i);
    }
    */
    myfile.close();
  }
  free(grad_0_alpha_cpu);
//   free(grad_1_alpha_cpu);
//   free(grad_2_alpha_cpu);
//   free(grad_0_beta0_cpu);
//   free(grad_1_beta0_cpu);
//   free(grad_2_beta0_cpu);
//   free(grad_0_beta1_cpu);
//   free(grad_1_beta1_cpu);
//   free(grad_2_beta1_cpu);
//   free(grad_0_beta2_cpu);
//   free(grad_1_beta2_cpu);
//   free(grad_2_beta2_cpu);
//   free(grad_0_B0_cpu);
//   free(grad_1_B0_cpu);
//   free(grad_2_B0_cpu);
//   free(grad_0_B1_cpu);
//   free(grad_1_B1_cpu);
//   free(grad_2_B1_cpu);
//   free(grad_0_B2_cpu);
//   free(grad_1_B2_cpu);
//   free(grad_2_B2_cpu);
//   free(grad_0_chi_cpu);
//   free(grad_1_chi_cpu);
//   free(grad_2_chi_cpu);
//   free(grad_0_Gt0_cpu);
//   free(grad_1_Gt0_cpu);
//   free(grad_2_Gt0_cpu);
//   free(grad_0_Gt1_cpu);
//   free(grad_1_Gt1_cpu);
//   free(grad_2_Gt1_cpu);
//   free(grad_0_Gt2_cpu);
//   free(grad_1_Gt2_cpu);
//   free(grad_2_Gt2_cpu);
//   free(grad_0_K_cpu);
//   free(grad_1_K_cpu);
//   free(grad_2_K_cpu);
//   free(grad_0_gt0_cpu);
//   free(grad_1_gt0_cpu);
//   free(grad_2_gt0_cpu);
//   free(grad_0_gt1_cpu);
//   free(grad_1_gt1_cpu);
//   free(grad_2_gt1_cpu);
//   free(grad_0_gt2_cpu);
//   free(grad_1_gt2_cpu);
//   free(grad_2_gt2_cpu);
//   free(grad_0_gt3_cpu);
//   free(grad_1_gt3_cpu);
//   free(grad_2_gt3_cpu);
//   free(grad_0_gt4_cpu);
//   free(grad_1_gt4_cpu);
//   free(grad_2_gt4_cpu);
//   free(grad_0_gt5_cpu);
//   free(grad_1_gt5_cpu);
//   free(grad_2_gt5_cpu);
//   free(grad_0_At0_cpu);
//   free(grad_1_At0_cpu);
//   free(grad_2_At0_cpu);
//   free(grad_0_At1_cpu);
//   free(grad_1_At1_cpu);
//   free(grad_2_At1_cpu);
//   free(grad_0_At2_cpu);
//   free(grad_1_At2_cpu);
//   free(grad_2_At2_cpu);
//   free(grad_0_At3_cpu);
//   free(grad_1_At3_cpu);
//   free(grad_2_At3_cpu);
//   free(grad_0_At4_cpu);
//   free(grad_1_At4_cpu);
//   free(grad_2_At4_cpu);
//   free(grad_0_At5_cpu);
//   free(grad_1_At5_cpu);
//   free(grad_2_At5_cpu);
//   free(grad2_0_0_gt0_cpu);
//   free(grad2_0_1_gt0_cpu);
//   free(grad2_0_2_gt0_cpu);
//   free(grad2_1_1_gt0_cpu);
//   free(grad2_1_2_gt0_cpu);
//   free(grad2_2_2_gt0_cpu);
//   free(grad2_0_0_gt1_cpu);
//   free(grad2_0_1_gt1_cpu);
//   free(grad2_0_2_gt1_cpu);
//   free(grad2_1_1_gt1_cpu);
//   free(grad2_1_2_gt1_cpu);
//   free(grad2_2_2_gt1_cpu);
//   free(grad2_0_0_gt2_cpu);
//   free(grad2_0_1_gt2_cpu);
//   free(grad2_0_2_gt2_cpu);
//   free(grad2_1_1_gt2_cpu);
//   free(grad2_1_2_gt2_cpu);
//   free(grad2_2_2_gt2_cpu);
//   free(grad2_0_0_gt3_cpu);
//   free(grad2_0_1_gt3_cpu);
//   free(grad2_0_2_gt3_cpu);
//   free(grad2_1_1_gt3_cpu);
//   free(grad2_1_2_gt3_cpu);
//   free(grad2_2_2_gt3_cpu);
//   free(grad2_0_0_gt4_cpu);
//   free(grad2_0_1_gt4_cpu);
//   free(grad2_0_2_gt4_cpu);
//   free(grad2_1_1_gt4_cpu);
//   free(grad2_1_2_gt4_cpu);
//   free(grad2_2_2_gt4_cpu);
//   free(grad2_0_0_gt5_cpu);
//   free(grad2_0_1_gt5_cpu);
//   free(grad2_0_2_gt5_cpu);
//   free(grad2_1_1_gt5_cpu);
//   free(grad2_1_2_gt5_cpu);
//   free(grad2_2_2_gt5_cpu);
//   free(grad2_0_0_chi_cpu);
//   free(grad2_0_1_chi_cpu);
//   free(grad2_0_2_chi_cpu);
//   free(grad2_1_1_chi_cpu);
//   free(grad2_1_2_chi_cpu);
//   free(grad2_2_2_chi_cpu);
//   free(grad2_0_0_alpha_cpu);
//   free(grad2_0_1_alpha_cpu);
//   free(grad2_0_2_alpha_cpu);
//   free(grad2_1_1_alpha_cpu);
//   free(grad2_1_2_alpha_cpu);
//   free(grad2_2_2_alpha_cpu);
//   free(grad2_0_0_beta0_cpu);
//   free(grad2_0_1_beta0_cpu);
//   free(grad2_0_2_beta0_cpu);
//   free(grad2_1_1_beta0_cpu);
//   free(grad2_1_2_beta0_cpu);
//   free(grad2_2_2_beta0_cpu);
//   free(grad2_0_0_beta1_cpu);
//   free(grad2_0_1_beta1_cpu);
//   free(grad2_0_2_beta1_cpu);
//   free(grad2_1_1_beta1_cpu);
//   free(grad2_1_2_beta1_cpu);
//   free(grad2_2_2_beta1_cpu);
//   free(grad2_0_0_beta2_cpu);
//   free(grad2_0_1_beta2_cpu);
//   free(grad2_0_2_beta2_cpu);
//   free(grad2_1_1_beta2_cpu);
//   free(grad2_1_2_beta2_cpu);
//   free(grad2_2_2_beta2_cpu);
 }
