#include <fstream>
using namespace std;
if (test == 1)
{
  ofstream myfile("output_cpu.txt");
  if ((myfile.is_open()))
  {
    for (int i = 0; i < n; i++)
    {
      myfile << *(grad_0_alpha + i);
    }

    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_B0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_B0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_B0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_B1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_B1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_B1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_B2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_B2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_B2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_Gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_Gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_Gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_Gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_Gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_Gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_Gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_Gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_Gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_K + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_K + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_K + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_At0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_At0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_At0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_At1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_At1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_At1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_At2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_At2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_At2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_At3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_At3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_At3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_At4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_At4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_At4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_0_At5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_1_At5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad_2_At5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_gt0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_gt1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_gt2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_gt3 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_gt4 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_gt5 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_chi + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_alpha + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_beta0 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_beta1 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_0_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_1_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_0_2_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_1_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_1_2_beta2 + i);
    // }
    // for (int i = 0; i < n; i++)
    // {
    //   myfile << *(grad2_2_2_beta2 + i);
    // }

    myfile.close();
  }
  test = 0;
}