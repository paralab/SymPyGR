#include <fstream>
ofstream gpuTestFile("ouput_cuda.txt");
if ((gpuTestFile.is_open()))
{
  for (int i = 0; i < total_points; i++)
  {
    gpuTestFile << *(grad_0_alpha_cpu + i) << "\n";
  }
  gpuTestFile.close();
}
free(grad_0_alpha_cpu);