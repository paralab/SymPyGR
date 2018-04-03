#include <fstream>
ofstream cpuTestFile("output_cpu.txt");
if ((cpuTestFile.is_open()))
{
  for (int i = 0; i < total_points; i++)
  {
    cpuTestFile << *(grad_0_alpha_cpu + i) << "\n";
  }
  cpuTestFile.close();
}