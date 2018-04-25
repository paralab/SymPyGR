#include "test_support.h"

using namespace std;

namespace test_file_write
{
    void writeToFile(const char filename[], double * array, int size)
    {
        ofstream testFile(filename);
        if ((testFile.is_open()))
        {
            for (int i = 0; i < size; i++)
            {
                testFile << array[i] << "\n";
            }
            testFile.close();
        }
    }
    void appendToFile(const char filename[], double * array, int size)
    {
        ofstream testFile(filename, std::ios_base::app);
        if ((testFile.is_open()))
        {
            for (int i = 0; i < size; i++)
            {
                testFile << array[i] << "\n";
            }
            testFile.close();
        }
    }
} 
