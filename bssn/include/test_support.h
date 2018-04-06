/**
 * Created on: April 6, 2018
 * 		Author: Akila
 **/

#ifndef TEST_SUPPORT_H_
#define TEST_SUPPORT_H_

#include <fstream>

namespace test_file_write
{
    void writeToFile(const char filename[], double * array, int size);
    void appendToFile(const char filename[], double * array, int size);
}

#endif
