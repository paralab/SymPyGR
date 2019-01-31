You can checkout code, and please use the master branch.

`git clone https://github.com/paralab/SymPyGR.git`

* bssn - constains all the generated code, and examples on executing the codes. 
* GR- folder contains the symbolic code generations framework which mainly uses python. 
* GR/dendro.py : Contains basic utility functions to write the Einstein equations using sympy, Code generation function using Common Sub Expression Elimination. 
* GR/bssn.py : Contains the equations written using sympy using dendro gr. 

## To build the generated code.

    The generated code is already put in to the bssn folder, go to bssn folder, make a build directory. Then `cd build` and `ccmake ../` then pres `c` to configure and then `g` to generate the make file.
    Then you should be able to make the code using, `make all -j4`

## To Run the example.

source code for the example is located at bssn/examples/include/rhsTest.h and bssn/examples/src/rhsTest.cpp
execute, `./rhsTest 0 4 10` the
par 1 : min level of the block.
par 2: max level of the block.
par 3: total number of blocks,
par 4: mode 0 for unstaged 1 to staged.  
This will create a Gaussian distribution of blocks where the levels are between min and max level, mean is (min+max)/2
This will execute both staged and unstaged version for the smooth initial data.

## To Generate the code from python.

* You should be able to generate the code, by executing the bssn.py and bssn_stages.py using python3, e.g. python3 bssn.py
* Note that the generated code is already in the bssn folder. 
