* Default configuration of execution in GPU set to HYBRID. You can change it from following,
    cuda_gr/include/cudaBSSN.h


If you are executing GPU code with out test(vaerifying the output against CPU) followings are the dependecies,
    bssn/include/def.h
    bssn/include/utils.h

Otherwise,
    bssn/include/def.h
    bssn/include/utils.h
    bssn/include/rhs.h

* In order to execute output verification against CPU, add following line to cmakelist.txt
    add_definitions(-DENABLE_CUDA_TEST)