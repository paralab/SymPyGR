If you are executing GPU code with out test(vaerifying the output against CPU) followings are the dependecies,
    bssn/include/def.h
    bssn/include/utils.h

Otherwise,
    bssn/include/def.h
    bssn/include/utils.h
    bssn/include/rhs.h

* Default configuration of execution in GPU set to HYBRID. You can change it from following,
    cuda_gr/include/cudaBSSN.h

* In order to ran output verification against CPU, change the status of ENABLE_CUDA_TEST as follows,
    option(ENABLE_CUDA_TEST "enable test execution" ON)