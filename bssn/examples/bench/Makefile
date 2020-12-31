#CXXFLAGS += -xCOMMON-AVX512 -g -O3 -Wall -std=c++11
CXXFLAGS += -g -O3 -Wall -std=c++11 -march=native -mprefer-vector-width=512                                   

#CXXFLAGS += -g -O3 -Wall -std=c++11 -march=native -mavx2
#CHPC env
LIKWID_INC=/uufs/chpc.utah.edu/sys/installdir/likwid/5.0.1/include
LIKWID_LIB=/uufs/chpc.utah.edu/sys/installdir/likwid/5.0.1/lib
RUN_PATH=-Wl,-rpath,/uufs/chpc.utah.edu/sys/installdir/likwid/5.0.1/lib

SOURCE=a.cpp
a: a.cpp
	$(CXX) -S $(CXXFLAGS) -DN=$(n) -DR=10000000 -o a.s $<
	$(CXX) $(CXXFLAGS) -DN=$(n) -DR=10000000 -I$(LIKWID_INC) -L$(LIKWID_LIB) $(RUN_PATH) -fopenmp -DLIKWID_PERFMON -o $@ $< -llikwid

clean:
	rm a



