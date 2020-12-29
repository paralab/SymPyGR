CXXFLAGS += -g -O3 -Wall -std=c++11 -march=native -mavx2

#CHPC env
LIKWID_INC=/uufs/chpc.utah.edu/sys/installdir/likwid/5.0.1/include
LIKWID_LIB=/uufs/chpc.utah.edu/sys/installdir/likwid/5.0.1/lib
RUN_PATH=-Wl,-rpath,/uufs/chpc.utah.edu/sys/installdir/likwid/5.0.1/lib

a: a.cpp
	$(CXX) $(CXXFLAGS) -DN=256 -DR=1000000 -I$(LIKWID_INC) -L$(LIKWID_LIB) $(RUN_PATH) -fopenmp -DLIKWID_PERFMON -o $@ $< -llikwid
clean:
	rm a	
