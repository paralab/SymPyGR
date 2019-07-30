1. GW & multimessenger astronomy is important. 1-2 sentences (15s) for motivation

   1. simulations to get G waveforms is essential 
   2. current detections have been for roughly equal mass ratio,
      1. eXisting frameworks do not scale/ handle large mass-ratios.  $q^3$ 
   3. Our framework is scalable and fully adaptive and demonstrate upto $q=100$ and 131k cores of Titan

2. Now, let me explain why this **challenging** and highlight our main **contributions**.

   1. complexity of Einstein Equations 
      1. Curved Space time, **FD** formulation (**adaptivity**)
         1. common: block-adaptivity, **we**: WAMR+Octrees
         2. new algorithms: **TreeSearch**
      2. Difficult to write and optimize code, **portability & extensibility**
         1. age of ET/Cactus, 
         2. SymPy + Optimized generation of code. 
            1. avx, openMP & CUDA
      3. Scale of problem: **scalability**
         1. highlight size of problem, need for supercomputing resources
         2. previous runs, took order of month to simulate. 
         3. during the observation runs, there is a need to do on-demand sims. Would like this to be < 24 hrs. 

3. Methods

4. Results/ illustration
