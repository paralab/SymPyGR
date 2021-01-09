for (unsigned int r=0; r < iter; r++)
{
    for(unsigned int k=3; k < nz-3; k++)
    {
        for(unsigned int j=3; j < ny-3; j++)
        {
            #pragma novector         
            for(unsigned int i=3; i < nx-3; i++)
            {
                const double x  =  i * hx;
                const double y  =  j * hy;
                const double z  =  k * hz;

                const double r_coord =sqrt(x*x + y*y + z*z);
                const double sigma=0.4;
                double eta=bssn::ETA_CONST;
                if (r_coord >= bssn::ETA_R0) {
                    eta *= pow( (bssn::ETA_R0/r_coord), bssn::ETA_DAMPING_EXP);
                }

                const unsigned int pp = k*ny*nx + j* nx + i;

                #include "../../src/bssneqs_eta_const_standard_gauge.cpp"
            }
        }
    }
    //prevent compiler from optimizing out the iter loop
    #pragma novector
    for (unsigned int pp=0; pp < NN*bssn::BSSN_NUM_VARS; pp+=10)
    {
        unzipOut[pp] += 1e-6;
        unzipIn[pp] +=1e-6;
    }
}