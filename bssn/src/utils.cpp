//
// Created by milinda on 5/3/18.
//

#include "utils.h"

void initial_data(double *u, const double *xi)
{


    const double pi = acos(-1.0);
    const double f1 = 31.0/17.0;
    const double f2 = 37.0/11.0;

    double x=xi[0];
    double y=xi[1];
    double z=xi[2];

    u[VAR::U_ALPHA] = 1.0 - 0.25*sin(f1*x);

    u[VAR::U_BETA0] = 4.0/17.0*sin(x)*cos(z);
    u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
    u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

    u[VAR::U_B0] = 31.0*x*cos(f1*z+y);
    u[VAR::U_B1] = 7.0*y*sin(f1*x+y) + 3.0*cos(z);
    u[VAR::U_B2] = 5.0*z*cos(f1*x+y) + 7.0*sin(z+y+x) + 1.0;

    u[VAR::U_GT0] = 5.0*cos(x)/(10.0*sin(x+z)+26.0-1.0*cos(x*z)*cos(x));
    u[VAR::U_GT1] = -5.0*sin(y)/(25.0+10.0*cos(y+z)+cos(y)*cos(y*z));
    u[VAR::U_GT2] = -5.0*sin(z)/(25.0+10.0*cos(y+x)+cos(y*x)*cos(z));

    u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));
    //u[F_CHI][pp] = 2.0;

    u[VAR::U_SYMGT0] = 1.00+0.2*sin(x+z)*cos(y);
    u[VAR::U_SYMGT3] = 1.00+0.2*cos(y)*cos(z+ x);
    u[VAR::U_SYMGT5] = 1.00 / ( u[VAR::U_SYMGT0] +  u[VAR::U_SYMGT3] );
    u[VAR::U_SYMGT1] = 0.7*cos(x*x + y*y);
    u[VAR::U_SYMGT2] = 0.3*sin(z)*cos(x);
    u[VAR::U_SYMGT4] = -0.5*sin(x*x)*cos(y)*cos(z);

    u[VAR::U_K] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
                      +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
                      +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
                       *exp(-4.0*cos(x)*sin(y))*cos(z);

    u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
    u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);

    u[VAR::U_SYMAT3] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);

    u[VAR::U_SYMAT5] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));

}