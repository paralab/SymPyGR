## EQTEST CODE

This code is designed to verify that the BSSN equations generated for
Dendro are equivalent to the BSSN equations currently used in the had
code.

The files
  1. main
  2. rhs.cpp : The Dendro RHS
  3. bssnrhs.F90 : The had RHS

The eqtest code

  1. Generates arbitrary initial data for all functions. 
  2. The code calls the had RHS routine.
  3. The code calls the Dendro RHS routine.
  4. The code then compares all RHS values.
