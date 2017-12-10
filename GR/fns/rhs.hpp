#ifndef HAVE_RHS_H
#define HAVE_RHS_H

void bssn_rhs(coll_point_t *point, const int gen);
void sommerfeld_boundary(coll_point_t *point, const int gen);
void enforce_bssn_constraints(double *upt);


#endif
