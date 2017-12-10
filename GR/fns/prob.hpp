#ifndef HAWAII_PROB_H
#define HAWAII_PROB_H

#include "param.h"
#include "types.h"

// Note: 
// 1. Points stored in the hash table might be visited more than once within a
// particular right-hande side computation stage. Care need to be taken to avoid
// redundant computation. 
// 
// 2. There are two attributes requring attention for the nonessential points
// involved in the right-hand side computation. (i) Between different right-hand
// side computation, the time stamp can be outdated because time stamp is
// incremented each time the right-hand side is computed. (ii) Within each
// right-hand side computation, the stage can be outdated. 
// 
// 3. The correct time stamp of the derivative stencil affects the correctness
// of the right-hand side computation. Inside status_update_helper, the time
// stamp is checked before applying the time integrator for the new time
// step. If the integrator involves only one generation, like Forward Euler
// integrator, the check performed in status_update_helper is
// sufficient. However, if the integrator involves more than one generation,
// such as RK4, the time stamp of the derivative stencil needs to be checked
// whenever a new generation of the right-hand side computation starts. This
// additional check needs to be incoporated in the rhs_stage1 function. 

void problem_init(void); 

void initial_ranges(double *ranges);

void initial_condition(const coord_t *coord, double *u); 
void init_data_puncture(const coord_t *coord, double *u); 
void init_data_kerrschild(const coord_t *coord, double *u); 
void init_data_fake(const coord_t *coord, double *u); 

void compute_rhs(const double t, const int gen); 

double get_local_dt(const coll_point_t *point); 

#endif 
