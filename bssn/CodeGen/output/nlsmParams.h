//This file was generated at 07:03PM on October 07, 2019, through the parameters object defined in cog
#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include <string.h>
#include <iostream>

namespace nlsm
{

	extern double DENDRO_VERSION;

	extern unsigned int NUM_VARS;

	//element order
	extern unsigned int ELE_ORDER;

	//Set to 1 to restore solver from a checkpoint. 0 otherwise
	extern unsigned int RESTORE_SOLVER;

	//frequency for VTU output
	extern unsigned int NLSM_IO_OUTPUT_FREQ;

	//frequency for remeshing test based on wavelets
	extern unsigned int NLSM_REMESH_TEST_FREQ;

	//frequency for checkpoint output
	extern unsigned int NLSM_CHECKPT_FREQ;

	//VTU file output gap. (Not currently used. Might be usefull in adaptive timestepping)
	extern unsigned int NLSM_IO_OUTPUT_GAP;

	//file prefix for the vtu files
	extern std::string NLSM_VTU_FILE_PREFIX;

	//file prefix for the checkpoint files
	extern std::string CHKPT_FILE_PREFIX;

	//file prefix for the intermediate profile files
	extern std::string NLSM_PROFILE_FILE_PREFIX;

	//number of variables (evolution) to output in vtu files
	extern unsigned int NLSM_NUM_EVOL_VARS_VTU_OUTPUT;

	//evolution variable ids
	extern unsigned int NLSM_VTU_OUTPUT_EVOL_INDICES[2];

	//grain size N/p , Where N number of total octants, p number of active cores
	extern unsigned int DENDRO_GRAIN_SZ;

	//variable group size for the asynchronous unzip operation
	extern unsigned int ASYNC_COMM_K;

	//dendro coarsening factor, corsent if computed wavelet tol < NLSM_DENDRO_AMR_FAC*NLSM_WAVELET_TOL 
	extern double NLSM_DENDRO_AMR_FAC;

	//dendro load imbalance tolerance for flexible partitioning
	extern double LOAD_IMB_TOL;

	//dimentionality of the octree, (meshing is supported only for 3D)
	extern unsigned int DIM;

	//maximum level of refinement of the mesh
	extern unsigned int MAXDEPTH;

	//Splitter fix value
	extern unsigned int SPLIT_FIX;

	//wavelet tolerance
	extern double WAVELET_TOL;

	//number of refinement variables
	extern unsigned int NUM_REFINE_VARS;

	//refinement variable IDs
	extern unsigned int REFINE_VARIABLE_INDICES[2];

	extern double GRID_MIN_X;

	extern double GRID_MAX_X;

	extern double GRID_MIN_Y;

	extern double GRID_MAX_Y;

	extern double GRID_MIN_Z;

	extern double GRID_MAX_Z;

	//simulation time begin
	extern unsigned int RK45_TIME_BEGIN;

	//simulation time end
	extern unsigned int RK45_TIME_END;

	//CFL Factor
	extern double CFL_FACTOR;

	//Should be an array defined from values of GRID_MIN
	extern double COMPD_MIN[3];

	//Should be an array defined from values of GRID_MAX
	extern double COMPD_MAX[3];

	//value should be calculated using the CFL Factor
	extern double RK45_TIME_STEP_SIZE;

	//used in adaptive time stepping (not currently used)
	extern double NLSM_RK45_DESIRED_TOL;

	//Kreiss-Oliger dissipation
	extern double KO_DISS_SIGMA;

	//Set to 1 disable AMR and use block adaptivity (not recomended).
	extern unsigned int ENABLE_BLOCK_ADAPTIVITY;

	extern double BLK_MIN_X;

	extern double BLK_MIN_Y;

	extern double BLK_MIN_Z;

	extern double BLK_MAX_X;

	extern double BLK_MAX_Y;

	extern double BLK_MAX_Z;

	//0-test data
	extern unsigned int NLSM_ID_TYPE;

	extern double NLSM_ID_AMP1;

	extern double NLSM_ID_R1;

	extern double NLSM_ID_DELTA1;

	extern double NLSM_ID_XC1;

	extern double NLSM_ID_YC1;

	extern double NLSM_ID_ZC1;

	extern double NLSM_ID_EPSX1;

	extern double NLSM_ID_EPSY1;

	extern double NLSM_ID_NU1;

	extern double NLSM_ID_OMEGA;

	extern double NLSM_ID_AMP2;

	extern double NLSM_ID_R2;

	extern double NLSM_ID_DELTA2;

	extern double NLSM_ID_XC2;

	extern double NLSM_ID_YC2;

	extern double NLSM_ID_ZC2;

	extern double NLSM_ID_EPSX2;

	extern double NLSM_ID_EPSY2;

	extern double NLSM_ID_NU2;

}
#endif //SFCSORTBENCH_PARAMETERS_H