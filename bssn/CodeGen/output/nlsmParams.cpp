//This file was generated at 07:03PM on October 07, 2019, through the parameters object defined in cog
#include "nlsmParams.h"
namespace nlsm
{

	double DENDRO_VERSION=5.0;

	unsigned int NUM_VARS=2;

	//element order
	unsigned int ELE_ORDER=4;

	//Set to 1 to restore solver from a checkpoint. 0 otherwise
	unsigned int RESTORE_SOLVER=0;

	//frequency for VTU output
	unsigned int NLSM_IO_OUTPUT_FREQ=400;

	//frequency for remeshing test based on wavelets
	unsigned int NLSM_REMESH_TEST_FREQ=10;

	//frequency for checkpoint output
	unsigned int NLSM_CHECKPT_FREQ=10000;

	//VTU file output gap. (Not currently used. Might be usefull in adaptive timestepping)
	unsigned int NLSM_IO_OUTPUT_GAP=1;

	//file prefix for the vtu files
	std::string NLSM_VTU_FILE_PREFIX="nlsm_gr";

	//file prefix for the checkpoint files
	std::string CHKPT_FILE_PREFIX="nlsm_cp";

	//file prefix for the intermediate profile files
	std::string NLSM_PROFILE_FILE_PREFIX="nlsm_prof";

	//number of variables (evolution) to output in vtu files
	unsigned int NLSM_NUM_EVOL_VARS_VTU_OUTPUT=2;

	//evolution variable ids
	unsigned int NLSM_VTU_OUTPUT_EVOL_INDICES[2]={0,1};

	//grain size N/p , Where N number of total octants, p number of active cores
	unsigned int DENDRO_GRAIN_SZ=100;

	//variable group size for the asynchronous unzip operation
	unsigned int ASYNC_COMM_K=2;

	//dendro coarsening factor, corsent if computed wavelet tol < NLSM_DENDRO_AMR_FAC*NLSM_WAVELET_TOL 
	double NLSM_DENDRO_AMR_FAC=1.0;

	//dendro load imbalance tolerance for flexible partitioning
	double LOAD_IMB_TOL=0.1;

	//dimentionality of the octree, (meshing is supported only for 3D)
	unsigned int DIM=3;

	//maximum level of refinement of the mesh
	unsigned int MAXDEPTH=6;

	//Splitter fix value
	unsigned int SPLIT_FIX=2;

	//wavelet tolerance
	double WAVELET_TOL=0.0001;

	//number of refinement variables
	unsigned int NUM_REFINE_VARS=2;

	//refinement variable IDs
	unsigned int REFINE_VARIABLE_INDICES[2]={0,1};

	double GRID_MIN_X=-200.0;

	double GRID_MAX_X=200.0;

	double GRID_MIN_Y=-200.0;

	double GRID_MAX_Y=200.0;

	double GRID_MIN_Z=-200.0;

	double GRID_MAX_Z=200.0;

	//simulation time begin
	unsigned int RK45_TIME_BEGIN=0;

	//simulation time end
	unsigned int RK45_TIME_END=10;

	//CFL Factor
	double CFL_FACTOR=0.1;

	//Should be an array defined from values of GRID_MIN
	double COMPD_MIN[3]={-200.0,-200.0,-200.0};

	//Should be an array defined from values of GRID_MAX
	double COMPD_MAX[3]={200.0,200.0,200.0};

	//value should be calculated using the CFL Factor
	double RK45_TIME_STEP_SIZE=0.625;

	//used in adaptive time stepping (not currently used)
	double NLSM_RK45_DESIRED_TOL=0.001;

	//Kreiss-Oliger dissipation
	double KO_DISS_SIGMA=0.1;

	//Set to 1 disable AMR and use block adaptivity (not recomended).
	unsigned int ENABLE_BLOCK_ADAPTIVITY=0;

	double BLK_MIN_X=-6.0;

	double BLK_MIN_Y=-6.0;

	double BLK_MIN_Z=-6.0;

	double BLK_MAX_X=6.0;

	double BLK_MAX_Y=6.0;

	double BLK_MAX_Z=6.0;

	//0-test data
	unsigned int NLSM_ID_TYPE=0;

	double NLSM_ID_AMP1=1.3;

	double NLSM_ID_R1=8.0;

	double NLSM_ID_DELTA1=3.0;

	double NLSM_ID_XC1=0.0;

	double NLSM_ID_YC1=0.0;

	double NLSM_ID_ZC1=0.0;

	double NLSM_ID_EPSX1=0.5;

	double NLSM_ID_EPSY1=1.0;

	double NLSM_ID_NU1=0.0;

	double NLSM_ID_OMEGA=0.4;

	double NLSM_ID_AMP2=0.0;

	double NLSM_ID_R2=0.0;

	double NLSM_ID_DELTA2=3.0;

	double NLSM_ID_XC2=0.0;

	double NLSM_ID_YC2=0.0;

	double NLSM_ID_ZC2=0.0;

	double NLSM_ID_EPSX2=2.0;

	double NLSM_ID_EPSY2=1.0;

	double NLSM_ID_NU2=0.0;

}