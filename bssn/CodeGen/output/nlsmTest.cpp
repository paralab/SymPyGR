//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Header file for the GR simulation.
*/
//Run: 'cog -o output/nlsmTest.cpp nlsm.cpp.in' to generate cpp file
//

#include "nlsm.h"
#include "nlsmUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rk4nlsm.h"
#include "octUtils.h"

int main (int argc, char** argv)
{
    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

	/*[[[cog
	import nlsmParams as nlsmParams

	nlsmParams.parameters.writeJson("output/nlsm.par.json")
	nlsmParams.parameters.writeCpp("output/nlsmParams.h", "output/nlsmParams.cpp", nlsmParams.namespace)

	]]]*/
	//[[[end]]]

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    nlsm::timer::initFlops();

    nlsm::timer::total_runtime.start();

    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    nlsm::readParamFile(argv[1],comm);

    if(rank==1|| npes==1)
    {
		std::cout << "parameters read: " << std::endl;

		/*[[[cog
		import cog
		import nlsmParams as nlsmParams

		for param in nlsmParams.parameters.values():
			if param.isComment():
				continue

			if param.arraySize is not None:
				cog.out("std::cout<<YLW<<\"\\t{0} : (\"<<".format(param.id))
				#may want to reconsider this formatting, cout arrays on multiple lines
				for i in range(len(param.value)):
					if i > 0:
						cog.out("<<\", \"<<")
					cog.out("{0}::{1}[{2}]".format(nlsmParams.namespace, param.id, i))
				cog.outl("<<\" )\"<<NRM<<std::endl;")
			else:
				cog.outl("std::cout<<YLW<<\"\\t{0} :\"<<{1}::{0}<<NRM<<std::endl;".format(param.id, nlsmParams.namespace))

		]]]*/
		std::cout<<YLW<<"\tDENDRO_VERSION :"<<nlsm::DENDRO_VERSION<<NRM<<std::endl;
		std::cout<<YLW<<"\tNUM_VARS :"<<nlsm::NUM_VARS<<NRM<<std::endl;
		std::cout<<YLW<<"\tELE_ORDER :"<<nlsm::ELE_ORDER<<NRM<<std::endl;
		std::cout<<YLW<<"\tRESTORE_SOLVER :"<<nlsm::RESTORE_SOLVER<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_IO_OUTPUT_FREQ :"<<nlsm::NLSM_IO_OUTPUT_FREQ<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_REMESH_TEST_FREQ :"<<nlsm::NLSM_REMESH_TEST_FREQ<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_CHECKPT_FREQ :"<<nlsm::NLSM_CHECKPT_FREQ<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_IO_OUTPUT_GAP :"<<nlsm::NLSM_IO_OUTPUT_GAP<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_VTU_FILE_PREFIX :"<<nlsm::NLSM_VTU_FILE_PREFIX<<NRM<<std::endl;
		std::cout<<YLW<<"\tCHKPT_FILE_PREFIX :"<<nlsm::CHKPT_FILE_PREFIX<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_PROFILE_FILE_PREFIX :"<<nlsm::NLSM_PROFILE_FILE_PREFIX<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_NUM_EVOL_VARS_VTU_OUTPUT :"<<nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_VTU_OUTPUT_EVOL_INDICES : ("<<nlsm::NLSM_VTU_OUTPUT_EVOL_INDICES[0]<<", "<<nlsm::NLSM_VTU_OUTPUT_EVOL_INDICES[1]<<" )"<<NRM<<std::endl;
		std::cout<<YLW<<"\tDENDRO_GRAIN_SZ :"<<nlsm::DENDRO_GRAIN_SZ<<NRM<<std::endl;
		std::cout<<YLW<<"\tASYNC_COMM_K :"<<nlsm::ASYNC_COMM_K<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_DENDRO_AMR_FAC :"<<nlsm::NLSM_DENDRO_AMR_FAC<<NRM<<std::endl;
		std::cout<<YLW<<"\tLOAD_IMB_TOL :"<<nlsm::LOAD_IMB_TOL<<NRM<<std::endl;
		std::cout<<YLW<<"\tDIM :"<<nlsm::DIM<<NRM<<std::endl;
		std::cout<<YLW<<"\tMAXDEPTH :"<<nlsm::MAXDEPTH<<NRM<<std::endl;
		std::cout<<YLW<<"\tSPLIT_FIX :"<<nlsm::SPLIT_FIX<<NRM<<std::endl;
		std::cout<<YLW<<"\tWAVELET_TOL :"<<nlsm::WAVELET_TOL<<NRM<<std::endl;
		std::cout<<YLW<<"\tNUM_REFINE_VARS :"<<nlsm::NUM_REFINE_VARS<<NRM<<std::endl;
		std::cout<<YLW<<"\tREFINE_VARIABLE_INDICES : ("<<nlsm::REFINE_VARIABLE_INDICES[0]<<", "<<nlsm::REFINE_VARIABLE_INDICES[1]<<" )"<<NRM<<std::endl;
		std::cout<<YLW<<"\tGRID_MIN_X :"<<nlsm::GRID_MIN_X<<NRM<<std::endl;
		std::cout<<YLW<<"\tGRID_MAX_X :"<<nlsm::GRID_MAX_X<<NRM<<std::endl;
		std::cout<<YLW<<"\tGRID_MIN_Y :"<<nlsm::GRID_MIN_Y<<NRM<<std::endl;
		std::cout<<YLW<<"\tGRID_MAX_Y :"<<nlsm::GRID_MAX_Y<<NRM<<std::endl;
		std::cout<<YLW<<"\tGRID_MIN_Z :"<<nlsm::GRID_MIN_Z<<NRM<<std::endl;
		std::cout<<YLW<<"\tGRID_MAX_Z :"<<nlsm::GRID_MAX_Z<<NRM<<std::endl;
		std::cout<<YLW<<"\tRK45_TIME_BEGIN :"<<nlsm::RK45_TIME_BEGIN<<NRM<<std::endl;
		std::cout<<YLW<<"\tRK45_TIME_END :"<<nlsm::RK45_TIME_END<<NRM<<std::endl;
		std::cout<<YLW<<"\tCFL_FACTOR :"<<nlsm::CFL_FACTOR<<NRM<<std::endl;
		std::cout<<YLW<<"\tCOMPD_MIN : ("<<nlsm::COMPD_MIN[0]<<", "<<nlsm::COMPD_MIN[1]<<", "<<nlsm::COMPD_MIN[2]<<" )"<<NRM<<std::endl;
		std::cout<<YLW<<"\tCOMPD_MAX : ("<<nlsm::COMPD_MAX[0]<<", "<<nlsm::COMPD_MAX[1]<<", "<<nlsm::COMPD_MAX[2]<<" )"<<NRM<<std::endl;
		std::cout<<YLW<<"\tRK45_TIME_STEP_SIZE :"<<nlsm::RK45_TIME_STEP_SIZE<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_RK45_DESIRED_TOL :"<<nlsm::NLSM_RK45_DESIRED_TOL<<NRM<<std::endl;
		std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<nlsm::KO_DISS_SIGMA<<NRM<<std::endl;
		std::cout<<YLW<<"\tENABLE_BLOCK_ADAPTIVITY :"<<nlsm::ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
		std::cout<<YLW<<"\tBLK_MIN_X :"<<nlsm::BLK_MIN_X<<NRM<<std::endl;
		std::cout<<YLW<<"\tBLK_MIN_Y :"<<nlsm::BLK_MIN_Y<<NRM<<std::endl;
		std::cout<<YLW<<"\tBLK_MIN_Z :"<<nlsm::BLK_MIN_Z<<NRM<<std::endl;
		std::cout<<YLW<<"\tBLK_MAX_X :"<<nlsm::BLK_MAX_X<<NRM<<std::endl;
		std::cout<<YLW<<"\tBLK_MAX_Y :"<<nlsm::BLK_MAX_Y<<NRM<<std::endl;
		std::cout<<YLW<<"\tBLK_MAX_Z :"<<nlsm::BLK_MAX_Z<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_TYPE :"<<nlsm::NLSM_ID_TYPE<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_AMP1 :"<<nlsm::NLSM_ID_AMP1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_R1 :"<<nlsm::NLSM_ID_R1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_DELTA1 :"<<nlsm::NLSM_ID_DELTA1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_XC1 :"<<nlsm::NLSM_ID_XC1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_YC1 :"<<nlsm::NLSM_ID_YC1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_ZC1 :"<<nlsm::NLSM_ID_ZC1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_EPSX1 :"<<nlsm::NLSM_ID_EPSX1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_EPSY1 :"<<nlsm::NLSM_ID_EPSY1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_NU1 :"<<nlsm::NLSM_ID_NU1<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_OMEGA :"<<nlsm::NLSM_ID_OMEGA<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_AMP2 :"<<nlsm::NLSM_ID_AMP2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_R2 :"<<nlsm::NLSM_ID_R2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_DELTA2 :"<<nlsm::NLSM_ID_DELTA2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_XC2 :"<<nlsm::NLSM_ID_XC2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_YC2 :"<<nlsm::NLSM_ID_YC2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_ZC2 :"<<nlsm::NLSM_ID_ZC2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_EPSX2 :"<<nlsm::NLSM_ID_EPSX2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_EPSY2 :"<<nlsm::NLSM_ID_EPSY2<<NRM<<std::endl;
		std::cout<<YLW<<"\tNLSM_ID_NU2 :"<<nlsm::NLSM_ID_NU2<<NRM<<std::endl;
		//[[[end]]]
    }

	/*[[[cog
	import cog
	import nlsmParams as nlsmParams

	cog.outl("_InitializeHcurve({0}::DIM);".format(nlsmParams.namespace))
	cog.outl("m_uiMaxDepth={0}::MAXDEPTH;".format(nlsmParams.namespace))

	cog.outl("if({0}::NUM_VARS%{0}::ASYNC_COMM_K!=0)".format(nlsmParams.namespace))
	cog.outl("{")
	cog.outl("\tif(!rank)")
	cog.outl("\t\tstd::cout<<\"[overlap communication error]: total NUM_VARS: \"<<{0}::NUM_VARS<<\" is not divisable by ASYNC_COMM_K: \"<<{0}::ASYNC_COMM_K<<std::endl;".format(nlsmParams.namespace))
	cog.outl("\texit(0);")
	cog.outl("}")

	]]]*/
	_InitializeHcurve(nlsm::DIM);
	m_uiMaxDepth=nlsm::MAXDEPTH;
	if(nlsm::NUM_VARS%nlsm::ASYNC_COMM_K!=0)
	{
		if(!rank)
			std::cout<<"[overlap communication error]: total NUM_VARS: "<<nlsm::NUM_VARS<<" is not divisable by ASYNC_COMM_K: "<<nlsm::ASYNC_COMM_K<<std::endl;
		exit(0);
	}
	//[[[end]]]

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){nlsm::initData(x,y,z,var);};
    std::function<void(double,double,double,double,double*)> u_x_t=[](double x,double y,double z,double t,double*var){nlsm::analyticalSol(x,y,z,t,var);};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){nlsm::KerrSchildData(x,y,z,var);};

	/*[[[cog
	import cog
	import nlsmParams as nlsmParams

	cog.outl("const unsigned int interpVars={0}::NUM_VARS;".format(nlsmParams.namespace))
	cog.outl("unsigned int varIndex[interpVars];")
	cog.outl("for(unsigned int i=0;i<{0}::NUM_VARS;i++)".format(nlsmParams.namespace))
	cog.outl("\tvarIndex[i]=i;")

	cog.outl("")

	cog.outl("DendroIntL localSz,globalSz;")
	cog.outl("double t_stat;")
	cog.outl("double t_stat_g[3];")
	cog.outl("{0}::timer::t_f2o.start();".format(nlsmParams.namespace))
	cog.outl("if({0}::ENABLE_BLOCK_ADAPTIVITY)".format(nlsmParams.namespace))
	cog.outl("{")

	cog.outl("\tif(!rank) std::cout<<YLW<<\"Using block adaptive mesh. AMR disabled \"<<NRM<<std::endl;")
	cog.outl("")

	cog.outl("\tconst Point pt_min({0}::BLK_MIN_X,{0}::BLK_MIN_Y,{0}::BLK_MIN_Z);".format(nlsmParams.namespace))
	cog.outl("\tconst Point pt_max({0}::BLK_MAX_X,{0}::BLK_MAX_Y,{0}::BLK_MAX_Z);".format(nlsmParams.namespace))

	cog.outl("\t{0}::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);".format(nlsmParams.namespace))
	cog.outl("}")
	cog.outl("else")
	cog.outl("{")

	cog.outl("\tif(!rank) std::cout<<YLW<<\"Using function2Octree. AMR enabled \"<<NRM<<std::endl;")
	cog.outl("")

	cog.outl("\tfunction2Octree(f_init,{0}::NUM_VARS,{0}::REFINE_VARIABLE_INDICES,{0}::NUM_REFINE_VARS,tmpNodes,m_uiMaxDepth,{0}::WAVELET_TOL,{0}::ELE_ORDER,comm);".format(nlsmParams.namespace))
	cog.outl("\tstd::cout<<\"f2o else end\"<<std::endl;")

	cog.outl("}")
	cog.outl("{0}::timer::t_f2o.stop();".format(nlsmParams.namespace))
	]]]*/
	const unsigned int interpVars=nlsm::NUM_VARS;
	unsigned int varIndex[interpVars];
	for(unsigned int i=0;i<nlsm::NUM_VARS;i++)
		varIndex[i]=i;

	DendroIntL localSz,globalSz;
	double t_stat;
	double t_stat_g[3];
	nlsm::timer::t_f2o.start();
	if(nlsm::ENABLE_BLOCK_ADAPTIVITY)
	{
		if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;

		const Point pt_min(nlsm::BLK_MIN_X,nlsm::BLK_MIN_Y,nlsm::BLK_MIN_Z);
		const Point pt_max(nlsm::BLK_MAX_X,nlsm::BLK_MAX_Y,nlsm::BLK_MAX_Z);
		nlsm::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
	}
	else
	{
		if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;

		function2Octree(f_init,nlsm::NUM_VARS,nlsm::REFINE_VARIABLE_INDICES,nlsm::NUM_REFINE_VARS,tmpNodes,m_uiMaxDepth,nlsm::WAVELET_TOL,nlsm::ELE_ORDER,comm);
		std::cout<<"f2o else end"<<std::endl;
	}
	nlsm::timer::t_f2o.stop();
	//[[[end]]]

    t_stat=nlsm::timer::t_f2o.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);

	/*[[[cog
	import cog
	import nlsmParams as nlsmParams

	cog.outl("const unsigned int grainSz={0}::DENDRO_GRAIN_SZ;".format(nlsmParams.namespace))

	]]]*/
	const unsigned int grainSz=nlsm::DENDRO_GRAIN_SZ;
	//[[[end]]]

    bool isActive;
    MPI_Comm commActive;
    const int p_npes_prev=binOp::getPrevHighestPowerOfTwo((globalSz/grainSz));
    const int p_npes_next=binOp::getNextHighestPowerOfTwo((globalSz/grainSz));

    int p_npes=globalSz/grainSz;
    (std::abs(p_npes_prev-p_npes)<=std::abs(p_npes_next-p_npes)) ? p_npes=p_npes_prev : p_npes=p_npes_next;

    if(p_npes>npes) p_npes=npes;
    // quick fix to enforce the npes>=2 for any given grain size.
    if(p_npes<=1 && npes>1) p_npes=2;

    if(p_npes==npes)
    {
        MPI_Comm_dup(comm,&commActive);
        isActive=true;

    }
	else
    {
        //isActive=(rank*grainSz<globalSz);
        isActive=isRankSelected(npes,rank,p_npes);
        par::splitComm2way(isActive,&commActive,comm);

    }

	/*[[[cog
	import cog
	import nlsmParams as nlsmParams

	cog.outl("shrinkOrExpandOctree(tmpNodes,{0}::LOAD_IMB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);".format(nlsmParams.namespace))

	]]]*/
	shrinkOrExpandOctree(tmpNodes,nlsm::LOAD_IMB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);
	//[[[end]]]

    if(!isActive)
        if(tmpNodes.size()!=0)
            std::cout<<" rank_g: "<<rank<<" isActive: "<<isActive<<" f2O octants: "<<tmpNodes.size()<<std::endl;

    std::vector<ot::TreeNode> balOct;
    localSz=0;
    if(isActive)
    {
        int rank_active,npes_active;

        MPI_Comm_size(commActive,&npes_active);
        MPI_Comm_rank(commActive,&rank_active);

        if(!rank_active) std::cout<<"[MPI_COMM_SWITCH]: "<<npes_active<<std::endl;

		/*[[[cog
		import cog
		import nlsmParams as nlsmParams

		cog.outl("ot::TreeNode root({0}::DIM,{0}::MAXDEPTH);".format(nlsmParams.namespace))
		cog.outl("std::vector<ot::TreeNode> tmpVec;")
		cog.outl("{0}::timer::t_cons.start();".format(nlsmParams.namespace))
		cog.outl("")

		cog.outl("SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,{0}::LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,{0}::SPLIT_FIX,commActive);".format(nlsmParams.namespace))
		cog.outl("std::swap(tmpNodes,tmpVec);")
		cog.outl("tmpVec.clear();")
		cog.outl("")

		cog.outl("SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,{0}::LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,{0}::SPLIT_FIX,commActive);".format(nlsmParams.namespace))
		cog.outl("std::swap(tmpNodes,tmpVec);")
		cog.outl("tmpVec.clear();")

		]]]*/
		ot::TreeNode root(nlsm::DIM,nlsm::MAXDEPTH);
		std::vector<ot::TreeNode> tmpVec;
		nlsm::timer::t_cons.start();

		SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,nlsm::LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,nlsm::SPLIT_FIX,commActive);
		std::swap(tmpNodes,tmpVec);
		tmpVec.clear();

		SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,nlsm::LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,nlsm::SPLIT_FIX,commActive);
		std::swap(tmpNodes,tmpVec);
		tmpVec.clear();
		//[[[end]]]

        nlsm::timer::t_cons.stop();
        t_stat=nlsm::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;

		/*[[[cog
		import cog
		import nlsmParams as nlsmParams

		cog.outl("{0}::timer::t_bal.start();".format(nlsmParams.namespace))
		cog.outl("")
		cog.outl("SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,{0}::LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,{0}::SPLIT_FIX,commActive);".format(nlsmParams.namespace))
		cog.outl("tmpNodes.clear();")

		cog.outl("")
		cog.outl("{0}::timer::t_bal.stop();".format(nlsmParams.namespace))
		cog.outl("t_stat={0}::timer::t_bal.seconds;".format(nlsmParams.namespace))

		]]]*/
		nlsm::timer::t_bal.start();

		SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,nlsm::LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,nlsm::SPLIT_FIX,commActive);
		tmpNodes.clear();

		nlsm::timer::t_bal.stop();
		t_stat=nlsm::timer::t_bal.seconds;
		//[[[end]]]

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        if(!rank_active) std::cout<<GRN<<" 2:1 balancing max (s): "<<t_stat_g[2]<<NRM<<std::endl;
        localSz=balOct.size();
    }
    MPI_Comm_free(&commActive);

    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;

	/*[[[cog
	import cog
	import nlsmParams as nlsmParams

	cog.outl("{0}::timer::t_mesh.start();".format(nlsmParams.namespace))
	cog.outl("")

	cog.outl("ot::Mesh * mesh=new ot::Mesh(balOct,1,{0}::ELE_ORDER,comm,true,ot::SM_TYPE::FDM,{0}::DENDRO_GRAIN_SZ,{0}::LOAD_IMB_TOL,{0}::SPLIT_FIX);".format(nlsmParams.namespace))

	cog.outl("")
	cog.outl("{0}::timer::t_mesh.stop();".format(nlsmParams.namespace))
	cog.outl("t_stat={0}::timer::t_mesh.seconds;".format(nlsmParams.namespace))

	]]]*/
	nlsm::timer::t_mesh.start();

	ot::Mesh * mesh=new ot::Mesh(balOct,1,nlsm::ELE_ORDER,comm,true,ot::SM_TYPE::FDM,nlsm::DENDRO_GRAIN_SZ,nlsm::LOAD_IMB_TOL,nlsm::SPLIT_FIX);

	nlsm::timer::t_mesh.stop();
	t_stat=nlsm::timer::t_mesh.seconds;
	//[[[end]]]

    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=mesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
    if(!rank) std::cout<<GRN<<" # of CG nodes (vertices) : "<<globalSz<<NRM<<std::endl;
    if(!rank)
    {
        std::cout<< GRN<<"Mesh generation time (max): "<<t_stat_g[2]<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;
    }

	/*[[[cog
	import cog
	import nlsmParams as nlsmParams

	cog.outl("ode::solver::RK4_NLSM rk_nlsm(mesh,{0}::RK45_TIME_BEGIN,{0}::RK45_TIME_END,{0}::RK45_TIME_STEP_SIZE);".format(nlsmParams.namespace))
	######ode::solver::RK3_NLSM rk_nlsm(mesh,nlsm::NLSM_RK45_TIME_BEGIN,nlsm::NLSM_RK45_TIME_END,nlsm::NLSM_RK45_TIME_STEP_SIZE);

	cog.outl("if({0}::RESTORE_SOLVER==1)".format(nlsmParams.namespace))
	cog.outl("\trk_nlsm.restoreCheckPoint({0}::NLSM_CHKPT_FILE_PREFIX.c_str(),comm);".format("nlsmParams.namespace"))
	cog.outl("")

	cog.outl("{0}::timer::t_rkSolve.start();".format(nlsmParams.namespace))
	cog.outl("rk_nlsm.rkSolve();")
	cog.outl("{0}::timer::t_rkSolve.stop()".format(nlsmParams.namespace))
	cog.outl("{0}::timer::total_runtime.stop();".format(nlsmParams.namespace))

	]]]*/
	ode::solver::RK4_NLSM rk_nlsm(mesh,nlsm::RK45_TIME_BEGIN,nlsm::RK45_TIME_END,nlsm::RK45_TIME_STEP_SIZE);
	if(nlsm::RESTORE_SOLVER==1)
		rk_nlsm.restoreCheckPoint(nlsmParams.namespace::NLSM_CHKPT_FILE_PREFIX.c_str(),comm);

	nlsm::timer::t_rkSolve.start();
	rk_nlsm.rkSolve();
	nlsm::timer::t_rkSolve.stop()
	nlsm::timer::total_runtime.stop();
	//[[[end]]]

    rk_nlsm.freeMesh();
    //nlsm::timer::profileInfo(nlsm::NLSM_PROFILE_FILE_PREFIX.c_str(),mesh);
    //delete mesh;
    MPI_Finalize();

    return 0;
}
