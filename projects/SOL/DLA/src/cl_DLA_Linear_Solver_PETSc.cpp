/*
 * cl_DLA_Linear_Solver_PETSc.cpp
 *
 *  Created on: Dez 12, 2018
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Solver_PETSc.hpp"
#include "cl_DLA_Preconditioner_PETSc.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "petscmat.h"

#include <string>

#include "cl_Tracer.hpp"

using namespace moris;
using namespace dla;

Linear_Solver_PETSc::Linear_Solver_PETSc()
{

    this->set_solver_parameters();
}

Linear_Solver_PETSc::Linear_Solver_PETSc( const moris::ParameterList aParameterlist )
: Linear_Solver_Algorithm( aParameterlist )
{

}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::Linear_Solver_PETSc( Linear_Problem * aLinearSystem )
{
    mLinearSystem = aLinearSystem;

    //FIXME add rest
    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::~Linear_Solver_PETSc()
{
    //KSPDestroy(&mPetscKSPProblem);
//    PCDestroy(&mpc);
}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_solver_parameters()
{
    // Create parameter list and set default values fo solver parameters

    // Set KSP type
    mParameterList.insert( "KSPType", std::string( KSPGMRES ) );

    // Set default preconditioner
    mParameterList.insert( "PCType", std::string( PCILU ) );

    // Sets maximal iters for KSP
    mParameterList.insert( "KSPMaxits", 1000 );

    // Sets KSP gmres restart
    mParameterList.insert( "KSPMGMRESRestart", 500 );

    // Sets tolerance for determining happy breakdown in GMRES, FGMRES and LGMRES
    mParameterList.insert( "KSPGMRESHapTol", 1e-10 );

    // Sets tolerance for KSP
    mParameterList.insert( "KSPTol", 1e-10 );

    // Sets the number of levels of fill to use for ILU
    mParameterList.insert( "ILUFill", 0 );

    // Sets drop tolerance for ilu
    mParameterList.insert( "ILUTol", 1e-6 );

    // Set multigrid levels
    mParameterList.insert( "MultigridLevels", 3 );
}

//----------------------------------------------------------------------------------------
moris::sint Linear_Solver_PETSc::solve_linear_system( )
{
    return 0;
}

//----------------------------------------------------------------------------------------
moris::sint Linear_Solver_PETSc::solve_linear_system(
        Linear_Problem * aLinearSystem,
        const moris::sint       aIter )
{
    Tracer tTracer( "LinearAlgorithm", "PETSc", "Solve" );

    mSolverInterface = aLinearSystem->get_solver_input();

    // Create KSP and PC
    KSPCreate( PETSC_COMM_WORLD, &mPetscKSPProblem );
    KSPGetPC( mPetscKSPProblem, &mpc );
    KSPSetOperators( mPetscKSPProblem,
                     aLinearSystem->get_matrix()->get_petsc_matrix(),
                     aLinearSystem->get_matrix()->get_petsc_matrix() );
    KSPGMRESSetOrthogonalization( mPetscKSPProblem, KSPGMRESModifiedGramSchmidtOrthogonalization );
//    KSPSetFromOptions( mPetscKSPProblem );

//    KSPSetUp( mPetscKSPProblem );

    this->set_solver_internal_parameters( );

    // build preconditiner class
    dla::Preconditioner_PETSc tPreconditioner( this );

    if ( ! strcmp(mParameterList.get< std::string >( "PCType" ).c_str(), "mg") )
    {
        tPreconditioner.build_multigrid_preconditioner( aLinearSystem );
    }
    else if ( ! strcmp(mParameterList.get< std::string >( "PCType" ).c_str(), "asm") )
    {
        // build schwarz preconditioner
        tPreconditioner.build_schwarz_preconditioner_petsc();
    }
    else if ( ! strcmp(mParameterList.get< std::string >( "PCType" ).c_str(), "mat") )
    {
        // build schwarz preconditioner
        tPreconditioner.build_schwarz_preconditioner( aLinearSystem );

        KSPSetOperators( mPetscKSPProblem,
                         aLinearSystem->get_matrix()->get_petsc_matrix(),
                         tPreconditioner.get_preconditioner_matrix()->get_petsc_matrix() );
    }

//    aLinearSystem->get_free_solver_LHS()->read_vector_from_HDF5( "Exact_Sol_petsc.h5" );
//    aLinearSystem->get_free_solver_LHS()->print();

//    aLinearSystem->get_solver_RHS()->save_vector_to_HDF5( "Res_vec.h5" );
//    aLinearSystem->get_solver_RHS()->print();

    this->set_solver_analysis_options();

    KSPSetFromOptions( mPetscKSPProblem );
    KSPSetUp( mPetscKSPProblem );

    // Solve System
    KSPSolve(
            mPetscKSPProblem,
            static_cast<Vector_PETSc*>(aLinearSystem->get_solver_RHS())->get_petsc_vector(),
            static_cast<Vector_PETSc*>(aLinearSystem->get_free_solver_LHS())->get_petsc_vector() );

    // Output
//    KSPView( mPetscKSPProblem, PETSC_VIEWER_STDOUT_WORLD );
    moris::sint Iter;
    KSPGetIterationNumber(mPetscKSPProblem, &Iter );
    std::cout<<Iter<<" Iterations"<<std::endl;


    mSolverInterface = nullptr;

    return 0;
}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_solver_analysis_options()
{
    if(true)
    {
        PetscViewer tViewerRes;
        PetscViewerCreate( PETSC_COMM_WORLD, &tViewerRes );
        PetscViewerSetType( tViewerRes, PETSCVIEWERASCII );
        PetscViewerFileSetName( tViewerRes, "Residual_Norms.txt" );

        PetscViewerAndFormat *tViewerAndFormatRes;
        PetscViewerAndFormatCreate( tViewerRes, PETSC_VIEWER_DEFAULT, &tViewerAndFormatRes );

        KSPMonitorSet( mPetscKSPProblem,
                       reinterpret_cast< int(*)( KSP, sint, real, void* ) >( KSPMonitorTrueResidualNorm ),
                       tViewerAndFormatRes,
                       NULL );
    }

    if(false)
    {
        KSPSetComputeSingularValues( mPetscKSPProblem, PETSC_TRUE );

        PetscViewer tViewerSV;
        PetscViewerCreate( PETSC_COMM_WORLD, &tViewerSV );
        PetscViewerSetType( tViewerSV, PETSCVIEWERASCII );
        PetscViewerFileSetName( tViewerSV, "Singular_Values.txt" );

        PetscViewerAndFormat *tViewerAndFormatSV;
        PetscViewerAndFormatCreate( tViewerSV, PETSC_VIEWER_DEFAULT, &tViewerAndFormatSV );

        //KSPMonitorSet(KSP ksp,PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void*),void *mctx,PetscErrorCode (*monitordestroy)(void**))
        // e.g. KSPMonitorSet(ksp1,MyKSPMonitor,NULL,0);
        // with      MyKSPMonitor(KSP,PetscInt,PetscReal,void*);

        KSPMonitorSet( mPetscKSPProblem,
                       reinterpret_cast< int(*)( KSP, sint, real, void* ) >( KSPMonitorSingularValue ),
                       tViewerAndFormatSV,
                       NULL );
    }

}

//----------------------------------------------------------------------------------------

void Linear_Solver_PETSc::set_solver_internal_parameters( )
{
        // Set KSP type
        KSPSetType( mPetscKSPProblem, mParameterList.get< std::string >( "KSPType" ).c_str() );
//        KSPSetInitialGuessNonzero( mPetscKSPProblem, PETSC_TRUE );

        // Set maxits and tolerance for ksp
        KSPSetTolerances( mPetscKSPProblem, mParameterList.get< moris::real >( "KSPTol" ), PETSC_DEFAULT, PETSC_DEFAULT, mParameterList.get< moris::sint >( "KSPMaxits" ) );

        // Set Gmres restart
        KSPGMRESSetRestart( mPetscKSPProblem, mParameterList.get< moris::sint >( "KSPMGMRESRestart" ) );

        // Sets tolerance for determining happy breakdown in GMRES, FGMRES and LGMRES.
        KSPGMRESSetHapTol( mPetscKSPProblem, mParameterList.get< moris::real >( "KSPGMRESHapTol" ) );

        // Set PC type
        PCSetType( mpc, mParameterList.get< std::string >( "PCType" ).c_str() );

        // Set levels of fill for ILU
        PCFactorSetLevels( mpc, mParameterList.get< moris::sint >( "ILUFill" ) );

        // Set drop tolerance for Ilu
        PCFactorSetDropTolerance( mpc, mParameterList.get< moris::real >( "ILUTol" ), PETSC_DEFAULT, PETSC_DEFAULT );

        PCSORSetOmega( mpc, 1 );

        PCSORSetIterations( mpc, 1 , 1 );
}


