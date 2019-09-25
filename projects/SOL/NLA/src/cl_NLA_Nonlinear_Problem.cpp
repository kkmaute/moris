/*
 * cl_NLA_Nonlinear_Problem.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_SOL_Warehouse.hpp"

#include <ctime>

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Enums.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

Nonlinear_Problem::Nonlinear_Problem(       SOL_Warehouse      * aNonlinDatabase,
                                            Solver_Interface   * aSolverInterface,
                                            Dist_Vector        * aFullVector,
                                      const moris::sint          aNonlinearSolverManagerIndex,
                                      const bool                 aBuildLinerSystemFlag,
                                      const enum MapType         aMapType) :     mFullVector( aFullVector ),
                                                                                 mBuildLinerSystemFlag( aBuildLinerSystemFlag ),
                                                                                 mMapType( aMapType ),
                                                                                 mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
    mSolverInterface = aSolverInterface;

    if( mMapType == MapType::Petsc )
    {
        // Initialize petsc solvers
        PetscInitializeNoArguments();
    }

    // delete pointers if they already exist
    this->delete_pointers();

    // create solver factory
    Solver_Factory  tSolFactory;

    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory( mMapType );

    // create map object FIXME ask liner problem for map
    mMap = tMatFactory.create_map( aSolverInterface->get_max_num_global_dofs(),
                                   aSolverInterface->get_my_local_global_map(),
                                   aSolverInterface->get_constr_dof());

    // create map object FIXME ask liner problem for map
    mMapFull = tMatFactory.create_map( aSolverInterface->get_my_local_global_overlapping_map() );


    // create solver object
    if ( mBuildLinerSystemFlag )
    {
        mLinearProblem = tSolFactory.create_linear_system( aSolverInterface,
                                                           mMap,
                                                           mMapFull,
                                                           mMapType );
    }

    //---------------------------arc-length vectors---------------------------------
    mFext          = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );

    mJacVals       = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mJacVals0      = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );

    mDTildeVec     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mDTilde0Vec    = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );

    mDK            = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mDSolve        = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mDSolveNMinus1 = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mDSolveNMinus2 = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );

    mGlobalRHS     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );

    mDFArcDDeltaD  = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );

    mDelLamNum     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mDelLamDen     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mDeltaD        = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    mdeltaD        = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FREE );
    //---------------------------arc-length matrices--------------------------------
    mJacobian  = tMatFactory.create_matrix( aSolverInterface, mMap );

    //------------------------------------------------------------------------------
}

Nonlinear_Problem::Nonlinear_Problem(       Solver_Interface * aSolverInterface,
                                      const moris::sint        aNonlinearSolverManagerIndex,
                                      const bool               aBuildLinerSystemFlag,
                                      const enum MapType       aMapType) :     mBuildLinerSystemFlag( aBuildLinerSystemFlag ),
                                                                               mMapType( aMapType ),
                                                                               mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
    mSolverInterface = aSolverInterface;

    if( mMapType == MapType::Petsc )
    {
        // Initialize petsc solvers
        PetscInitializeNoArguments();
    }

    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory( mMapType );

    // create map object FIXME ask liner problem for map
    mMap = tMatFactory.create_map( aSolverInterface->get_max_num_global_dofs(),
                                   aSolverInterface->get_my_local_global_map(),
                                   aSolverInterface->get_constr_dof(),
                                   aSolverInterface->get_my_local_global_overlapping_map());

    // full vector
    mFullVector = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mDummyFullVector = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );       // FIXME delete
    mDummyFullVector->vec_put_scalar( 0.0 );
    aSolverInterface->set_solution_vector_prev_time_step(mDummyFullVector);

    //---------------------------arc-length vectors---------------------------------
    mFext          = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mJacVals       = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mJacVals0      = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mDTildeVec     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mDTilde0Vec    = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mDK            = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mDSolve        = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mDSolveNMinus1 = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mDSolveNMinus2 = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mGlobalRHS     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mDFArcDDeltaD  = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mDelLamNum     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mDelLamDen     = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mDeltaD        = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mdeltaD        = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    //---------------------------arc-length matrices--------------------------------
    mJacobian  = tMatFactory.create_matrix( aSolverInterface, mMap );

    //------------------------------------------------------------------------------

//    mFullForDiag = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mFullVector->vec_put_scalar( 0.0 );

    // delete pointers if they already exist
    this->delete_pointers();

    // create solver factory
    Solver_Factory  tSolFactory;

    // create solver object
    if ( mBuildLinerSystemFlag )
    {
        MORIS_LOG_INFO( "Build linear problem with index %-5i \n", mNonlinearSolverManagerIndex );

        mLinearProblem = tSolFactory.create_linear_system( aSolverInterface, mMapType );
    }

    // set flag that interface has been set
    mIsMasterSystem = true;
}

void Nonlinear_Problem::set_interface( Solver_Interface * aSolverInterface )
{
}

Nonlinear_Problem::~Nonlinear_Problem()
{
    this->delete_pointers();

    if( mMap != nullptr )
    {
        delete( mMap );
    }

    if( mMapFull != nullptr )
    {
        delete( mMapFull );
    }

    if( mIsMasterSystem )
    {
        delete( mFullVector );
        delete( mFullVector );

        delete( mJacVals );
        delete( mJacVals0 );
        delete( mDTildeVec );
        delete( mDTilde0Vec );
        delete( mDK );
        delete( mDSolve );
        delete( mDSolveNMinus1 );
        delete( mDSolveNMinus2 );
        delete( mGlobalRHS );
        delete( mDFArcDDeltaD );
        delete( mDelLamNum );
        delete( mDelLamDen );
        delete( mDeltaD );
        delete( mdeltaD );
        delete( mFext );

        delete( mJacobian );
    }

    if ( mMapType == MapType::Petsc)
    {
        PetscFinalize();
    }
}

void Nonlinear_Problem::delete_pointers()
{
    if( mLinearProblem != nullptr )
    {
        delete( mLinearProblem );
    }
}

void Nonlinear_Problem::build_linearized_problem( const bool & aRebuildJacobian,
                                                  const sint aNonLinearIt )
{
    // Set VectorFreeSol and LHS
    mLinearProblem->set_free_solver_LHS( mFullVector );

    this->print_sol_vec( aNonLinearIt );

    if( aRebuildJacobian )
    {
        mLinearProblem->assemble_jacobian( mFullVector );
    }

    mLinearProblem->assemble_residual( mFullVector );
}

void Nonlinear_Problem::build_linearized_problem( const bool & aRebuildJacobian,
                                                  const sint aNonLinearIt,
                                                  const sint aRestart )
{
    delete( mFullVector );

    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory;
    mFullVector = tMatFactory.create_vector();

    this->restart_from_sol_vec( aRestart );

    // Set VectorFreeSol and LHS
    mLinearProblem->set_free_solver_LHS( mFullVector );

    if( aRebuildJacobian )
    {
        mLinearProblem->assemble_jacobian( mFullVector );
    }

    mLinearProblem->assemble_residual( mFullVector );
}

Dist_Vector * Nonlinear_Problem::get_full_vector()
{
    return mFullVector;
}

void Nonlinear_Problem::extract_my_values( const moris::uint         & aNumIndices,
                                       const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                       const moris::uint             & aBlockRowOffsets,
                                             moris::Matrix< DDRMat > & LHSValues )
{
    mFullVector->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
}

void Nonlinear_Problem::print_sol_vec( const sint aNonLinearIt )
{
    char NonLinNum[100];
    std::sprintf( NonLinNum, "NonLIt.%04u", aNonLinearIt );

    char SolVector[100];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector,".h5\0");

    mFullVector->save_vector_to_HDF5( SolVector );
}

void Nonlinear_Problem::restart_from_sol_vec( const sint aRestart )
{
    char NonLinNum[100];
    std::sprintf( NonLinNum, "NonLIt.%04u", aRestart );

    char SolVector[100];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector,".h5\0");

    mFullVector->read_vector_from_HDF5( SolVector );
}

//--------------------------------------------------------------------------------------------------
void Nonlinear_Problem::set_time_value( const moris::real & aLambda,
                                              moris::uint   aPos )
{
    mSolverInterface->set_time_value( aLambda, aPos );
}
