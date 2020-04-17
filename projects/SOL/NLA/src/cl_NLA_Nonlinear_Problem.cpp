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
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_Tracer.hpp"
#include "cl_Tracer_Enums.hpp"


using namespace moris;
using namespace NLA;
using namespace dla;

Nonlinear_Problem::Nonlinear_Problem(       sol::SOL_Warehouse * aNonlinDatabase,
                                            Solver_Interface   * aSolverInterface,
                                            Dist_Vector        * aFullVector,
                                      const moris::sint          aNonlinearSolverManagerIndex,
                                      const bool                 aBuildLinerSystemFlag,
                                      const enum sol::MapType    aMapType) :     mFullVector( aFullVector ),
                                                                                 mBuildLinerSystemFlag( aBuildLinerSystemFlag ),
                                                                                 mMapType( aMapType ),
                                                                                 mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
    mSolverInterface = aSolverInterface;

//    if( mMapType == sol::MapType::Petsc )
//    {
//        // Initialize petsc solvers
//        PetscInitializeNoArguments();
//    }

    moris::Cell< enum MSI::Dof_Type > tRequesedDofTypes = mSolverInterface->get_requested_dof_types();

    // delete pointers if they already exist
    this->delete_pointers();

    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory( mMapType );

    // create map object FIXME ask linear problem for map
    mMap = tMatFactory.create_map( aSolverInterface->get_my_local_global_map( tRequesedDofTypes ) );

    // create map object FIXME ask linear problem for map
    mMapFull = tMatFactory.create_map( aSolverInterface->get_my_local_global_overlapping_map() );

    // create solver object
    if ( mBuildLinerSystemFlag )
    {
        // create solver factory
        Solver_Factory  tSolFactory;

        mLinearProblem = tSolFactory.create_linear_system( aSolverInterface,
                                                           mMap,
                                                           mMapFull,
                                                           mMapType );
    }
}

Nonlinear_Problem::Nonlinear_Problem(       Solver_Interface * aSolverInterface,
                                      const moris::sint        aNonlinearSolverManagerIndex,
                                      const bool               aBuildLinerSystemFlag,
                                      const enum sol::MapType  aMapType ) :     mBuildLinerSystemFlag( aBuildLinerSystemFlag ),
                                                                                mMapType( aMapType ),
                                                                                mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
    mSolverInterface = aSolverInterface;

    std::cout<<"build linear system"<<std::endl;
    if( mMapType == sol::MapType::Petsc )
    {
        std::cout<<"build linear system"<<std::endl;
        // Initialize petsc solvers
        PetscInitializeNoArguments();
    }

    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory( mMapType );

    mMap = tMatFactory.create_map( aSolverInterface->get_my_local_global_overlapping_map());

    uint tNumRHMS = aSolverInterface->get_num_rhs();

    // full vector
    mFullVector = tMatFactory.create_vector( aSolverInterface, mMap, tNumRHMS );

//    mDummyFullVector = tMatFactory.create_vector( aSolverInterface, mMap, tNumRHMS );       // FIXME delete
//    mDummyFullVector->vec_put_scalar( 0.0 );
//    aSolverInterface->set_solution_vector_prev_time_step(mDummyFullVector);

    mFullVector->vec_put_scalar( 0.0 );

    // delete pointers if they already exist
    this->delete_pointers();

    // create solver object
    if ( mBuildLinerSystemFlag )
    {
        MORIS_LOG_INFO( "Build linear problem with index %-5i", mNonlinearSolverManagerIndex );

        // create solver factory
        Solver_Factory  tSolFactory;

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
//        delete( mDummyFullVector );
    }

    if(mIsMasterSystem)
    {
        if ( mMapType == sol::MapType::Petsc)
        {
            std::cout<<"selete linear system"<<std::endl;
            PetscFinalize();
            std::cout<<"selete linear system"<<std::endl;
        }
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
    Tracer tTracer(EntityBase::NonLinearProblem, EntityType::NoType, EntityAction::Build);

    // Set VectorFreeSol and LHS
    mLinearProblem->set_free_solver_LHS( mFullVector );

//    this->print_sol_vec( aNonLinearIt );

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
    Tracer tTracer(EntityBase::NonLinearProblem, EntityType::NoType, EntityAction::Build);

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

void Nonlinear_Problem::extract_my_values( const moris::uint                            & aNumIndices,
                                           const moris::Matrix< DDSMat >                & aGlobalBlockRows,
                                           const moris::uint                            & aBlockRowOffsets,
                                                 moris::Cell< moris::Matrix< DDRMat > > & LHSValues )
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
