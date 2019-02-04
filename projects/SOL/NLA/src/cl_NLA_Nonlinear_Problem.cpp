/*
 * cl_NLA_Nonlinear_Problem.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Nonlinear_Problem.hpp"

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

Nonlinear_Problem::Nonlinear_Problem(            Solver_Interface * aSolverInterface,
                                                 Dist_Vector      * aFullVector,
                                      const moris::sint             aNonlinearSolverManagerIndex,
                                      const bool                    aBuildLinerSystemFlag,
                                      const enum MapType            aMapType) :     mFullVector( aFullVector ),
                                                                                    mMapType( aMapType ),
                                                                                    mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
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

    mBuildLinerSystemFlag = aBuildLinerSystemFlag;
    // create solver factory
    this->set_interface( aSolverInterface );
}

Nonlinear_Problem::Nonlinear_Problem(            Solver_Interface * aSolverInterface,
                                      const moris::sint             aNonlinearSolverManagerIndex,
                                      const bool                    aBuildLinerSystemFlag,
                                      const enum MapType            aMapType) :     mMapType( aMapType ),
                                                                                    mNonlinearSolverManagerIndex( aNonlinearSolverManagerIndex )
{
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

    mFullVector->vec_put_scalar( 0.0 );

    mBuildLinerSystemFlag = aBuildLinerSystemFlag;
    // create solver factory
    this->set_interface( aSolverInterface );
}

void Nonlinear_Problem::set_interface( Solver_Interface * aSolverInterface )
{
    // delete pointers if they already exist
    this->delete_pointers();

    // create solver factory
    Solver_Factory  tSolFactory;

    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory( mMapType );

    // create solver object
    if ( mBuildLinerSystemFlag )
    {
        MORIS_LOG_INFO( "Build linear problem with index %-5i \n", mNonlinearSolverManagerIndex );

        mLinearProblem = tSolFactory.create_linear_system( aSolverInterface, mMapType );
    }

    // set flag that interface has been set
    mIsMasterSystem = true;

}

Nonlinear_Problem::~Nonlinear_Problem()
{
    this->delete_pointers();
	
	    if( mMap != nullptr )
    {
        delete( mMap );
    }

    if( mIsMasterSystem )
    {
        delete( mFullVector );
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

void Nonlinear_Problem::build_linearized_problem( const bool & aRebuildJacobian, const sint aNonLinearIt )
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


void Nonlinear_Problem::build_linearized_problem( const bool & aRebuildJacobian, const sint aNonLinearIt, const sint aRestart )
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
    char NonLinNum[10];
    std::sprintf( NonLinNum, "NonLIt.%04u", aNonLinearIt );

    char SolVector[100];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector,".h5\0");

    mFullVector->save_vector_to_HDF5( SolVector );
}

void Nonlinear_Problem::restart_from_sol_vec( const sint aRestart )
{
    char NonLinNum[10];
    std::sprintf( NonLinNum, "NonLIt.%04u", aRestart );

    char SolVector[100];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector,".h5\0");

    mFullVector->read_vector_from_HDF5( SolVector );
}

