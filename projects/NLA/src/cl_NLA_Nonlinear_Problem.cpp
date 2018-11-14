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

Nonlinear_Problem::Nonlinear_Problem( Solver_Interface * aSolverInterface )
{
    // create solver factory
   this->set_interface( aSolverInterface );
}

void Nonlinear_Problem::set_interface( Solver_Interface * aSolverInterface )
{
    // delete pointers if they already exist
    this->delete_pointers();

    // create solver factory
    Solver_Factory  tSolFactory;

    // create solver object
    mLinearProblem = tSolFactory.create_linear_system( aSolverInterface, MapType::Epetra );

    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory;

    // create map object
    mMap = tMatFactory.create_map( aSolverInterface->get_num_my_dofs(),
                                   aSolverInterface->get_my_local_global_map(),
                                   aSolverInterface->get_constr_dof(),
                                   aSolverInterface->get_my_local_global_overlapping_map());

    // Build free and full vector
    mVectorFullSol = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );
    mPrevVectorFullSol = tMatFactory.create_vector( aSolverInterface, mMap, VectorType::FULL_OVERLAPPING );

    mVectorFullSol->vec_put_scalar( 0.0 );

    // set flag that interface has been set
    mHasSolverInterface = true;

}

Nonlinear_Problem::~Nonlinear_Problem()
{
    this->delete_pointers();
}

void Nonlinear_Problem::delete_pointers()
{
    // test if interface has been set
    if( mHasSolverInterface )
    {
        delete( mVectorFullSol );
        delete( mPrevVectorFullSol );
        delete( mMap );
        mHasSolverInterface = false;
    };
}

void Nonlinear_Problem::build_linearized_problem( const bool & aRebuildJacobian, const sint aNonLinearIt )
{
    // Set VectorFreeSol and LHS
    mLinearProblem->set_free_solver_LHS( mVectorFullSol );

    this->print_sol_vec( aNonLinearIt );


    if( aRebuildJacobian )
    {
        mLinearProblem->assemble_jacobian( mVectorFullSol );
    }

    mLinearProblem->assemble_residual( mVectorFullSol );
}


void Nonlinear_Problem::build_linearized_problem( const bool & aRebuildJacobian, const sint aNonLinearIt, const sint aRestart )
{
    delete( mVectorFullSol );


    // Build Matrix vector factory
    Matrix_Vector_Factory tMatFactory;
    mVectorFullSol = tMatFactory.create_vector();

    this->restart_from_sol_vec( aRestart );

    // Set VectorFreeSol and LHS
    mLinearProblem->set_free_solver_LHS( mVectorFullSol );

    if( aRebuildJacobian )
    {
        mLinearProblem->assemble_jacobian( mVectorFullSol );
    }

    mLinearProblem->assemble_residual( mVectorFullSol );
}

Dist_Vector * Nonlinear_Problem::get_full_vector()
{
    return mVectorFullSol;
}

void Nonlinear_Problem::extract_my_values( const moris::uint         & aNumIndices,
                                       const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                       const moris::uint             & aBlockRowOffsets,
                                             moris::Matrix< DDRMat > & LHSValues )
{
    mVectorFullSol->save_vector_to_HDF5( "aaa" );

    mVectorFullSol->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
}

void Nonlinear_Problem::print_sol_vec( const sint aNonLinearIt )
{
    char NonLinNum[10];
    std::sprintf( NonLinNum, "NonLIt.%04u", aNonLinearIt );

    char SolVector[100];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector,".h5\0");

    mVectorFullSol->save_vector_to_HDF5( SolVector );
}

void Nonlinear_Problem::restart_from_sol_vec( const sint aRestart )
{
    char NonLinNum[10];
    std::sprintf( NonLinNum, "NonLIt.%04u", aRestart );

    char SolVector[100];
    std::strcpy( SolVector, "SolVector." );
    std::strcat( SolVector, NonLinNum );
    std::strcat( SolVector,".h5\0");

    mVectorFullSol->read_vector_from_HDF5( SolVector );
}

