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
    MORIS_ASSERT( !mHasSolverInterface, "Interface has already been set.");
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

    mHasSolverInterface = true;

}

Nonlinear_Problem::~Nonlinear_Problem()
{
	if( mHasSolverInterface )
	{
		delete( mVectorFullSol );
		delete( mPrevVectorFullSol );
		delete( mMap );
	}
}

void Nonlinear_Problem::build_linearized_problem()
{
    // Set VectorFreeSol and LHS
    mLinearProblem->set_free_solver_LHS( mVectorFullSol );

    mLinearProblem->assemble_residual_and_jacobian( mVectorFullSol );
    //mLinearProblem->assemble_residual( mVectorFullSol );
    //mLinearProblem->assemble_jacobian( mVectorFullSol );
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
    mVectorFullSol->extract_my_values( aNumIndices, aGlobalBlockRows, aBlockRowOffsets, LHSValues );
}
