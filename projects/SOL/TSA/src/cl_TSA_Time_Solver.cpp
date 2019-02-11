/*
 * cl_TSA_Time_Solver.cpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#include <ctime>

#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

using namespace moris;
using namespace tsa;

//-------------------------------------------------------------------------------

Time_Solver::Time_Solver( const enum MapType aMapType )
{

}

//-------------------------------------------------------------------------------

Time_Solver::~Time_Solver()
{
    if ( mIsMasterTimeSolver )
    {
        delete( mFullVector );
        delete( mPrevFullVector );
    }
}
//-------------------------------------------------------------------------------

void Time_Solver::get_full_solution( moris::Matrix< DDRMat > & LHSValues )
{
    mFullVector->extract_copy( LHSValues );
}

//-------------------------------------------------------------------------------

moris::real Time_Solver::calculate_time_needed( const clock_t aTime )
{
    moris::real tDeltaTime = (moris::real) ( clock() - aTime ) / CLOCKS_PER_SEC;

    moris::real tDeltaTimeMax   = tDeltaTime;

    MPI_Allreduce( &tDeltaTime, &tDeltaTimeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

    return tDeltaTimeMax;
}

//-------------------------------------------------------------------------------

void Time_Solver::set_dof_type_list( const moris::Cell< enum MSI::Dof_Type > aStaggeredDofTypeList,
                                     const moris::sint                       aLevel )
{
    mMyDofTypeList.push_back( aStaggeredDofTypeList );
}

//-------------------------------------------------------------------------------

void Time_Solver::finalize()
{
    mSolverInterface = mDatabase->get_solver_interface();

    // create map object FIXME ask liner problem for map
    Matrix_Vector_Factory tMatFactory( MapType::Epetra );
    mFullMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map());

    mPrevFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, VectorType::FREE );

    mPrevFullVector->vec_put_scalar( 0.0 );

    if ( mIsMasterTimeSolver )
    {
        // create map object FIXME ask liner problem for map
        //Matrix_Vector_Factory tMatFactory( MapType::Epetra );
        //mFullMap = tMatFactory.create_map( mSolverInterface->get_my_local_global_overlapping_map());

        // full vector and prev full vector
        mFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, VectorType::FREE );
        //mPrevFullVector = tMatFactory.create_vector( mSolverInterface, mFullMap, VectorType::FREE );

        mFullVector->vec_put_scalar( 0.0 );
        //mPrevFullVector->vec_put_scalar( 0.0 );
    }
}
