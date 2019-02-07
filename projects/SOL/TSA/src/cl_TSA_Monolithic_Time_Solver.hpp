/*
 * cl_TSA_Monolithic_Time_Solver.hpp
 *
 *  Created on: Feb 02, 2019
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_
#define MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_

#include <iostream>

#include "cl_TSA_Time_Solver.hpp"
#include "cl_Vector.hpp"

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Matrix_Vector_Factory.hpp"


// MORIS header files.

namespace moris
{
class Map_Class;
class Dist_Vector;
class Solver_Interface;

namespace tsa
{
    class Monolithic_Time_Solver : public Time_Solver
    {
    private:

    protected:
        //! Pointer to my nonlinear solver manager


    public:
        //-------------------------------------------------------------------------------

        Monolithic_Time_Solver(){};

        //-------------------------------------------------------------------------------

        ~Monolithic_Time_Solver(){};

        //-------------------------------------------------------------------------------

        void finalize()
        {
            mSolverInterface = mDatabase->get_solver_interface();

            mFullVector = mDatabase->get_full_vector();

            Matrix_Vector_Factory tMatFactory( MapType::Epetra );

            mPrevFullVector = tMatFactory.create_vector( mSolverInterface, mDatabase->get_list_of_maps( 1 ), VectorType::FREE );          //FIXME
            mPrevFullVector->vec_put_scalar( 0.0 );
        }

        //-------------------------------------------------------------------------------

        void solve()
        {
            uint tTimeSteps = 1000;
            moris::real tTime = 0;

            for ( uint Ik = 0; Ik < tTimeSteps; Ik++ )
            {
                tTime = tTime + 0.01;
                mSolverInterface->set_time( tTime );

                mSolverInterface->set_solution_vector_prev_time_step( mPrevFullVector );

                mDatabase->solve();

                mPrevFullVector->vec_plus_vec( 1.0, *mFullVector, 0.0);

                this->perform_mapping();
            }
        };

        //-------------------------------------------------------------------------------

        void perform_mapping()
        {
            //mFullVector->vec_put_scalar( 0.0 );

            Matrix< DDSMat > tGlobalRows( 1, 1, 0 );
            Matrix< DDRMat > tMat;
            mPrevFullVector->extract_my_values( 1, tGlobalRows, 0 , tMat );

            Matrix< DDSMat > tGlobalRows1( 1, 1, 1 );
            mPrevFullVector->sum_into_global_values( 1, tGlobalRows1, tMat );

            mPrevFullVector->vector_global_asembly();
        };

    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_ */
