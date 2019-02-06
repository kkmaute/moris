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

        //! Full Vector
        Dist_Vector * mFullVector = nullptr;

        //! Full Vector
        Dist_Vector * mPrevFullVector = nullptr;


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

        void solve()
        {
            uint tTimeSteps = 100;
            moris::real tTime = 0;
            Matrix< DDRMat > tMat (tTimeSteps, 1 , 0);
            for ( uint Ik = 0; Ik < tTimeSteps; Ik++ )
            {
                tTime = tTime + 0.1;
                mSolverInterface->set_time( tTime );

                mSolverInterface->set_solution_vector_prev_time_step( mPrevFullVector );

                mDatabase->solve();

                mPrevFullVector->vec_plus_vec( 1.0, *mFullVector, 0.0);

                //-------------------------------------------------------------------------------------
                Matrix< DDRMat > mMySolVec;
                mPrevFullVector->extract_copy( mMySolVec );
                tMat(Ik,0) = mMySolVec(0,0);
            }
            print(tMat,"mMySolVec");
        };

    };
}
}
#endif /* MORIS_DISTLINALG_CL_TSA_MONOLITHIC_TIME_SOLVER_HPP_ */
