/*
 * cl_DLA_Linear_Solver_Algorithm.hpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_DLA_LINEAR_SOLVER_ALGORITHM_HPP_
#define MORIS_DISTLINALG_CL_DLA_LINEAR_SOLVER_ALGORITHM_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Enums.hpp"

#include "cl_Param_List.hpp" // CON/src

namespace moris
{
namespace dla
{
    class Linear_Problem;
    class Linear_Solver_Algorithm
    {
    private:

    protected:

        moris::real mCondEstimate;

        moris::uint mSolNumIters;
        moris::real mSolTrueResidual;
        moris::real mSolScaledResidual;
        moris::real mSolTime;
        moris::real mSymFactTime;
        moris::real mNumFactTime;
        moris::real mPreCondTime;

        Param_List< boost::variant< bool, sint, real, std::string > > mParameterList; // The Algorithm specific parameter list

    public:
        Linear_Solver_Algorithm( ){};

        virtual ~Linear_Solver_Algorithm(){};

        virtual void set_linear_problem( Linear_Problem * aLinearSystem ) = 0;

        virtual moris::sint solve_linear_system() = 0;

        virtual moris::sint solve_linear_system( Linear_Problem * aLinearSystem, const moris::sint aIter = 1 ) = 0;

//        Dist_Vector * get_solver_LHS()
//        {
//            return mFreeVectorLHS;
//        };
//
//        Dist_Vector * get_solver_RHS()
//        {
//            return mVectorRHS;
//        };

//        auto get_solver_input() const ->decltype( mInput )
//        {
//            return mInput;
//        };

//        virtual void get_solution( moris::Matrix< DDRMat > & LHSValues ) =0;

        boost::variant< bool, sint, real, std::string > &  set_param( const std::string & aKey )
        {
            return mParameterList( aKey );
        }
    };
}
}
#endif /* MORIS_DISTLINALG_CL_DLA_LINEAR_SOLVER_ALGORITHM_HPP_ */

