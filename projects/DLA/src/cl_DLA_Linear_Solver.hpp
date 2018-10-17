/*
 * cl_DLA_Linear_Problem.hpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_DLA_Linear_Solver_HPP_
#define MORIS_DISTLINALG_CL_DLA_Linear_Solver_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

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
    class Linear_Solver
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

        Param_List< boost::variant< bool, sint, real  > > mParameterList; // The Algorithm specific parameter list

    public:
        Linear_Solver( )
        {};

        virtual ~Linear_Solver(){};

        virtual void set_linear_problem( std::shared_ptr< Linear_Problem > aLinearSystem ) = 0;

        virtual moris::sint solve_linear_system() = 0;

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

        boost::variant< bool, sint, real > &  set_param( char const* aKey)
        {
            return mParameterList(aKey);
        }
    };
}
}
#endif /* MORIS_DISTLINALG_CL_DLA_Linear_Solver_HPP_ */

