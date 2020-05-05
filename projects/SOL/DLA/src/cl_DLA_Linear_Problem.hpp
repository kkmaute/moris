/*
 * cl_DLA_Linear_System.hpp
 *
 *  Created on: May 16, 2018
 *      Author: schmidt
 */
#ifndef MORIS_DISTLINALG_CL_DLA_LINEAR_SYSTEM_HPP_
#define MORIS_DISTLINALG_CL_DLA_LINEAR_SYSTEM_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
class Dist_Matrix;
class Dist_Vector;
class Dist_Map;
class Solver_Interface;
namespace dla
{
    class Linear_Problem
    {
    private:

    protected:
        Dist_Matrix   * mMat;
        Dist_Vector   * mVectorRHS;
        Dist_Vector   * mFreeVectorLHS;
        Dist_Vector   * mFullVectorLHS;
        Dist_Map      * mMap = nullptr;
        Dist_Map      * mMapFree= nullptr;

        Solver_Interface * mSolverInterface = nullptr;

        moris::real mCondEstimate;

    public:
        Linear_Problem( Solver_Interface * aInput ) : mMat(NULL),
                                                      mVectorRHS(NULL),
                                                      mFreeVectorLHS(nullptr),
                                                      mMap(NULL),
                                                      mSolverInterface( aInput )
        {};

        virtual ~Linear_Problem(){};

        void assemble_residual_and_jacobian(  );
        void assemble_residual_and_jacobian( Dist_Vector * aFullSolutionVector );
        void assemble_residual();
        void assemble_jacobian();

        virtual moris::sint solve_linear_system() = 0;

        Dist_Vector * get_free_solver_LHS() { return mFreeVectorLHS; };

        void set_free_solver_LHS( Dist_Vector * aFullSolVector);

        Dist_Vector * get_full_solver_LHS();

        Dist_Vector * get_solver_RHS() { return mVectorRHS; };

        Dist_Matrix * get_matrix() { return mMat; };

        Solver_Interface * get_solver_input() const { return mSolverInterface; };

        virtual void get_solution( moris::Matrix< DDRMat > & LHSValues ) =0;
    };
}
}
#endif /* MORIS_DISTLINALG_CL_DLA_LINEAR_SYSTEM_HPP_ */

