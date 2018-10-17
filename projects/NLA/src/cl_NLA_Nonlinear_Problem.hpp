#ifndef MORIS_DISTLINALG_CL_NLA_NONLINEAR_PROBLEM_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLINEAR_PROBLEM_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"
#include "cl_DLA_Linear_Problem.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
class Map_Class;
class Dist_Vector;
class Solver_Interface;
namespace dla
{
    class Linear_Problem;
}
namespace NLA
{
    class Nonlinear_Problem
    {
    private:

    protected:
        Dist_Vector * mVectorFullSol;
        Dist_Vector * mPrevVectorFullSol;

        Map_Class   * mMap;

        std::shared_ptr< dla::Linear_Problem > mLinearProblem;

    public:
        Nonlinear_Problem( Solver_Interface * aSolverInterface );

        ~Nonlinear_Problem();

        void build_linearized_problem();

        std::shared_ptr< dla::Linear_Problem > get_linearized_problem(){ return mLinearProblem; };

        Dist_Vector * get_full_vector();

        void extract_my_values( const moris::uint             & aNumIndices,
                                const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                const moris::uint             & aBlockRowOffsets,
                                      moris::Matrix< DDRMat > & LHSValues );

    };
}
}
#endif /* MORIS_DISTLINALG_CL_NLA_NONLINEAR_PROBLEM_HPP_ */
