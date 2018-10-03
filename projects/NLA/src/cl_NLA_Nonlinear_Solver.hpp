#ifndef MORIS_DISTLINALG_CL_NLA_NONLinearSolver_HPP_
#define MORIS_DISTLINALG_CL_NLA_NONLinearSolver_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
class Map_Class;
class Dist_Vector;
class Linear_Solver;
namespace NLA
{
    class Nonlinear_Solver
    {
    private:

    protected:
        Dist_Vector * mVectorFullSol;
        Dist_Vector * mVectorFreeSol;
        Dist_Vector * mPrevVectorFreeSol;
        Dist_Vector * mPrevVectorFullSol;

        Map_Class   * mMap;

        std::shared_ptr< Linear_Solver > mLinearSolver;

        Param_List< boost::variant< bool, sint, real, const char* > > mParameterListNonlinearSolver;

    public:
        virtual ~Nonlinear_Solver(){};

        virtual void set_linear_solver( std::shared_ptr< Linear_Solver > aLinearSolver ) = 0;

        virtual void solver_nonlinear_system() = 0;

        virtual Dist_Vector * get_full_sol_vec() = 0;

        virtual void get_full_solution( moris::Matrix< DDRMat > & LHSValues ) = 0;

        virtual void get_solution( moris::Matrix< DDRMat > & LHSValues ) =0;

        virtual void extract_my_values( const moris::uint             & aNumIndices,
                                        const moris::Matrix< DDSMat > & aGlobalBlockRows,
                                        const moris::uint             & aBlockRowOffsets,
                                              moris::Matrix< DDRMat > & LHSValues ) = 0;

        virtual boost::variant< bool, sint, real, const char* > & set_param( char const* aKey ) = 0;
    };
}
}
#endif /* MORIS_DISTLINALG_CL_NLA_NONLinearSolver_HPP_ */
