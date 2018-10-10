#ifndef MORIS_DISTLINALG_CL_LinearSolver_HPP_
#define MORIS_DISTLINALG_CL_LinearSolver_HPP_

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Matrix_Vector_Factory.hpp"
#include "cl_DistLinAlg_Enums.hpp"

#include "cl_Param_List.hpp" // CON/src
#include "cl_DistLinAlg_Enums.hpp"

namespace moris
{
class Sparse_Matrix;
class Dist_Vector;
class Map_Class;
class Solver_Interface;
class Linear_Solver
{
private:

protected:
    Sparse_Matrix * mMat;
    Dist_Vector   * mVectorRHS;
    Dist_Vector   * mVectorLHS;
    Dist_Vector   * mVectorLHSOverlapping;
    Map_Class     * mMap;

    Solver_Interface * mInput;

    moris::real mCondEstimate;

    moris::uint mSolNumIters;
    moris::real mSolTrueResidual;
    moris::real mSolScaledResidual;
    moris::real mSolTime;
    moris::real mSymFactTime;
    moris::real mNumFactTime;
    moris::real mPreCondTime;

public:
    Linear_Solver( Solver_Interface *  aInput ) : mMat(NULL),
                                              mVectorRHS(NULL),
                                              mVectorLHS(NULL),
                                              mMap(NULL),
                                              mInput( aInput )
    {};

    virtual ~Linear_Solver(){};

//    virtual void build_linear_system( Epetra_FECrsMatrix*       aEpetraMat,
//                                      Epetra_FEVector*          aEpetraVector_x,
//                                      Epetra_FEVector*          aEpetraVector_b ) = 0;

    virtual void assemble_residual_and_jacobian( Dist_Vector * aFullSolutionVector ) = 0;

    virtual void assemble_residual_and_jacobian(  ) = 0;

    virtual void build_linear_system() = 0;

    virtual moris::sint solve_linear_system() = 0;

    virtual void solve_eigenvalues() = 0;

    Dist_Vector * get_solver_LHS()
    {
        return mVectorLHS;
    };

    Dist_Vector * get_solver_RHS()
    {
        return mVectorRHS;
    };

    auto get_solver_input() const ->decltype( mInput )
    {
        return mInput;
    };

    virtual void get_solution( moris::Matrix< DDRMat > & LHSValues ) =0;

    virtual boost::variant< bool, sint, real, const char* > & set_param( char const* aKey ) = 0;
};
}
#endif /* MORIS_DISTLINALG_CL_LinearSolver_HPP_ */

