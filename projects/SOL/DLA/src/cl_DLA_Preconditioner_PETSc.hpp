/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_PETSc.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_PRECONDITIONER_PETSC_HPP_
#define SRC_DISTLINALG_CL_PRECONDITIONER_PETSC_HPP_

#include <memory>
//#include "cl_DLA_Linear_Solver_Algorithm.hpp"
//#include "cl_VectorPETSc.hpp"
//#include "cl_MatrixPETSc.hpp"
//
//#include "cl_SOL_Matrix_Vector_Factory.hpp"
//#include "cl_DLA_Solver_Interface.hpp"
//
//#include "cl_DLA_Linear_Problem.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map;
        class Dist_Vector;
        class Dist_Matrix;
    }

namespace dla
{
class Linear_Solver_PETSc;
class Linear_Problem;
class Preconditioner_PETSc
{
    private:

        Linear_Solver_PETSc * mLinearSolverAlgoritm;

        sol::Dist_Matrix   * mPreconMat = nullptr;
        sol::Dist_Map* mMapFree = nullptr;

    protected:

    public:
        Preconditioner_PETSc( Linear_Solver_PETSc * aLinearSolverAlgoritm );

        ~Preconditioner_PETSc()
        {
//            if( mPreconMat != nullptr )
//            {
//                delete mPreconMat;
//            }
//            if( mMapFree != nullptr )
//            {
//                delete mMapFree;
//            }
        };

        sol::Dist_Matrix * get_preconditioner_matrix(){ return mPreconMat; };

        void build_multigrid_preconditioner( Linear_Problem * aLinearSystem );

        void build_schwarz_preconditioner_petsc();

        void build_schwarz_preconditioner( Linear_Problem * aLinearSystem );

};
}
}

#endif /* SRC_DISTLINALG_CL_PRECONDITIONER_PETSC_HPP_ */

