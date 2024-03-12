/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_PETSc.hpp
 *
 */

#pragma once

#include <memory>

#include "cl_MatrixPETSc.hpp"


#include "cl_DLA_Preconditioner.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map;
        class Dist_Vector;
        class Dist_Matrix;
    }    // namespace sol

    namespace dla
    {
        class Geometric_Multigrid;
        class Linear_Solver_PETSc;
        class Linear_Problem;
        class Preconditioner_PETSc : public Preconditioner
        {
          private:
            PC                mpc;
            sol::Dist_Matrix* mPreconMat = nullptr;
            sol::Dist_Map*    mMapFree   = nullptr;

            dla::Geometric_Multigrid* mGeoMultigrid = nullptr;

          protected:

          public:
            Preconditioner_PETSc();
            // Preconditioner_PETSc( Linear_Solver_PETSc* aLinearSolverAlgoritm );

            ~Preconditioner_PETSc(){
                //            if( mPreconMat != nullptr )
                //            {
                //                delete mPreconMat;
                //            }
                //            if( mMapFree != nullptr )
                //            {
                //                delete mMapFree;
                //            }
            };

            Preconditioner_PETSc(
                    moris::ParameterList* aParameterlist );

            sol::Dist_Matrix*
            get_preconditioner_matrix()
            {
                return mPreconMat;
            };

            //-----------------------------------------------------------------------------------

            void build_ilu_preconditioner( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------------------

            void build_multigrid_preconditioner( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------------------

            void build_schwarz_preconditioner_petsc( Linear_Problem* aLinearSystem, KSP aPetscKSPProblem );

            //-----------------------------------------------------------------------------------

            void build_schwarz_preconditioner( Linear_Problem* aLinearSystem );

            //-----------------------------------------------------------------------------------

            /**
             * @brief compute the preconditioner based on the type of preconditioner
             *
             * @param aLinearSystem linear problem conraining the matrix
             * @param aPetscKSPProblem KSP problem conext used in petsc
             */

            void build_preconditioner( Linear_Problem* aLinearSystem, KSP aPetscKSPProblem );

            //-----------------------------------------------------------------------------------
        };
    }    // namespace dla
}    // namespace moris
