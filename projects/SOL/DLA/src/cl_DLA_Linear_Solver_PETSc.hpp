/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_PETSc.hpp
 *
 */

#pragma once

#include "core.hpp"
#include "cl_DLA_Linear_Solver_Algorithm_Petsc.hpp"
#include "cl_Vector_PETSc.hpp"
#include "cl_MatrixPETSc.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_DLA_Linear_Problem.hpp"

namespace moris::dla
{
    class Linear_Solver_PETSc : public moris::dla::Linear_Solver_Algorithm_Petsc

    {
      private:
        KSP mPetscKSPProblem;

        PC mpc;

        moris::Cell< KSP > tKSPBlock;

        friend class Preconditioner_PETSc;

      protected:

      public:
        //------------------------------------------------------------------------------

        Linear_Solver_PETSc();

        //------------------------------------------------------------------------------

        Linear_Solver_PETSc( const moris::ParameterList aParameterlist );

        //------------------------------------------------------------------------------

        Linear_Solver_PETSc( moris::Solver_Interface* aInput );

        //------------------------------------------------------------------------------

        Linear_Solver_PETSc( Linear_Problem* aLinearSystem );

        //------------------------------------------------------------------------------

        ~Linear_Solver_PETSc();

        //------------------------------------------------------------------------------

        void set_solver_parameters();

        //------------------------------------------------------------------------------

        void construct_solver_and_preconditioner( Linear_Problem* aLinearSystem );

        //------------------------------------------------------------------------------

        moris::sint solve_linear_system();

        //------------------------------------------------------------------------------

        moris::sint solve_linear_system(
                Linear_Problem*   aLinearSystem,
                const moris::sint aIter );

        //------------------------------------------------------------------------------

        void set_solver_analysis_options();

        //------------------------------------------------------------------------------

        /**
         * @brief compute the n first eigen values and eigen vectors
         *        the eigenvector needs to be defined through an IQI
         *
         * @param aLinearSystem
         */
        void compute_eigenspectrum( Linear_Problem* aLinearSystem );
    };
}    // namespace moris::dla


