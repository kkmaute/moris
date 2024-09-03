/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Eigen_Solver_SLEPc.hpp
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
#include <slepceps.h>

#include "fn_PRM_SOL_Parameters.hpp"

namespace moris::dla
{
    class Eigen_Solver_SLEPc : public moris::dla::Linear_Solver_Algorithm_Petsc

    {
      private:
        EPS mEps;    // mEps is the eigenvalue solver object
        PC  mPc;
        ST  mSt;
        KSP mKsp;

        std::map< std::string, EPSWhich >           mStringToEPSWhich;
        std::map< std::string, STType >             mMorisSTToSTType;
        std::map< std::string, EPSPowerShiftType >  mMorisToPowerType;
        std::map< EPSConvergedReason, std::string > mConvergenceReasonToString;

        // linear solver option in eigen solver
        const Parameter_List* mSubSolverParameterlist;
        const Parameter_List* mSubSolverPreconditionerParameterlist;

        EPSConvergedReason reason;

        Vector<real> mEigenValues;

      protected:

      public:
        //------------------------------------------------------------------------------

        // Eigen_Solver_SLEPc();

        //------------------------------------------------------------------------------

        Eigen_Solver_SLEPc( const moris::Parameter_List& aParameterlist = prm::create_slepc_algorithm_parameter_list() );

        //------------------------------------------------------------------------------

        Eigen_Solver_SLEPc( moris::Solver_Interface* aInput );

        //------------------------------------------------------------------------------

        Eigen_Solver_SLEPc( Linear_Problem* aLinearSystem );

        //------------------------------------------------------------------------------

        ~Eigen_Solver_SLEPc() override;

        //------------------------------------------------------------------------------

        void set_solver_parameters();

        //------------------------------------------------------------------------------

        void construct_solver_and_preconditioner( Linear_Problem* aLinearSystem );

        //------------------------------------------------------------------------------

        moris::sint solve_linear_system() override;

        //------------------------------------------------------------------------------

        moris::sint solve_linear_system(
                Linear_Problem*   aLinearSystem,
                const moris::sint aIter ) override;

        //------------------------------------------------------------------------------

        EPSProblemType determine_problem_type( Linear_Problem* aLinearSystem );

        //------------------------------------------------------------------------------

        void set_sublinear_solver_options( const Parameter_List* aParameterlistsubSolver, const Parameter_List* aParameterlistPreconditioner ) override;

        //------------------------------------------------------------------------------

        void set_sublinear_solver_and_preconditioner( Linear_Problem* aLinearSystem );

        //------------------------------------------------------------------------------

        void build_preconditioner( Linear_Problem* aLinearSystem );

        //------------------------------------------------------------------------------

        void
        set_eps_type_and_params();

        //------------------------------------------------------------------------------

        void
        print_slepc_determined_solver_paramaters();

        //------------------------------------------------------------------------------

        Vector<real> const &
        get_eigenvalues() const;
    };
}    // namespace moris::dla
