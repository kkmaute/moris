/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Algorithm.hpp
 *
 */

#ifndef MORIS_DISTLINALG_CL_DLA_LINEAR_SOLVER_ALGORITHM_HPP_
#define MORIS_DISTLINALG_CL_DLA_LINEAR_SOLVER_ALGORITHM_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_Parameter_List.hpp"    // CON/src

namespace moris
{
    class Solver_Interface;

    namespace dla
    {
        class Linear_Problem;
        class Preconditioner;
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

            Solver_Interface* mSolverInterface = nullptr;

            moris::Parameter_List mParameterList;    // The Algorithm specific parameter list

            Linear_Problem* mLinearSystem = nullptr;

          public:
            //-----------------------------------------------------------------------------------
            Linear_Solver_Algorithm( const moris::Parameter_List& aParameterlist )
                    : mParameterList( aParameterlist ) {};

            //-----------------------------------------------------------------------------------

            virtual ~Linear_Solver_Algorithm() {};

            //-----------------------------------------------------------------------------------

            virtual moris::sint solve_linear_system() = 0;

            //-----------------------------------------------------------------------------------

            virtual moris::sint solve_linear_system( Linear_Problem* aLinearSystem,
                    const moris::sint                                aIter = 1 ) = 0;

            //-----------------------------------------------------------------------------------
            virtual void
            set_preconditioner( Preconditioner* aPreconditioner ) {};

            //-----------------------------------------------------------------------------------

            /**
             * This function is only used for eigensolver where the left hand side matrix (K/A) can be altered.
             *
             * @param aPreconditioner
             */

            virtual void set_left_hand_side_preconditioner( Preconditioner* aPreconditioner ) {};

            //-----------------------------------------------------------------------------------

            /**
             * @brief Set the sublinear solver options object
             *
             * @param aParameterlistsubSolver linear solver needed for inverse
             * @param aParameterlistPreconditioner precondioner for the linear solver
             */

            virtual void set_sublinear_solver_options( const Parameter_List* aParameterlistsubSolver, const Parameter_List* aParameterlistPreconditioner ) {};

            //-----------------------------------------------------------------------------------

            /**
             * @brief compute the condition number of the operator with moris arma/eigen
             *
             */

            virtual void compute_operator_condition_number_with_moris( const std::string& tComputationMode ) {};

            //-----------------------------------------------------------------------------------

            /**
             * @brief compute the condition number of the preconditioned operator with moris arma/eigen
             *
             */

            virtual void compute_preconditioned_operator_condition_number_with_moris( const std::string& tComputationMode ) {};
        };
    }    // namespace dla
}    // namespace moris
#endif /* MORIS_DISTLINALG_CL_DLA_LINEAR_SOLVER_ALGORITHM_HPP_ */
