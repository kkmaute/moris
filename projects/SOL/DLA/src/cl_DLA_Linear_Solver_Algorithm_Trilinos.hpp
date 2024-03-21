/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp
 *
 */

#pragma once

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_Parameter_List.hpp"    // CON/src
#include "cl_DLA_Preconditioner_Trilinos.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"

namespace moris::dla
{
    class Linear_Problem;
    class Linear_Solver_Algorithm_Trilinos : public Linear_Solver_Algorithm
    {
      private:

      protected:
        Preconditioner_Trilinos* mPreconditioner = nullptr;

      public:

        //-----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Linear_Solver_Algorithm_Trilinos object
         *
         * @param aParameterlist
         */

        Linear_Solver_Algorithm_Trilinos( const moris::Parameter_List& aParameterlist )
                : Linear_Solver_Algorithm( aParameterlist ){};

        //-----------------------------------------------------------------------------------

        /**
         * @brief Destroy the Linear_Solver_Algorithm_Trilinos object
         *
         */

        virtual ~Linear_Solver_Algorithm_Trilinos(){};

        //-----------------------------------------------------------------------------------

        /**
         * @brief solve linear system
         *
         * @return moris::sint
         */

        virtual moris::sint solve_linear_system() = 0;

        //-----------------------------------------------------------------------------------

        /**
         * @brief solves linear system x=k^-1*b
         *
         * @param aLinearSystem
         * @param aIter
         * @return moris::sint
         */

        virtual moris::sint solve_linear_system( Linear_Problem* aLinearSystem,
                const moris::sint                                aIter = 1 ) = 0;

        //-----------------------------------------------------------------------------------

        /**
         * @brief Set the preconditioner for Trilinos
         *
         * @param aPreconditioner
         */

        // cast the preconditioner to the correct type(trillons and assign it to the object)
        virtual void
        set_preconditioner( Preconditioner* aPreconditioner ) override;

        //-----------------------------------------------------------------------------------

        /**
         * @brief compute the condition number of the operator with arma/eigen
         *
         */

        virtual void compute_operator_condition_number_with_moris( std::string tComputationMode ) override;

        //-----------------------------------------------------------------------------------

        /**
         * @brief compute the condition number of the preconditioned operator with arma/eigen
         *
         */

        virtual void compute_preconditioned_operator_condition_number_with_moris( std::string tComputationMode ) override;
    };
}    // namespace moris::dla
