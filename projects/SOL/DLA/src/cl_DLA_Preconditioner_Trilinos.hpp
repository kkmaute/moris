/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_Trilinos.hpp
 *
 */

#pragma once

// TPL header files
#include "Epetra_ConfigDefs.h"

#include "cl_DLA_Linear_Solver_Algorithm.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"

#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_preconditioner.h"

#include "cl_DLA_Preconditioner.hpp"

namespace moris::dla
{
    class Preconditioner_Trilinos : public Preconditioner
    {
      private:
        // possible trillions preconditioner
        Teuchos::RCP< Ifpack_Preconditioner >               mIfPackPrec;
        Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > mMlPrec;

        //-------------------------------------------------------------------------------

        moris::sint build_ifpack_preconditioner();

        //-------------------------------------------------------------------------------

        moris::sint build_ml_preconditioner();

        //-------------------------------------------------------------------------------

        moris::sint compute_ifpack_preconditioner( bool tRecompute = false );

        //-------------------------------------------------------------------------------

        moris::sint compute_ml_preconditioner( bool tRecompute = false );

        //-------------------------------------------------------------------------------

      public:

        //-------------------------------------------------------------------------------

        explicit Preconditioner_Trilinos(
                const Parameter_List& aParameterList );

        //-------------------------------------------------------------------------------

        /*
         * initialize preconditioner by setting parameter list and linear system
         */
        void initialize();

        //-------------------------------------------------------------------------------

        /*
         * build and compute preconditioner
         *
         *  @param[in] iteration index - used to decided whether preconditioner needs to
         *                               be build and computed or just recomputed
         */
        void build( Linear_Problem* aProblem, const sint& aIter = 1 ) override;

        //-------------------------------------------------------------------------------

        /*
         * returns true if a preconditioner has been built
         */
        bool exists() override;

        //-------------------------------------------------------------------------------

        /*
         * accessor to underling Epetra operator
         */
        Teuchos::RCP< Epetra_Operator > get_operator();

        //-------------------------------------------------------------------------------

        /*
         * returns Ifpack Preconditioner
         */
        Teuchos::RCP< Ifpack_Preconditioner >&
        get_ifpack_prec()
        {
            return mIfPackPrec;
        };

        //-------------------------------------------------------------------------------

        /*
         * returns Multilevel preconditioner
         */
        Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner >&
        get_ml_prec()
        {
            return mMlPrec;
        };

        //-------------------------------------------------------------------------------
    };
}    // namespace moris
