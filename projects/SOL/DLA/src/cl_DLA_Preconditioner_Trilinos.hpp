/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_Trilinos.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_PRECONDITIONER_TRILINOS_HPP_
#define SRC_DISTLINALG_CL_PRECONDITIONER_TRILINOS_HPP_

// TPL header files
#include "Epetra_ConfigDefs.h"

#include "cl_DLA_Linear_Solver_Algorithm.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"

#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra_preconditioner.h"

namespace moris
{
    namespace dla
    {
        class Preconditioner_Trilinos
        {
            private:

                bool mIsInitialized = false;

                Linear_Problem   * mLinearSystem =  nullptr;

                moris::ParameterList mParameterList;

                Teuchos::RCP< Ifpack_Preconditioner >  mIfPackPrec;

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

                Preconditioner_Trilinos();

                //-------------------------------------------------------------------------------

                Preconditioner_Trilinos(
                        const moris::ParameterList   aParameterlist,
                        Linear_Problem             * aLinearSystem  );

                //-------------------------------------------------------------------------------

                ~Preconditioner_Trilinos();

                //-------------------------------------------------------------------------------

                /*
                 * initialize preconditioner by setting parameter list and linear system
                 */
                void initialize(
                        const moris::ParameterList   aParameterlist,
                        Linear_Problem             * aLinearSystem);

                //-------------------------------------------------------------------------------

                /*
                 * build and compute preconditioner
                 *
                 *  @param[in] iteration index - used to decided whether preconditioner needs to
                 *                               be build and computed or just recomputed
                 */
                void build( const sint & aIter = 1 );

                //-------------------------------------------------------------------------------

                /*
                 * returns true if a preconditioner has been built
                 */
                bool exists();

                //-------------------------------------------------------------------------------

                /*
                 * accessor to underling Eptra operator
                 */
                Teuchos::RCP< Epetra_Operator > get_operator();

                //-------------------------------------------------------------------------------
        };
    }
}

#endif /* SRC_DISTLINALG_CL_PRECONDITIONER_TRILINOS_HPP_ */

