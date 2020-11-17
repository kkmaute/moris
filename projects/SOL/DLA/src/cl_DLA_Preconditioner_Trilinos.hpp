/*
 * cl_DLA_Preconditioner_Trilinos.hpp
 *
 *  Created on: Jan 06, 2020
 *      Author: schmidt
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

                Linear_Problem   * mLinearSystem =  nullptr;

                Teuchos::RCP< Teuchos::ParameterList > mMyPl;

                Teuchos::RCP< Ifpack_Preconditioner > mPreconditioner;

                moris::ParameterList mParameterList;

                Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > mMlPrec;
                Teuchos::ParameterList                mlParams;

            protected:
            public:

                Preconditioner_Trilinos( const moris::ParameterList aParameterlist,
                        Linear_Problem * aLinearSystem  );

                ~Preconditioner_Trilinos();

                moris::sint build_ifpack_preconditioner();

                moris::sint build_ml_preconditioner();

                Teuchos::RCP< Ifpack_Preconditioner > get_ifpack_prec()
                {
                    return mPreconditioner;
                };

                Teuchos::RCP< ML_Epetra::MultiLevelPreconditioner > get_ml_prec()
                {
                    return mMlPrec;
                };
        };
    }
}

#endif /* SRC_DISTLINALG_CL_PRECONDITIONER_TRILINOS_HPP_ */
