/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Belos.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_BELOS_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_BELOS_HPP_

// TPL header files
#include "Epetra_ConfigDefs.h"

#include "cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosGCRODRSolMgr.hpp"

#include "Teuchos_ParameterList.hpp"

namespace moris
{
namespace dla
{
class Linear_Solver_Belos : public Linear_Solver_Algorithm_Trilinos
{
private:

    Teuchos::RCP< Teuchos::ParameterList > mMyPl;

protected:
public:
    Linear_Solver_Belos( const moris::Parameter_List& aParameterlist = prm::create_linear_algorithm_parameter_list_belos() );

    Linear_Solver_Belos( Linear_Problem * aLinearSystem  );

    ~Linear_Solver_Belos();

    //int SetSystemMatrix ( bool aUseTranspose );

    moris::sint solve_linear_system(){ return 0; };

    moris::sint solve_linear_system(      Linear_Problem * aLinearSystem,
                                     const moris::sint     aIter );

    void set_solver_internal_parameters();
};
}
}

#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_BELOS_HPP_ */

