/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Solver_Factory.hpp
 *
 */

#pragma once

#include <memory>

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

namespace moris
{
    class Solver_Interface;
    namespace NLA
    {
        class Nonlinear_Algorithm;
        class Nonlinear_Solver_Factory
        {
        private:

        protected:

        public:
            Nonlinear_Solver_Factory() = default;

            ~Nonlinear_Solver_Factory() = default;

            std::shared_ptr< Nonlinear_Algorithm > create_nonlinear_solver( const Parameter_List& aParameterList = prm::create_nonlinear_algorithm_parameter_list() );
        };
    }
}
