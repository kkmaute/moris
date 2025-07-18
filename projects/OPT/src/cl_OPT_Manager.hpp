/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Manager.hpp
 *
 */

#pragma once

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

// forward declarations
namespace moris
{
    class Module_Parameter_Lists;
    class Submodule_Parameter_Lists;

    namespace opt
    {
      class Algorithm;
      class Problem;
      class Criteria_Interface;
    }
}

namespace moris::opt
{
    class Manager
    {
      private:
        Vector< std::shared_ptr< Algorithm > > mAlgorithms;    // Cell of pointers of optimization algorithms
        std::shared_ptr< Problem >             mProblem;

      public:
        /**
         * Constructor with both parameter list and a cell of criteria interfaces
         *
         * @param aParameterLists parameter lists for defining an optimization problem
         * @param aCriteriaInterface criteria interfaces, in addition to any specified in the parameter lists
         */
        Manager( const Module_Parameter_Lists&             aParameterLists,
                const Vector< std::shared_ptr< Criteria_Interface > >& aInterfaces = Vector< std::shared_ptr< Criteria_Interface > >( 0 ) );

        /**
         * Constructor with problem class already created
         *
         * @param aAlgorithmParameterLists parameter lists for defining just the optimization algorithms
         */
        Manager( const Submodule_Parameter_Lists& aAlgorithmParameterLists,
                std::shared_ptr< Problem >       aProblem );

        /**
         * Destructor
         */
        ~Manager();

        /**
         * Solve the Cell of optimization algorithms
         */
        void perform();

        void reinitialize(
                uint aI,
                uint aOptIteration );

        /**
         * Gets the avs from the problem
         *
         * @return Vector of current advs
         */
        Vector< real > get_advs();

        /**
         * Gets the objectives from the problem
         *
         * @return Vector of objectives, aka a single objective value
         */
        Matrix< DDRMat > get_objectives();
    };
}    // namespace moris::opt
