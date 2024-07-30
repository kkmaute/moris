/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner.hpp
 *
 */

#pragma once

// TPL header files
#include "cl_DLA_Linear_Solver_Algorithm.hpp"

namespace moris::dla
{

    class Preconditioner
    {
      protected:
        bool mIsInitialized = false;
        Linear_Problem* mLinearSystem = nullptr;
        const Parameter_List& mParameterList;

      public:

        /**
         * Constructor
         */
        explicit Preconditioner( const Parameter_List& aParameterList )
                : mParameterList( aParameterList )
        {
        }

        /**
         * Destructor
         */
        virtual ~Preconditioner() = default;

        /*
         * build and compute preconditioner
         *
         *  @param[in] iteration index - used to decided whether preconditioner needs to
         *                               be build and computed or just recomputed
         */
        virtual void build( Linear_Problem* aProblem, const sint& aIter = 1 ) {};

        //-------------------------------------------------------------------------------

        /**
         * returns true if a preconditioner has been built
         **/

        virtual bool exists() { return mIsInitialized; };

        //-------------------------------------------------------------------------------
    };
}    // namespace dla
