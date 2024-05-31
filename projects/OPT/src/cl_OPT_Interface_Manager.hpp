/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Interface_Manager.hpp
 *
 */

#pragma once

#include "cl_OPT_Criteria_Interface.hpp"
#include "cl_Parameter_List.hpp"

namespace moris::opt
{
    class Interface_Manager : public Criteria_Interface
    {
      private:
        Vector< std::shared_ptr< Criteria_Interface > > mInterfaces;

        Matrix< DDUMat > mNumADVsPerInterface;
        Matrix< DDUMat > mNumCriteriaPerInterface;
        Matrix< DDSMat > mProcessorBoundaries;

        uint mNumInterfaces;

        bool mSharedADVs;
        bool mParallel;

      public:
        /**
         * Constructor
         */
        Interface_Manager(
                Parameter_List                                  aParameterList,
                Vector< std::shared_ptr< Criteria_Interface > > aInterfaces );

        /**
         * Destructor
         */
        ~Interface_Manager()
        {
        }

        /**
         * Sets the individual interfaces based on a cell of parameter lists
         */
        void set_interfaces();

        /**
         * Initializes the vectors of ADV values, lower bounds, and upper bounds
         *
         * @param aGlobalADVs Global ADVs
         * @param aGlobalLowerBounds Global lower bounds
         * @param aGlobalUpperBounds Global upper bounds
         */
        void initialize(
                Vector< real >& aGlobalADVs,
                Vector< real >& aGlobalLowerBounds,
                Vector< real >& aGlobalUpperBounds,
                Matrix< IdMat >&  aIjklIds );

        /**
         * Gets the criteria values
         *
         * @return vector of criteria
         */
        Vector< real > perform( Vector< real >& aNewADVs );

        /**
         * Gets the derivative of the criteria with respect to the advs
         *
         * @return matrix d(criteria)_i/d(adv)_j
         */
        Matrix< DDRMat > compute_dcriteria_dadv();

        /**
         * Gets the local advs based on whether or not they are shared
         */
        Vector< real > get_local_advs(
                Vector< real > aGlobalADVs,
                uint           tInterfaceIndex );
    };
}    // namespace moris::opt
