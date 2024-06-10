/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Criteria_Interface.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "cl_Parameter_List.hpp"

namespace moris::opt
{
    class Criteria_Interface
    {
      private:
        bool             mEvaluated = false;
        Matrix< DDRMat > mSensitivities;

      protected:
        bool mInitializeOptimizationRestart = false;

        uint mReinitializeIterInterval = MORIS_UINT_MAX;
        uint mReinitializeIterFirst    = MORIS_UINT_MAX;

        uint mIter = 0;

      public:
        /**
         * Constructor
         */
        Criteria_Interface()
        {
        }

        /**
         * Destructor
         */
        virtual ~Criteria_Interface()
        {
        }

        /**
         * Initializes the vectors of ADV values, lower bounds, and upper bounds
         *
         * @param aADVs Initial ADVs to be filled.
         * @param aLowerBounds Lower ADV bounds to be filled.
         * @param aUpperBounds Upper ADV bounds to be filled.
         */
        virtual void initialize(
                Vector< real >& aADVs,
                Vector< real >& aLowerBounds,
                Vector< real >& aUpperBounds,
                Matrix< IdMat >&  aIjklIds ) = 0;

        /**
         * Gets the criteria values given a new set of ADVs
         *
         * @return vector of criteria
         */
        Vector< real > get_criteria( Vector< real >& aNewADVs );

        /**
         * Gets the derivatives of the criteria with respect to the advs, and computes if not already done
         *
         * @return matrix d(criteria)_i/d(adv)_j
         */
        Matrix< DDRMat > get_dcriteria_dadv();

        /**
         * indicates a request to restart the optimization
         *
         */
        bool
        get_restart_optimization()
        {
            return mInitializeOptimizationRestart;
        }

        void
        set_reinitialize_iter( uint aReinitializeIter )
        {
            mReinitializeIterInterval = aReinitializeIter;
        }

        void
        set_first_reinitialize_iter( uint aReinitializeIter )
        {
            mReinitializeIterFirst = aReinitializeIter;
        }

      private:
        /**
         * Gets the criteria values given a new set of ADVs
         *
         * @return vector of criteria
         */
        virtual Vector< real > perform( Vector< real >& aNewADVs ) = 0;

        /**
         * Computes the derivatives of the criteria with respect to the advs
         *
         * @return matrix d(criteria)_i/d(adv)_j
         */
        virtual Matrix< DDRMat > compute_dcriteria_dadv() = 0;
    };
}
