/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Interface_User_Defined.hpp
 *
 */

#pragma once

#include "cl_OPT_Criteria_Interface.hpp"

// Forward declarations
namespace moris
{
    class Parameter_List;
    class Library_IO;
}


namespace moris::opt
{
    // User-defined function types
    typedef void ( *Criteria_Initialize_Function ) (
            Vector< real >& aAdvs,
            Vector< real >& aLowerBounds,
            Vector< real >& aUpperBounds);

    typedef Vector< real > ( *Criteria_Function ) (const moris::Vector< real >& );
    typedef Matrix< DDRMat > ( *Criteria_Gradient_Function ) (const moris::Vector< real >& );

    class Interface_User_Defined : public Criteria_Interface
    {
    private:

        Vector< real > mADVs;
        std::shared_ptr<Library_IO> mLibrary;

    public:

        /**
         * Constructor
         *
         * @param aParameterList Parameter list containing parameters for a user-defined interface.
         */
      explicit Interface_User_Defined( const Parameter_List& aParameterList );

      /**
       * Alternate constructor where the user-defined functions are provided directly. Used in the OPT tests.
       *
       * @param aInitializationFunction Function for initializing ADVs and lower/upper bounds.
       * @param aCriteriaEvaluationFunction Function for evaluating the criteria vector.
       * @param aCriteriaGradientFunction Function for evaluating the gradient of the criteria vector wrt ADVs.
       */
      Interface_User_Defined(
              Criteria_Initialize_Function aInitializationFunction,
              Criteria_Function            aCriteriaEvaluationFunction,
              Criteria_Gradient_Function   aCriteriaGradientFunction );

      /**
       * Initializes the vectors of ADV values, lower bounds, and upper bounds
       *
       * @param aADVs Initial ADVs to be filled.
       * @param aLowerBounds Lower ADV bounds to be filled.
       * @param aUpperBounds Upper ADV bounds to be filled.
       */
      void initialize(
              Vector< real >&  aADVs,
              Vector< real >&  aLowerBounds,
              Vector< real >&  aUpperBounds,
              Matrix< IdMat >& aIjklIds ) override;

      /**
       * Gets the criteria values.
       *
       * @return vector of criteria
       */
      Vector< real > perform( Vector< real >& aNewADVs ) override;

      /**
       * Gets the derivative of the criteria with respect to the advs.
       *
       * @return matrix d(criteria)_i/d(adv)_j
       */
      Matrix< DDRMat > compute_dcriteria_dadv() override;

    private:

        // Loaded user-defined functions
        Criteria_Initialize_Function initialize_user_defined  = nullptr;
        Criteria_Function get_criteria_user_defined           = nullptr;
        Criteria_Gradient_Function   compute_dcriteria_dadv_user_defined = nullptr;
    };
}
