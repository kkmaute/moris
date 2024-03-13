/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Interface_User_Defined.cpp
 *
 */

#include "cl_OPT_Interface_User_Defined.hpp"
#include "cl_Library_Factory.hpp"

namespace moris
{
    namespace opt
    {
        //--------------------------------------------------------------------------------------------------------------

        Interface_User_Defined::Interface_User_Defined( Parameter_List aParameterList )
        {
            // Load library
            moris::Library_Factory tLibraryFactory;
            mLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
            std::string tLibraryName = aParameterList.get< std::string >( "library" );
            mLibrary->load_parameter_list( tLibraryName, File_Type::SO_FILE );
            mLibrary->finalize();

            // Set user-defined functions
            initialize_user_defined             = mLibrary->load_function< Criteria_Initialize_Function >( "initialize" );
            get_criteria_user_defined           = mLibrary->load_function< Criteria_Function >( "get_criteria" );
            compute_dcriteria_dadv_user_defined = mLibrary->load_function< Criteria_Function >( "get_dcriteria_dadv" );
        }

        //--------------------------------------------------------------------------------------------------------------

        Interface_User_Defined::Interface_User_Defined(
                Criteria_Initialize_Function aInitializationFunction,
                Criteria_Function            aCriteriaEvaluationFunction,
                Criteria_Function            aCriteriaGradientFunction )
        {
            // Set user-defined functions
            initialize_user_defined             = aInitializationFunction;
            get_criteria_user_defined           = aCriteriaEvaluationFunction;
            compute_dcriteria_dadv_user_defined = aCriteriaGradientFunction;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Interface_User_Defined::initialize(
                Matrix< DDRMat >& aADVs,
                Matrix< DDRMat >& aLowerBounds,
                Matrix< DDRMat >& aUpperBounds,
                Matrix< IdMat >& )
        {
            initialize_user_defined( aADVs, aLowerBounds, aUpperBounds );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Interface_User_Defined::perform( Matrix< DDRMat >& aNewADVs )
        {
            mADVs = aNewADVs;

            return this->get_criteria_user_defined( mADVs );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Interface_User_Defined::compute_dcriteria_dadv()
        {
            return this->compute_dcriteria_dadv_user_defined( mADVs );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace opt
}    // namespace moris
