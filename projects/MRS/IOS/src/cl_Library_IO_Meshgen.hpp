/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Meshgen.hpp
 *
 */

#pragma once

#include "cl_Library_IO.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    /**
     * Class providing standard input parameters for generating just the mesh information
     * and additional parameters specified through an input file
     */
    class Library_IO_Meshgen : public Library_IO
    {
        // -----------------------------------------------------------------------------

      private:

        // store number of spatial dimensions
        uint mNumSpatialDims = 0;

        // determine mesh refinements in background mesh parameterlist, but apply in GEN, hence store it
        Vector< uint > mGenRefineMeshIndices;
        Vector< uint > mGenNumRefinements;

        // -----------------------------------------------------------------------------

      public:
        // -----------------------------------------------------------------------------

        /**
         * Default constructor
         */
        Library_IO_Meshgen();

        // -----------------------------------------------------------------------------

        /**
         * Default destructor
         */
        ~Library_IO_Meshgen() override;

        /**
         * @brief loads parameters from the simplified xml file and overwrites any previously specified parameters by it
         */
        void
        load_parameters_from_xml() override;

        // -----------------------------------------------------------------------------

        void
        load_HMR_parameters_from_xml(
                std::string const & aHmrPath,
                std::string const & aXtkPath );

        void
        load_XTK_parameters_from_xml(
                std::string const & aXtkPath,
                std::string const & aHmrPath );

        void
        load_GEN_parameters_from_xml(
                std::string const & aGenPath,
                std::string const & aHmrPath,
                std::string const & aXtkPath );

        // -----------------------------------------------------------------------------

        /**
         * @brief fills the member parameter lists with the standard parameters for all modules
         */
        void
        load_all_standard_parameters() override;

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the standard parameters for the specified module in the mesh generation only workflow
         * note: this function deletes all previous entries in the parameter list passed to it to fill
         *
         * @param aParamListType module for which the standard parameters are to be provided
         * @param aParameterList address to the parameter list to fill with standard parameters
         */
        void
        create_standard_parameter_list_for_module(
                Parameter_List_Type  aParamListType,
                ModuleParameterList& aParameterList );

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the standard parameters for OPT in the mesh generation only workflow
         * note: this function deletes all previous entries in the parameter list passed to it to fill
         *
         * @param aParameterList address to the parameter list to fill with standard parameters
         */
        void
        create_standard_OPT_parameter_list( ModuleParameterList& aParameterList );

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the standard parameters for XTK in the mesh generation only workflow
         * note: this function deletes all previous entries in the parameter list passed to it to fill
         *
         * @param aParameterList address to the parameter list to fill with standard parameters
         */
        void
        create_standard_XTK_parameter_list( ModuleParameterList& aParameterList );

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the standard parameters for HMR in the mesh generation only workflow
         * note: this function deletes all previous entries in the parameter list passed to it to fill
         *
         * @param aParameterList address to the parameter list to fill with standard parameters
         */
        void
        create_standard_HMR_parameter_list( ModuleParameterList& aParameterList );

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the standard parameters for GEN in the mesh generation only workflow
         * note: this function deletes all previous entries in the parameter list passed to it to fill
         *
         * @param aParameterList address to the parameter list to fill with standard parameters
         */
        void
        create_standard_GEN_parameter_list( ModuleParameterList& aParameterList );

        // -----------------------------------------------------------------------------

    };    // class Library_IO_Meshgen

    // -----------------------------------------------------------------------------

}    // namespace moris
