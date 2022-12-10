/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Meshgen.hpp
 *
 */

#ifndef MORIS_CL_LIBRARY_IO_MESHGEN_HPP
#define MORIS_CL_LIBRARY_IO_MESHGEN_HPP

#include <set>
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

      protected:

        // -----------------------------------------------------------------------------

        // list of supported parameter list types
        std::set< Parameter_List_Type > mSupportedParamListTypes = { Parameter_List_Type::HMR, Parameter_List_Type::XTK, Parameter_List_Type::GEN };

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
        ~Library_IO_Meshgen();

        // -----------------------------------------------------------------------------

        /**
         * @brief finishes the initialization of the library and locks it from modification
         */
        virtual
        void
        finalize();

        // -----------------------------------------------------------------------------

        /**
         * @brief fills the member parameter lists with the standard parameters for all modules
         */
        virtual
        void
        load_all_standard_parameters();

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the standard parameters for the specified module in the mesh generation only workflow
         * note: this function deletes all previous entries in the parameter list passed to it to fill
         * 
         * @param aParamListType module for which the standar parameters are to be provided
         * @param aParameterList address to the parameter list to fill with standard parameters
         */
        void
        create_standard_parameter_list_for_module( 
                Parameter_List_Type aParamListType,
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

        void
        create_standard_geometry_parameter_list( ParameterList & aParameterList );

        void
        create_standard_user_defined_geometry_parameter_list( ParameterList & aParameterList );

        void
        create_standard_voxel_field_parameter_list( ParameterList & aParameterList );

        void
        create_standard_image_sdf_field_parameter_list( ParameterList & aParameterList );

        // -----------------------------------------------------------------------------

    }; // class Library_IO_Meshgen

// -----------------------------------------------------------------------------

}    // namespace moris

#endif //MORIS_CL_LIBRARY_IO_MESHGEN_HPP

