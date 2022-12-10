/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Standard.hpp
 *
 */

#ifndef MORIS_CL_LIBRARY_IO_STANDARD_HPP
#define MORIS_CL_LIBRARY_IO_STANDARD_HPP

#include <set>
#include "cl_Library_IO.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    /**
     * Wrapper class for shared object file
     */
    class Library_IO_Standard : public Library_IO
    {
        // -----------------------------------------------------------------------------

      protected:

        // -----------------------------------------------------------------------------

        // list of supported parameter list types
        std::set< Parameter_List_Type > mSupportedParamListTypes = {
            Parameter_List_Type::OPT,
            Parameter_List_Type::HMR,
            Parameter_List_Type::STK,
            Parameter_List_Type::XTK,
            Parameter_List_Type::GEN,
            Parameter_List_Type::FEM,
            Parameter_List_Type::SOL,
            Parameter_List_Type::MSI,
            Parameter_List_Type::VIS,
            Parameter_List_Type::MIG,
            Parameter_List_Type::WRK,
            Parameter_List_Type::MORISGENERAL };

        // -----------------------------------------------------------------------------

      public:

        // -----------------------------------------------------------------------------

        /**
         * Default constructor
         */
        Library_IO_Standard();

        // -----------------------------------------------------------------------------

        /**
         * Default destructor
         */
        virtual ~Library_IO_Standard();

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

    }; // class Library_IO_Standard

// -----------------------------------------------------------------------------

} // namespace moris

#endif    // MORIS_CL_LIBRARY_IO_STANDARD_HPP
