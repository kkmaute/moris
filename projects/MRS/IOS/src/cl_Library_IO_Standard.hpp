/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Standard.hpp
 *
 */

#pragma once

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
        ~Library_IO_Standard() override;

        // -----------------------------------------------------------------------------

        /**
         * Gets a list of enums representing all of the modules being used, and requiring filled parameter lists
         *
         * @return Vector of used module enums
         */
        bool is_module_supported( Module_Type aModuleType ) override;

        // -----------------------------------------------------------------------------

    };    // class Library_IO_Standard

    // -----------------------------------------------------------------------------

}    // namespace moris
