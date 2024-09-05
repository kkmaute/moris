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
         * @brief finishes the initialization of the library and locks it from modification
         */
        void
        finalize( std::string aFilePath );

        // -----------------------------------------------------------------------------

        /**
         * @brief fills the member parameter lists with the standard parameters for all modules
         */
        void
        load_all_standard_parameters() override;

        // -----------------------------------------------------------------------------

    };    // class Library_IO_Standard

    // -----------------------------------------------------------------------------

}    // namespace moris

#endif    // MORIS_CL_LIBRARY_IO_STANDARD_HPP
