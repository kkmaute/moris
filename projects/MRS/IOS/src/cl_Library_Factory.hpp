/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_Library_Factory.hpp  
 * 
 */
#ifndef SRC_cl_Library_Factory
#define SRC_cl_Library_Factory

#include "cl_Library_Enums.hpp"
#include "cl_Library_IO.hpp"
#include "cl_Library_IO_Meshgen.hpp"
#include "cl_Library_IO_Standard.hpp"

namespace moris
{

    // -----------------------------------------------------------------------------

    /**
     * @brief factory for creating Library objects for the various input file options
     */
    class Library_Factory
    {
      public:

        // -----------------------------------------------------------------------------

        /**
         * @brief Default constructor doing nothing.
         * 
         */
        Library_Factory(){};

        // -----------------------------------------------------------------------------

        /**
         * @brief Default destructor just deleting the object
         * 
         */
        ~Library_Factory(){};

        // -----------------------------------------------------------------------------

        /**
         * @brief Create a shared pointer to a specialized Library 
         * casted into the type of the parent class
         * 
         * @param aLibraryType library type to be created
         * @return std::shared_ptr< Library_IO > 
         */
        std::shared_ptr< Library_IO >
        create_Library( Library_Type aLibraryType )
        {
            switch ( aLibraryType )
            {
                case Library_Type::STANDARD :
                    return std::make_shared< Library_IO_Standard >();

                case Library_Type::MESHGEN :
                    return std::make_shared< Library_IO_Meshgen >();

                default :
                    MORIS_ERROR( false, "Library_Factory::create_Library() - Library type specified is not defined. " );
                    return nullptr;
            }
        }

        // -----------------------------------------------------------------------------

    }; // class Library_Factory

} // namespace moris

#endif /* cl_Library_Factory.hpp */
