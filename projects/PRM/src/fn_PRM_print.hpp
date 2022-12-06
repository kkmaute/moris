/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_SOL_Parameters.hpp
 *
 */

#ifndef MORIS_FN_PRM_PRINT_HPP
#define MORIS_FN_PRM_PRINT_HPP

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {
        //------------------------------------------------------------------------------

        // types in parameter list
        // bool, sint, real, const char*, std::string, uint, std::pair< std::string, std::string >
        inline void
        print( ParameterList& aParamList )
        {
            for ( auto it = aParamList.begin(); it != aParamList.end(); ++it )
            {
                std::cout << std::setw( 40 ) << std::left << it->first << " | ";

                if ( boost::get< bool >( &( it->second ) ) != nullptr )
                {
                    std::cout << boost::get< bool >( it->second ) << "\n";
                }

                else if ( boost::get< sint >( &( it->second ) ) != nullptr )
                {
                    std::cout << boost::get< sint >( it->second ) << "\n";
                }

                else if ( boost::get< real >( &( it->second ) ) != nullptr )
                {
                    std::cout << boost::get< real >( it->second ) << "\n";
                }

                else if ( boost::get< const char* >( &( it->second ) ) != nullptr )
                {
                    std::cout << boost::get< const char* >( it->second ) << "\n";
                }

                else if ( boost::get< std::string >( &( it->second ) ) != nullptr )
                {
                    std::cout << "\"" << boost::get< std::string >( it->second ) << "\""
                              << "\n";
                }

                else if ( boost::get< uint >( &( it->second ) ) != nullptr )
                {
                    std::cout << boost::get< uint >( it->second ) << "\n";
                }

                else
                {
                    std::cout << " not printable" << std::endl;
                }
            }
        }

        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif    // MORIS_FN_PRM_PRINT_HPP
