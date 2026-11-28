/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#ifndef MORIS_FN_EXA_TEST_UTILITIES_HPP
#define MORIS_FN_EXA_TEST_UTILITIES_HPP

#include <filesystem>
#include <initializer_list>
#include <string>
#include <system_error>

#include "cl_Communication_Tools.hpp"
#include "fn_assert.hpp"

namespace moris::exa
{
    // Remove only the Exodus/HDF5 outputs explicitly owned by the upcoming
    // workflow. Names use the serial form; each rank expands its own suffix.
    // Restart inputs are deliberately absent from the restart stage's list.
    inline void remove_test_outputs( std::initializer_list< std::string > aOutputs )
    {
        bool tFailed = false;
        for ( const std::string& tOutput : aOutputs )
        {
            std::string tPath = tOutput;
            if ( par_size() > 1 )
            {
                if ( std::filesystem::path( tOutput ).extension() == ".hdf5" )
                {
                    const auto tExtension = tOutput.find_last_of( '.' );
                    tPath                 = tOutput.substr( 0, tExtension ) + "_" + std::to_string( par_size() )
                          + "." + std::to_string( par_rank() ) + tOutput.substr( tExtension );
                }
                else
                {
                    std::string tRank = std::to_string( par_rank() );
                    tRank.insert( 0, std::to_string( par_size() ).size() - tRank.size(), '0' );
                    tPath += "." + std::to_string( par_size() ) + "." + tRank;
                }
            }

            std::error_code tError;
            std::filesystem::remove( tPath, tError );
            tFailed = tFailed || static_cast< bool >( tError );
        }

        // All ranks finish removing their own output before any starts writing.
        // Report errors collectively so no rank waits in the next workflow alone.
        MORIS_ERROR( max_all( static_cast< sint >( tFailed ) ) == 0,
                "EXA test could not remove an explicitly owned output file." );
    }
}    // namespace moris::exa

#endif
