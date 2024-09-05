/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Phase_Table.cpp
 *
 */

#include "cl_GEN_Phase_Table.hpp"
#include "fn_linspace.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Phase_Table::Phase_Table(
            uint                         aNumGeometries,
            const Matrix< DDUMat >&      aBulkPhases,
            const Vector< std::string >& aPhaseNames )
            : mNumGeometries( aNumGeometries )
            , mBulkPhases( aBulkPhases )
            , mPhaseNames( aPhaseNames )
    {
        if ( mBulkPhases.numel() > 0 )
        {
            // Number of phases
            mNumPhases = mBulkPhases.max() + 1;

            // Default phase names
            if ( mPhaseNames.size() == 0 )
            {
                this->set_default_phase_names();
            }

            // Check for phase table size
            MORIS_ERROR( mBulkPhases.length() == std::pow( 2, mNumGeometries ),
                    "Must provide bulk phase information for each of the 2^n geometry combinations." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Phase_Table::Phase_Table(
            uint                         aNumGeometries,
            const Vector< std::string >& aPhaseNames )
            : Phase_Table(
                      aNumGeometries,
                      linspace< uint >( 0, std::pow( 2, aNumGeometries ) - 1, std::pow( 2, aNumGeometries ) ),
                      aPhaseNames )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Phase_Table::Phase_Table(
            PHASE_FUNCTION               aPhaseFunction,
            uint                         aNumPhases,
            const Vector< std::string >& aPhaseNames )
            : Phase_Table( 0, { { aNumPhases - 1 } }, aPhaseNames )
    {
        mPhaseFunction = aPhaseFunction;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Phase_Table::get_num_phases()
    {
        return mNumPhases;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Phase_Table::get_phase_index( const Geometry_Bitset& aGeometrySigns )
    {
        if ( mPhaseFunction )
        {
            return mPhaseFunction( aGeometrySigns );
        }
        else
        {
            uint tPhaseEntry = 0;
            for ( uint tGeometryIndex = 0; tGeometryIndex < mNumGeometries; tGeometryIndex++ )
            {
                tPhaseEntry += std::pow( 2, mNumGeometries ) / std::pow( 2, tGeometryIndex + 1 ) * aGeometrySigns.test( tGeometryIndex );
            }

            return mBulkPhases( tPhaseEntry );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::string
    Phase_Table::get_phase_name( uint aPhaseIndex )
    {
        MORIS_ASSERT( aPhaseIndex < mPhaseNames.size(), "Phase index out of bounds" );
        return mPhaseNames( aPhaseIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Phase_Table::print()
    {
        // Print table information
        std::cout << "Phase Table Info:" << '\n';
        std::cout << "  Number of Geometries:  " << mNumGeometries << '\n';
        std::cout << "  Number of Bulk Phases: " << mNumPhases << '\n';

        // Print table header
        std::cout << std::setw( 8 ) << "Entry"
                  << " | "
                  << std::setw( 8 ) << "BP Index"
                  << " | "
                  << std::setw( 8 ) << "BP Name"
                  << " | ";

        // Print geometry index
        for ( uint tGeometryIndex = 0; tGeometryIndex < mNumGeometries; tGeometryIndex++ )
        {
            std::cout << std::setw( 8 ) << "Geom " + std::to_string( tGeometryIndex ) << " | ";
        }
        std::cout << '\n';

        // Print bulk phase information
        for ( uint tPhaseEntry = 0; tPhaseEntry < mBulkPhases.length(); tPhaseEntry++ )
        {
            // Print phase entry
            std::cout << std::setw( 8 ) << tPhaseEntry << " | "
                      << std::setw( 8 ) << mBulkPhases( tPhaseEntry ) << " | "
                      << std::setw( 8 ) << mPhaseNames( mBulkPhases( tPhaseEntry ) ).substr( 0, 8 ) << " | ";

            // Geometry signs (reversed)
            Geometry_Bitset tGeometrySigns( tPhaseEntry );

            // Print geometry signs
            for ( uint tGeometryIndex = mNumGeometries; tGeometryIndex-- > 0; )
            {
                if ( tGeometrySigns.test( tGeometryIndex ) )
                {
                    std::cout << std::setw( 8 ) << "+"
                              << " | ";
                }
                else
                {
                    std::cout << std::setw( 8 ) << "-"
                              << " | ";
                }
            }
            std::cout << '\n';
        }

        std::cout << "  +   -> greater than threshold" << '\n';
        std::cout << "  -   -> less than threshold" << '\n';
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Phase_Table::set_default_phase_names()
    {
        mPhaseNames       = Vector< std::string >( mNumPhases );
        std::string tBase = "p_";

        for ( uint tPhaseIndex = 0; tPhaseIndex < mNumPhases; tPhaseIndex++ )
        {
            mPhaseNames( tPhaseIndex ) = tBase + std::to_string( tPhaseIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Phase_Table::set_phase_function(
            PHASE_FUNCTION               aPhaseFunction,
            uint                         aNumPhases,
            const Vector< std::string >& aPhaseNames )
    {
        mNumGeometries = 0;
        mBulkPhases    = { { aNumPhases - 1 } };
        mPhaseNames    = aPhaseNames;
        mPhaseFunction = aPhaseFunction;
    }

}
