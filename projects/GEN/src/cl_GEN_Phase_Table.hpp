/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Phase_Table.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_Library_IO.hpp"

namespace moris::gen
{
    // User-defined phase function
    typedef uint ( *PHASE_FUNCTION )( const gen::Geometry_Bitset& aGeometrySigns );

    class Phase_Table
    {

      private:
        uint mNumGeometries = 0;    // Total number of geometries
        uint mNumPhases     = 0;    // Number of bulk phases

        Matrix< DDUMat >    mBulkPhases;    // Geometric sign to bulk phase
        Vector< std::string > mPhaseNames;    // Phase names

        PHASE_FUNCTION mPhaseFunction = nullptr;

      public:
        /**
         * Constructor for using a given phase table with the standard 2^n structure.
         *
         * @param aNumGeometries Number of geometries
         * @param aBulkPhases Geometric index to bulk phase map
         * @param aPhaseNames (optional) Phase names
         */
        Phase_Table(
                uint                aNumGeometries,
                Matrix< DDUMat >    aBulkPhases,
                Vector< std::string > aPhaseNames = {} );

        /**
         * Create a phase table with 2^n structure using a number of bulk phases. Delegating constructor.
         *
         * @param aNumGeometries Number of geometries
         * @param aPhaseNames (optional) Phase names
         */
        Phase_Table(
                uint                aNumGeometries,
                Vector< std::string > aPhaseNames = {} );

        /**
         * Create a phase table where the phase indices are decided by a user-defined function.
         *
         * @param aPhaseFunction User-defined phase function
         * @param aNumPhases Number of different bulk phases that the phase function can return
         * @param aPhaseNames (optional) Phase names
         */
        Phase_Table(
                PHASE_FUNCTION      aPhaseFunction,
                uint                aNumPhases,
                Vector< std::string > aPhaseNames = {} );

        /**
         * Get the number of phases
         *
         * @return Number of phases
         */
        uint get_num_phases();

        /**
         * Get phase index based on entity phase info
         *
         * @param aGeometrySigns Geometry sign info
         * @return Phase index
         */
        uint get_phase_index( const Geometry_Bitset& aGeometrySigns );

        /**
         * Gets the name of a requested phase
         *
         * @param aPhaseIndex The index of the requested phase
         * @return Phase name
         */
        std::string get_phase_name( uint aPhaseIndex );

        /*!
         * Set the phase table without a library
         */
        void
        set_phase_function( PHASE_FUNCTION aPhaseFunction,
                uint                       aNumPhases,
                Vector< std::string >        aPhaseNames = {} );

        /*!
         * Print information for setting up phase table
         */
        void print();

      private:
        /**
         * Set all of the phases to have default names (p_i)
         */
        void set_default_phase_names();
    };
}    // namespace moris
