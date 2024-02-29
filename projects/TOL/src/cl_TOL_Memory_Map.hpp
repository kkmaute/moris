/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TOL_Memory_Map.hpp
 *
 */

#ifndef SRC_TOOLS_CL_TOL_MEMORY_MAP_HPP_
#define SRC_TOOLS_CL_TOL_MEMORY_MAP_HPP_

#include <unordered_map>

#include "cl_Matrix.hpp"

#include "fn_TOL_Capacities.hpp"

namespace moris
{
    /*!
    * This class is very slim. It is a safe guard class
    * so I can define operators without them being
    * directly defined
    */
    class Memory_Map
    {
    public:
        // ----------------------------------------------------------------------------------

        Memory_Map();

        // ----------------------------------------------------------------------------------

        Memory_Map( Vector<std::string> const & aKeys,
                   Matrix<DDSTMat>   const & aVals);

        // ----------------------------------------------------------------------------------

        ~Memory_Map();

        // ----------------------------------------------------------------------------------
       /*!
        * @brief Print the memory usage of the map
        */
        void
        print(std::string const & aTitle);

        // ----------------------------------------------------------------------------------

       /*!
        * @brief Parallel print the memory usage of the map
        */
        void
        par_print( std::string const & aTitle );

        // ----------------------------------------------------------------------------------

        /*
        * @brief Add Memory maps together. Data with same key is combined
        */
        Memory_Map
        operator+( const Memory_Map& aMemMapB );

        // ----------------------------------------------------------------------------------

        /*
        * @brief Sum up the memory in memory map map
        */
        size_t
        sum( );

        // ----------------------------------------------------------------------------------

        std::unordered_map<std::string, size_t> mMemoryMapData;

        // ----------------------------------------------------------------------------------

        private:

        // ----------------------------------------------------------------------------------

        /*!
         * @brief Gather all memory maps to root processor.
         * @param[out] aGatheredMemMap Gathered Memory Maps;
         */
        void
        gather_all( Vector<Memory_Map> & aGatheredMemMap);

        // ----------------------------------------------------------------------------------
        /*!
         * @brief Serialize memory map
         * mpi purposes.
         * @param[out] aKeyCell Cell of keys (on root proc)
         * @param[out] aValCell Mat of vals (on root proc)
         */
        void
        serialize( Vector<std::string> & aKeyCell,
                   Matrix<DDSTMat>   & aValCell);

        /*!
         * @brief Deserialize memory maps on root processor
         * @param[out]
         */
        void
        deserialize( Vector< Vector<std::string>> & aGatheredKeyCells,
                Vector<Matrix<DDSTMat>>   & aGatheredValCells,
                Vector<Memory_Map>        & aGatheredMemMaps);

    };
} // namespace moris

#endif /* SRC_TOOLS_CL_TOL_MEMORY_MAP_HPP_ */

