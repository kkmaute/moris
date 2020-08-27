/*
 * cl_TOL_Memory_Map.hpp
 *      Author: doble
 */

#ifndef SRC_TOOLS_CL_TOL_MEMORY_MAP_HPP_
#define SRC_TOOLS_CL_TOL_MEMORY_MAP_HPP_

#include <unordered_map>

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

        ~Memory_Map();

        // ----------------------------------------------------------------------------------
       /*! 
        * @brief Print the memory usage of the map
        */
        void
        print_memory_map();

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
    };
} // namespace moris

#endif /* SRC_TOOLS_CL_TOL_MEMORY_MAP_HPP_ */
