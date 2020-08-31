
#include "cl_TOL_Memory_Map.hpp"
#include <iostream>
#include <iomanip>
#include "typedefs.hpp"
namespace moris
{
    Memory_Map::Memory_Map()
    {
    }

    // ----------------------------------------------------------------------------------

    Memory_Map::~Memory_Map()
    {
    }

    // ----------------------------------------------------------------------------------
    void
    Memory_Map::print()
    {
        std::cout<<"\n----------------------------------------------------------------------------------\n";
        size_t tTotal = 0;
        size_t tWidth = 0;
        moris::real tTotalPercent= 0.0;

        // iterate and collect the total size
        for (auto it = mMemoryMapData.begin(); it != mMemoryMapData.end(); it++)
        {
            tTotal = tTotal + it->second;
            if (it->first.length() > tWidth)
            {
                tWidth = it->first.length();
            }
        }

        for (auto it = mMemoryMapData.begin(); it != mMemoryMapData.end(); it++)
        {
            std::cout << std::left
                      << std::setw(tWidth + 1)
                      << it->first // string (key)
                      << " | "
                      << std::right
                      << std::setw(15)
                      << (moris::real) it->second / (1024*8)  // string's value
                      << " KiB | "
                      << std::setw(12)
                      << (moris::real)it->second / (moris::real)tTotal * 100
                      << "%"
                      << std::endl;
            tTotalPercent = tTotalPercent + (moris::real)it->second / (moris::real)tTotal * 100;
        }

        std::cout     << std::left
                      << std::setw(tWidth + 1)
                      <<" "
                      << " | "
                      << std::right
                      << std::setw(15)
                      << tTotal/ (1024*8) 
                      << " KiB | "
                      << std::setw(12)
                      << tTotalPercent
                      << "%"
                      << std::endl;
    }

    // ----------------------------------------------------------------------------------
 
    Memory_Map
    Memory_Map::operator+(const Memory_Map& aMemMapB)
    {
        Memory_Map tCombinedMap;
        tCombinedMap.mMemoryMapData = mMemoryMapData;


        for (auto it = aMemMapB.mMemoryMapData.begin(); it != aMemMapB.mMemoryMapData.end(); it++)
        {
            if(tCombinedMap.mMemoryMapData.find(it->first) == tCombinedMap.mMemoryMapData.end())
            {
                tCombinedMap.mMemoryMapData[it->first] = 0;
            }

            tCombinedMap.mMemoryMapData[it->first] = tCombinedMap.mMemoryMapData[it->first] + it->second;
        }
        return tCombinedMap;
    }

    // ----------------------------------------------------------------------------------
 
    size_t
    Memory_Map::sum( )
    {
        size_t tTotal = 0;
        for (auto it = mMemoryMapData.begin(); it != mMemoryMapData.end(); it++)
        {
            tTotal = tTotal + it->second;
        }

        return tTotal;
    }

    // ----------------------------------------------------------------------------------
    
 
} // namespace moris