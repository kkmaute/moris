/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TOL_Memory_Map.cpp
 *
 */

#include "cl_TOL_Memory_Map.hpp"
#include <iostream>
#include <iomanip>
#include "moris_typedefs.hpp"
#include "cl_Communication_Tools.hpp"
#include "assert.hpp"
namespace moris
{
    Memory_Map::Memory_Map()
    {
    }

    // ----------------------------------------------------------------------------------

    Memory_Map::Memory_Map(Cell<std::string> const & aKeys,
                           Matrix<DDSTMat>   const & aVals)
    {
        MORIS_ERROR(aKeys.size() == aVals.numel(), "Dimension mismatch on allocation.");
        for(moris::uint i = 0; i < aKeys.size() ; i++)
        {
            MORIS_ERROR(mMemoryMapData.find(aKeys(i)) == mMemoryMapData.end(),"Duplicate key");
            mMemoryMapData[aKeys(i)] = aVals(i);
        }
    }

    // ----------------------------------------------------------------------------------

    Memory_Map::~Memory_Map()
    {
    }

    // ----------------------------------------------------------------------------------
    void
    Memory_Map::print(std::string const & aTitle)
    {
        std::cout<<"\n----------------------------------------------------------------------------------\n";
        std::cout<<" Memory Map Name: "<<aTitle<<"\n";
        size_t tTotal = 0;
        size_t tWidth = 0;
        moris::real tTotalPercent = 0.0;

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
                      << (moris::real) it->second / (1000)  // string's value
                      << " KiB | "
                      << std::setw(12)
                      << (moris::real)it->second / (moris::real)tTotal * 100
                      << "%"
                      << std::endl;
            tTotalPercent = tTotalPercent + (moris::real)it->second / (moris::real)tTotal * 100;
        }
        std::cout<<"----------------------------------------------------------------------------------\n";
        std::cout     << std::left
                      << std::setw(tWidth + 1)
                      <<" "
                      << " | "
                      << std::right
                      << std::setw(15)
                      << tTotal/ (1000)
                      << " KiB | "
                      << std::setw(12)
                      << tTotalPercent
                      << "%"
                      << std::endl;
    }

    // ----------------------------------------------------------------------------------

    void
    Memory_Map::par_print(std::string const & aTitle)
    {
        // get all the memory maps onto this proc
        Cell<Memory_Map> tGatheredMM;
        this->gather_all(tGatheredMM);

        // Combined memory map
        Memory_Map tFullMM;

        if(par_rank() == 0)
        {
            // total memory
            size_t tTotalMem = 0;
            // determine memory on each proc
            Cell<size_t> tTotalMemPerProc(par_size());
            for (int i = 0; i < par_size(); i++)
            {
                tTotalMemPerProc(i) = tGatheredMM(i).sum();
                tTotalMem += tTotalMemPerProc(i);
                tFullMM = tFullMM + tGatheredMM(i);
            }
            moris::real tTotalPercent = 0.0;

            // print the full map
            tFullMM.print(aTitle);

            // Header
            std::cout << std::setw(16) << "Proc Rank"
                      << " | "
                      << std::setw(16) << "Memory (KiB)"
                      << " | "
                      << std::setw(12) << "     " << std::endl;

            for (moris::uint i = 0; i < tTotalMemPerProc.size(); i++)
            {
                std::cout << std::setw(16) << i
                          << " | "
                          << std::setw(16) << tTotalMemPerProc(i) / 1000
                          << " | "
                          << std::setw(12) << (real)tTotalMemPerProc(i) / (real)tTotalMem * 100
                          << "%"
                          << std::endl;
                tTotalPercent += (real)tTotalMemPerProc(i) / (real)tTotalMem * 100;
            }

            std::cout << "----------------------------------------------------------------------------------\n";
            std::cout << std::left
                      << std::setw(16)
                      << " "
                      << " | "
                      << std::right
                      << std::setw(12)
                      << tTotalMem / 1000
                      << " KiB | "
                      << std::setw(12)
                      << tTotalPercent
                      << "%"
                      << std::endl;
        }
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

    void
    Memory_Map::gather_all(Cell<Memory_Map> & aGatheredMemMap)
    {
        // serialize my map as a cell of keys (str) and a cell of size_t (mem)
        Cell<std::string> tKeyCell;
        Matrix<DDSTMat>   tVals;
        this->serialize(tKeyCell,tVals);

        Cell<Cell<std::string>> tGatheredKeys;
        moris_index tTag = 7000; // no particular reason for number
        all_gather_cell_of_str(tKeyCell,tGatheredKeys,tTag);

        // all gather the values
        Cell<Matrix<DDSTMat>> tGatheredVals;
        tTag = 7001; // no particular reason for number
        all_gather_vector(tVals,tGatheredVals,tTag,0,0);

        if(par_rank() == 0)
        {
            this->deserialize(tGatheredKeys,tGatheredVals,aGatheredMemMap);
        }

    }

    // ----------------------------------------------------------------------------------

    void
    Memory_Map::serialize(Cell<std::string> & aKeyCell,
                           Matrix<DDSTMat>  & aValCell)
    {
        // size the serialized vars
        size_t tNumEntries = mMemoryMapData.size();
        aKeyCell.reserve(tNumEntries);
        aValCell.resize(1,tNumEntries);

        // iterate through map and grab the key and val
        uint tCount = 0;
        for(auto it = mMemoryMapData.begin(); it != mMemoryMapData.end(); it++)
        {
            aKeyCell.push_back(it->first);
            aValCell(tCount++) = it->second;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Memory_Map::deserialize(Cell<Cell<std::string>> & aGatheredKeyCells,
                            Cell<Matrix<DDSTMat>>   & aGatheredValCells,
                            Cell<Memory_Map>        & aGatheredMemMaps)
    {
        aGatheredMemMaps.resize(par_size());
        for(moris::uint i = 0; i < aGatheredKeyCells.size(); i++)
        {
            aGatheredMemMaps(i) = Memory_Map(aGatheredKeyCells(i),aGatheredValCells(i));
        }
    }

} // namespace moris

