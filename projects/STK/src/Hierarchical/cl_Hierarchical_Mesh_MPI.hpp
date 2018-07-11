/*
 * cl_Hierarchical_Mesh_MPI.hpp
 *
 *  Created on: Jan 4, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_MPI_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_MPI_HPP_

#include "mpi.h"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_BoostBitset.hpp" // CON/src
#include <boost/dynamic_bitset.hpp>

namespace moris
{

    class Hierarchical_Mesh_MPI
    {
    public:
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_MPI()
    {
    }
        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_MPI() = default;
        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha only those, where both vectors have at the same position a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aBitset       Bitset with active elements or basis functions
         *
         * @param[out] aBitset      Bitset with active elements or basis functions after comparing with all procs
         *
         */
        void broadcast_bitset_logical_and(
                BoostBitset & aBitset);

        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha all bits with a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aBitset       Bitset with active elements or basis functions
         *
         * @param[out] aBitset      Bitset with active elements or basis functions after comparing with all procs
         *
         */
        void broadcast_bitset_logical_or(
                BoostBitset & aBitset);

        /**
         * Gather a message to proc 0, take the max of the proc values and broadcast the max.
         *
         * @param[in] aElementList       aElementList.Level
         *
         * @param[out] aElementList      aElementList.Level (updated)
         *
         */
        void gather_value_and_bcast_max(
                uint & aMessage);
    };
}
#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_MPI_HPP_ */
