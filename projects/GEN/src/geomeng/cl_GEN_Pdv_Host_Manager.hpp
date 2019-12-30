/*
 * cl_GEN_Pdv_Host_Manager.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_

#include <unordered_map>

#include "cl_GEN_Pdv_Host.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv_Host_Manager
        {

        public:
            Pdv_Host_Manager(){};

            ~Pdv_Host_Manager(){};
            //------------------------------------------------------------------------------
            /*
             * Stores pdv host objects in the pdv host manager associated with nodes
             */
            void store_pdv_hosts( moris::Matrix< moris::IndexMat > const & aNodeIndices,
                             moris::Cell< GEN_Pdv_Host >      const & aGeometryObjects )
            {
                moris::size_t tNumExistingPdvHosts = mPdvHosts.size();
                moris::size_t tNumNewPdvHpsts = aNodeIndices.n_cols();

                MORIS_ASSERT(tNumNewPdvHpsts == aGeometryObjects.size(),"cl_GEN_Pdv_Host_Manager()::store_pdv_hosts() - number of pdv hosts does not match number of node indices provided");

                // append the geometry object cell
                mPdvHosts.append(aGeometryObjects);

                // Add to map
                for(moris::size_t i = 0; i< tNumNewPdvHpsts; i++)
                {
                    mNodeToPdvHostMap[aNodeIndices(0,i)] = i + tNumExistingPdvHosts;
                }
            }
            //------------------------------------------------------------------------------
            /*
             * returns the pdv host associated with the specified node index
             */
            GEN_Pdv_Host & get_pdv_host_from_manager( moris::moris_index const & aNodeIndex )
            {
                MORIS_ASSERT( mNodeToPdvHostMap.find(aNodeIndex) != mNodeToPdvHostMap.end(), "cl_GEN_Pdv_Host_Manager()::get_pdv_host_from_manager() - node index does not have an associated geometry object" );

                moris::moris_index tGOIndex = mNodeToPdvHostMap[aNodeIndex];

                return mPdvHosts(tGOIndex);
            }
            //------------------------------------------------------------------------------
            /*
             * returns the  pdv host associated with the specified node index
             * const version of above
             */
            GEN_Pdv_Host const & get_pdv_host_from_manager(moris::moris_index const & aNodeIndex) const
            {
              auto  tIter = mNodeToPdvHostMap.find(aNodeIndex);
              MORIS_ASSERT( tIter != mNodeToPdvHostMap.end(), "cl_GEN_Pdv_Host_Manager()::get_pdv_host_from_manager() - node index does not have an associated geometry object");

              moris::moris_index tGOIndex = tIter->second;

              return mPdvHosts(tGOIndex);
            }
            //------------------------------------------------------------------------------
             void link_to_node_to_another_nodes_pdv_host( moris::moris_index aNodeIndexWithPdvHost,
                                                          moris::moris_index aNodeIndexToLink )
             {
                 // pdv host index
                 moris::moris_index tPdvHostIndex = mNodeToPdvHostMap[ aNodeIndexWithPdvHost ];

                 // link new node by putting it's index in map with same tPdvHostIndex as the aNodeIndexWithPdvHost has
                 mNodeToPdvHostMap[aNodeIndexToLink] = tPdvHostIndex;
             }
             //------------------------------------------------------------------------------

        private:
            // pdv hosts
            moris::Cell< GEN_Pdv_Host > mPdvHosts;

            // Node to pdv host map (key - node index, val - geometry object index)
            std::unordered_map<moris::moris_index, moris::moris_index> mNodeToPdvHostMap;
        };
    }   // end ge namespace
}       // end ,moris namepspace



#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_ */
