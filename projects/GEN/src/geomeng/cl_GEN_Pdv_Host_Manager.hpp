/*
 * cl_GEN_Pdv_Host_Manager.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_

//#include <unordered_map>

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
            /*
             * NOTE:
             * The user needs to be sure that the entry in aPdvList is associated with the same entry in aPropertyList.
             *
             * e.g.
             * aPdvList      = { GEN_PDV::DENSITY, GEN_PDV::MODULUS }
             * aPropertyList = { property 1,             property 2 }
             *
             * property 1 must correspond to DENSITY and property 2 must correspond to MODULUS
             */
            ~Pdv_Host_Manager(){};
            //------------------------------------------------------------------------------
            /*
             * Stores pdv host objects in the pdv host manager associated with nodes
             */
            void store_pdv_hosts( moris::Matrix< moris::IndexMat > const & aNodeIndices,
                                  moris::Cell< GEN_Pdv_Host >   const & aPdvHosts )
            {
                moris::size_t tNumExistingPdvHosts = mPdvHosts.size();

                moris::size_t tNumNewPdvHosts = aNodeIndices.n_cols();

                MORIS_ASSERT(tNumNewPdvHosts == aPdvHosts.size(),"cl_GEN_Pdv_Host_Manager()::store_pdv_hosts() - number of pdv hosts does not match number of node indices provided");

                // append the pdv object cell
                mPdvHosts.append( aPdvHosts );

                size_t tTotalPdvHosts = tNumExistingPdvHosts + tNumNewPdvHosts;
                mNodeToPdvHostMap.resize( tTotalPdvHosts, 2 );

                // Add to map
                for(moris::size_t i = 0; i< tNumNewPdvHosts; i++)
                {
                    mNodeToPdvHostMap( i, 0 ) = aNodeIndices(0,i);
                    mNodeToPdvHostMap( i, 1 ) = i + tNumExistingPdvHosts;
                }
            }

            //------------------------------------------------------------------------------
            /*
             * returns the pdv host associated with the specified node index
             */
            GEN_Pdv_Host & get_pdv_host_from_manager( moris::moris_index const & aNodeIndex )
            {
                MORIS_ASSERT( (moris_index)mNodeToPdvHostMap.n_rows() >= aNodeIndex, "cl_GEN_Pdv_Host_Manager()::get_pdv_host_from_manager() - node index does not have an associated pdv object" );

                moris_index tPDVIndex = mNodeToPdvHostMap( aNodeIndex, 1 );

                return mPdvHosts(tPDVIndex);
            }
            //------------------------------------------------------------------------------
            /*
             * returns the  pdv host associated with the specified node index
             * const version of above
             */
            GEN_Pdv_Host const & get_pdv_host_from_manager(moris::moris_index const & aNodeIndex) const
            {
                MORIS_ASSERT( (moris_index)mNodeToPdvHostMap.n_rows() >= aNodeIndex, "cl_GEN_Pdv_Host_Manager()::get_pdv_host_from_manager() - node index does not have an associated pdv object" );

                moris_index tPDVIndex = mNodeToPdvHostMap( aNodeIndex, 1 );

                return mPdvHosts(tPDVIndex);
            }
            //------------------------------------------------------------------------------
             void link_to_node_to_another_nodes_pdv_host( moris::moris_index aNodeIndexWithPdvHost,
                                                          moris::moris_index aNodeIndexToLink )
             {

                 mNodeToPdvHostMap( aNodeIndexWithPdvHost, 1 ) = aNodeIndexToLink;

//                 // pdv host index
//                 moris::moris_index tPdvHostIndex = mNodeToPdvHostMap[ aNodeIndexWithPdvHost ];
//                 // link new node by putting it's index in map with same tPdvHostIndex as the aNodeIndexWithPdvHost has
//                 mNodeToPdvHostMap[ aNodeIndexToLink ] = tPdvHostIndex;
             }
             //------------------------------------------------------------------------------
             uint get_num_pdv_hosts()
             {
                 return mPdvHosts.size();
             }
             //------------------------------------------------------------------------------
             /*
              * @brief returns a cell of the pdv types associated with the objects
              */
             Cell< enum GEN_PDV > get_pdv_types(  )
             {
                 return mPdvList;
             }
             //------------------------------------------------------------------------------
             void get_all_node_indices( Matrix< IndexMat > & aNodeIndices )
             {
                 uint tNumRows = mNodeToPdvHostMap.n_rows();

                 aNodeIndices.resize( tNumRows, 1 );

                 aNodeIndices.get_column( 0 ) = mNodeToPdvHostMap.get_column( 1 );
             }
             //------------------------------------------------------------------------------
             void create_association( Cell< enum GEN_PDV >  aPdvList,
                                      Cell< GEN_Property* > aPropertyList )
             {
                 mPdvList.append(aPdvList);
                 mProperties.append(aPropertyList);

                 uint tNumPdvs = aPdvList.size();
                 MORIS_ASSERT( tNumPdvs == aPropertyList.size(),"cl_GEN_Pdv_Host() - input pdv list does not match the size of the input property list" );

                 for( uint i=0; i<tNumPdvs; i++ )
                 {
                     mPdvToPropertyMap[ aPdvList(i) ] = i;
                 }

                 uint tNumHosts = mPdvHosts.size();
                 for(uint j=0; j<tNumHosts; j++)
                 {
                     this->get_pdv_host_from_manager( j ).set_property_list( mProperties );
                 }
             };
             //------------------------------------------------------------------------------
             void get_pdv_values( moris_index        aNodeIndex,
                                  const enum GEN_PDV aPdvType,
                                  Matrix< DDRMat > & aPdvValueMatrix )
             {
                 auto tSearch = mPdvToPropertyMap.find( aPdvType );

                 moris::moris_index tPropIndex = tSearch->second;
                 this->get_pdv_host_from_manager( aNodeIndex ).eval_pdv_property( tPropIndex, aPdvValueMatrix );
             }
             //------------------------------------------------------------------------------

        private:
            // pdv hosts
            moris::Cell< GEN_Pdv_Host > mPdvHosts;

            // node to host map ( column 0 - node index, column 1 - pdv object index )
            moris::Matrix< IndexMat > mNodeToPdvHostMap;
            //------------------------------------------------------------------------------
            Cell< enum GEN_PDV >  mPdvList;
            Cell< GEN_Property* > mProperties;

            // pdv to property map ( key - pdv enum, val - property index )
            std::unordered_map< enum GEN_PDV, moris::moris_index > mPdvToPropertyMap;
        };
    }   // end ge namespace
}       // end ,moris namepspace



#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_ */
