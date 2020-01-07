/*
 * cl_GEN_Pdv_Host.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_

#include <unordered_map>

#include "../additional/cl_GEN_Enums.hpp"
#include "../property/cl_GEN_Property.hpp"

// XTK includes
#include "cl_XTK_Topology.hpp"

namespace moris
{
    namespace ge
    {

        class GEN_Pdv_Host
        {
        private:
            //------------------------------------------------------------------------------
            moris::moris_index mParentEntityIndex;
            moris::moris_index mPhaseValIndex;

            // Parent topology
            std::shared_ptr< xtk::Topology > mParentTopology = nullptr;

            Cell< enum GEN_PDV > mPdvList;

            Cell< GEN_Property* > mProperties;

            // pdv to property map ( key - pdv enum, val - property index )
            std::unordered_map< enum GEN_PDV, moris::moris_index > mPdvToPropertyMap;

            //------------------------------------------------------------------------------
        public:
            //------------------------------------------------------------------------------
            GEN_Pdv_Host()
            {};

            /*
             * NOTE:
             * The user needs to be sure that the entry in aPdvList is associated with the same entry in aPropertyList.
             *
             * e.g.
             * aPdvList = { GEN_PDV::DENSITY, GEN_PDV::MODULUS }
             * aPropertyList = { property 1, property 2 }
             *
             * property 1 must correspond to DENSITY and property 2 must correspond to MODULUS
             */

            GEN_Pdv_Host( Cell< enum GEN_PDV >  aPdvList,
                          Cell< GEN_Property* > aPropertyList ):
                              mPdvList( aPdvList ),
                              mProperties( aPropertyList )
            {
                uint tNumPdvs = aPdvList.size();
                MORIS_ASSERT( tNumPdvs == aPropertyList.size(),"cl_GEN_Pdv_Host() - input pdv list does not match the size of the input property list" );

                for( uint i=0; i<tNumPdvs; i++ )
                {
                    mPdvToPropertyMap[ aPdvList(i) ] = i;
                }
            };
            //------------------------------------------------------------------------------
            virtual ~GEN_Pdv_Host(){};


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
            };

            //------------------------------------------------------------------------------

            void set_parent_entity_index( moris::moris_index aEntityIndex )
            {
                mParentEntityIndex = aEntityIndex;
            }

            moris::moris_index const & get_parent_entity_index( )
            {
                return mParentEntityIndex;
            }
            //------------------------------------------------------------------------------

            void set_parent_entity_topology( std::shared_ptr<xtk::Topology> tParentTopo )
            {
                MORIS_ASSERT(mParentTopology == nullptr, "cl_GEN_Pdv_Host::set_parent_entity_topology() - pdv host parent entity topology has already been set");
                mParentTopology = tParentTopo;
            }
            //------------------------------------------------------------------------------

            xtk::Topology const & get_parent_entity_topology( )
            {
                MORIS_ASSERT( mParentTopology != nullptr,
                        "cl_GEN_Pdv_Host::get_parent_entity_topology() - pdv host parent entity topology has not been set, either this is not an interface pdv host or set_parent_entity_topology was not called" );

                return (*mParentTopology);
            }
            //------------------------------------------------------------------------------
            void set_phase_val_row( moris::moris_index aPhaseValRowIndex )
            {
                mPhaseValIndex = aPhaseValRowIndex;
            }
            //------------------------------------------------------------------------------
            /*
             * @brief return the number of pdvs associated with the object
             */
            uint get_number_of_pdvs(  )
            {
                return mPdvList.size();
            }
            //------------------------------------------------------------------------------
            /*
             * @brief returns a cell of the pdv types associated with the object
             */
            Cell< enum GEN_PDV > get_pdv_types(  )
            {
                return mPdvList;
            }
            //------------------------------------------------------------------------------
            /*
             * @brief returns the value of the pdv defined by aPdvType
             */
            void get_pdv_values( const enum GEN_PDV aPdvType,
                                 Matrix< DDRMat > & aPdvValueMatrix )
            {
                auto tSearch = mPdvToPropertyMap.find( aPdvType );

                moris::moris_index tPDVIndex = tSearch->second;

                aPdvValueMatrix = mProperties( tPDVIndex )->val();
            }
            //------------------------------------------------------------------------------
        };

    }   // end ge namepsace
}       // end moris namespace





#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_ */
