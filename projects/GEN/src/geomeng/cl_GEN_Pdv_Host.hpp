/*
 * cl_GEN_Pdv_Host.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_

#include "../additional/cl_GEN_Enums.hpp"
#include "../property/cl_GEN_Property.hpp"

// XTK includes
#include "cl_XTK_Topology.hpp"

namespace moris
{
    namespace ge
    {
//        struct Pdv
//        {
//
//        };

        class GEN_Pdv_Host
        {
            /*
             * host has list of property pointers
             * manager knows the association
             * manager tells host which property in the list to evaluate
             */

        private:
            //------------------------------------------------------------------------------
            moris::moris_index mParentEntityIndex;
            moris::moris_index mPhaseValIndex;

            // Parent topology
            std::shared_ptr< xtk::Topology > mParentTopology = nullptr;

            //------------------------------------------------------------------------------
            Cell< GEN_Property* > mProperties;

//            Cell< >

            //------------------------------------------------------------------------------
        public:
            //------------------------------------------------------------------------------
            GEN_Pdv_Host()
            {};

            //------------------------------------------------------------------------------
            ~GEN_Pdv_Host(){};

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
             * @brief set the list of property pointers
             */
            void set_property_list( Cell< GEN_Property* > & aPropList )
            {
                mProperties = aPropList;
            }
            //------------------------------------------------------------------------------
            /*
             * @brief returns the value of the pdv defined by aPdvType
             */
            void eval_pdv_property( const moris_index  aPropIndex,
                                    Matrix< DDRMat > & aPdvValueMatrix )
            {
                aPdvValueMatrix = mProperties( aPropIndex )->val();
            }
            //------------------------------------------------------------------------------
        };

    }   // end ge namepsace
}       // end moris namespace





#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_HPP_ */
