/*
 * cl_GEN_Design_Variable_Interface.hpp
 *
 *  Created on: Jan 15, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_

#include "cl_GEN_Pdv_Host_Manager.hpp"

#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris
{
namespace ge
{
    class GEN_Design_Variable_Interface   :   MSI::Design_Variable_Interface
    {
    private:
        Pdv_Host_Manager* mManager;

    public:
//------------------------------------------------------------------------------
        GEN_Design_Variable_Interface( Pdv_Host_Manager* aManager ) : mManager(aManager)
        {  };
//------------------------------------------------------------------------------
        ~GEN_Design_Variable_Interface()
        {};
//------------------------------------------------------------------------------
        /**
         * @brief Function providing pdv values for requested vertex indices and Dv types
         *
         * @param[in] aIntegrationMeshSetIndex  Integration Mesh index
         * @param[in] aDvTypes                  List of Dv types
         *
         */
        void get_dv_types_for_set( const moris::moris_index            aIntegrationMeshSetIndex,
                                         Cell< Cell< enum GEN_DV > > & aDvTypes )
        {
            // this may change, hold off for now...
            MORIS_ASSERT( false,"GEN_Design_Variable_Interface::get_dv_types_for_set() - not implemented" );
        }
//------------------------------------------------------------------------------
        void get_unique_dv_types_for_set( const moris::moris_index           aIntegrationMeshSetIndex,
                                                       Cell< enum GEN_DV > & aDvTypes )
        {
            MORIS_ASSERT( false,"GEN_Design_Variable_Interface::get_unique_dv_types_for_set() - not implemented" );
        }
//------------------------------------------------------------------------------
        /**
         * @brief Function providing pdv values for requested vertex indices and Dv types
         *
         * @param[in] aNodeIndices   List of node indices
         * @param[in] aDvType        List of Dv types
         * @param[in] aDvValues      List of Dv values (DvType)(vertexIndex)
         * @param[in] aIsActive      List of active design variables (vertexIndex)(DvType)
         */

        void get_pdv_value( const Matrix< IndexMat >                & aNodeIndices,
                            const moris::Cell< enum GEN_DV >        & aDvTypes,
                            moris::Cell< moris::Matrix< DDRMat > >  & aDvValues,
                            moris::Cell< moris::Matrix< DDSMat > >  & aIsActiveDv )
        {
            uint tNumIndices = aNodeIndices.length();
            uint tNumTypes   = aDvTypes.size();

            aDvValues.resize( tNumTypes );
            aIsActiveDv.resize( tNumIndices );

            for( uint iInd=0; iInd<tNumIndices; iInd++ )
            {
                aIsActiveDv( iInd ).resize( tNumTypes, 1 );

                for( uint iType=0; iType<tNumTypes; iType++ )
                {
                    aDvValues( iType ).resize( tNumIndices, 1 );

                    aIsActiveDv( iInd )( iType ) = mManager->check_for_active_types( aNodeIndices(iInd), aDvTypes(iType) );
                    if( aIsActiveDv( iInd )( iType ) == 1 )
                    {
                        aDvValues( iType )( iInd ) = mManager->get_pdv_by_type_and_index( aNodeIndices(iInd), aDvTypes(iType) )->get_val()(0,0);
                    }
                }
            }
        }

        void get_pdv_value( const Matrix< IndexMat >                     & aNodeIndices,
                            const moris::Cell< enum GEN_DV >             & aDvTypes,
                                  moris::Cell< moris::Matrix< DDRMat > > & aDvValues )
        {
            MORIS_ASSERT(false, "Design_Variable_Interface::get_pdv_value() - overload not implemented on GE side");
        }
//------------------------------------------------------------------------------
        /**
         * @brief return local to global dv type map
         *
         */
        Matrix< DDSMat > get_my_local_global_map(  )
        {
            return mManager->get_global_map();
        }
//------------------------------------------------------------------------------
        /**
         * @brief return local to global DV type map
         *
         * @param[in] aVertexIndex   List of vertex indices
         * @param[in] aDvType        List of Dv types
         * @param[in] aDvIds         List of Dv Ids
         *
         */
        void get_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index > & aNodeIndices,
                                          const Cell< enum GEN_DV >               & aDvTypes,
                                                Cell< moris::Matrix< IdMat > >    & aDvIds )
        {
            /*
             * - each cell is a row vector of global IDs per each type
             * - return the global ids of the dv type on a specified vertex
             */

            uint tNumIndices = aNodeIndices.size();
            uint tNumTypes   = aDvTypes.size();

            moris::Cell< uint > tCounter( tNumTypes, 0 );

            aDvIds.resize( tNumTypes );
            for( uint Ik = 0; Ik <tNumTypes; Ik++)
            {
                aDvIds( Ik ).set_size( tNumIndices, 1);
            }

            for( uint iType=0; iType<tNumTypes; iType++ )
            {
                for( uint iInd=0; iInd<tNumIndices; iInd++ )
                {

                    bool tDvTypeExists = mManager->check_for_active_types( aNodeIndices(iInd), aDvTypes(iType) );   // flag for if the DV type exists on the current host

                    if( tDvTypeExists )
                    {
                        aDvIds( iType )( iInd ) = mManager->get_global_index_for_dv_type( aNodeIndices(iInd), aDvTypes(iType) );
                        tCounter( iType )++;
                    }
                }
            }

            for( uint Ik = 0; Ik <tNumTypes; Ik++)
            {
                aDvIds( Ik ).resize( tCounter( Ik ), 1);
            }
        }
//------------------------------------------------------------------------------

    };

}   // end ge namespace
}   // end moris namespace



#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_ */
