/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_PDV_Host_Manager.cpp
 *
 */

#include "cl_MSI_QI_Manager_STK.hpp"

// detailed logging
#include "cl_Tracer.hpp"

namespace moris::MSI
{
    //--------------------------------------------------------------------------------------------------------------

    QI_Manager_STK::QI_Manager_STK( std::shared_ptr< mtk::Integration_Mesh > aIgMesh, Vector< std::string >& aRequestedQIs )
            : Design_Variable_Interface( aRequestedQIs )
    {
        // Store number of nodes and pdv ids for use in local to global map
        mNumNodes = aIgMesh->get_num_entities( mtk::EntityRank::NODE );

        // Assign PDV IDs to all spatial directions and all nodes
        mPDVIds.set_size( mNumNodes * aIgMesh->get_spatial_dim(), 1 );
        uint tDim = aIgMesh->get_spatial_dim();
        for ( uint iNode = 0; iNode < mNumNodes; iNode++ )
        {
            for ( uint iDim = 0; iDim < tDim; iDim++ )
            {
                mPDVIds( iNode * tDim + iDim ) = iNode * tDim + iDim;
            }
        }
    }

    /**
     * get unique dv types for set
     * @param[ in ] aIntegrationMeshSetIndex
     * @param[ in ] aDvTypes
     */
    void QI_Manager_STK::get_ip_unique_dv_types_for_set(
            const moris_index             aIntegrationMeshSetIndex,
            Vector< enum gen::PDV_Type >& aDvTypes ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    void QI_Manager_STK::get_ig_unique_dv_types_for_set(
            const moris_index             aIntegrationMeshSetIndex,
            Vector< enum gen::PDV_Type >& aDvTypes ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    //------------------------------------------------------------------------------

    void QI_Manager_STK::get_ip_dv_types_for_set(
            const moris_index                       aIntegrationMeshSetIndex,
            Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    //------------------------------------------------------------------------------

    void QI_Manager_STK::get_ig_dv_types_for_set(
            const moris_index                       aIntegrationMeshSetIndex,
            Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    //------------------------------------------------------------------------------

    void QI_Manager_STK::get_ig_pdv_value(
            const Vector< moris_index >&        aNodeIndices,
            const Vector< enum gen::PDV_Type >& aDvTypes,
            Vector< Matrix< DDRMat > >&         aDvValues,
            Vector< Vector< bool > >&           aIsActiveDv ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    //------------------------------------------------------------------------------

    void QI_Manager_STK::get_ip_pdv_value(
            const Matrix< IndexMat >&           aNodeIndices,
            const Vector< enum gen::PDV_Type >& aDvTypes,
            Vector< Matrix< DDRMat > >&         aDvValues ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDSMat >& QI_Manager_STK::get_my_local_global_map()
    {
        // FIXME BRENDAN: THIS WILL NOT WORK IN PARALLEL AND ASSUMES THAT ALL NODES HAVE THE SAME NUMBER OF PDVS IN THE SAME ORDER.
        //  THIS FUNCTION WAS NOT INTENDED TO BE USED ROBUSTLY AND SHOULD BE REWORKED IF IT IS TO BE USED IN ANY REAL CAPACITY
        return mPDVIds;
    }

    //------------------------------------------------------------------------------

    void QI_Manager_STK::get_ip_dv_ids_for_type_and_ind(
            const Matrix< IndexMat >&           aNodeIndices,
            const Vector< enum gen::PDV_Type >& aDvTypes,
            Vector< Matrix< IdMat > >&          aDvIds ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    void QI_Manager_STK::get_ig_dv_ids_for_type_and_ind(
            const Vector< moris_index >&        aNodeIndices,
            const Vector< enum gen::PDV_Type >& aDvTypes,
            Vector< Vector< moris_index > >&    aDvIds ) const
    {
        // FIXME BRENDAN THIS FUNCTION IS NOT CONSISTENT WITH THE STORED PDV IDS CREATED BY THE CONSTRUCTOR
        // THIS ASSUMES THAT ALL PDV TYPES ARE REQUESTED FOR ALL NODES, WHICH MAY NOT BE THE CASE.
        // THIS FUNCTION ALSO ASSUMES A PARTICULAR ORDERING OF THE PDV IDS (ALL TYPES FOR NODE 1, THEN ALL TYPES FOR NODE 2, ETC) WHICH MAY NOT BE THE CASE
        // I CANNOT EMPHASIZE ENOUGH HOW MUCH THIS FUNCTION WAS NOT INTENDED TO BE USED ROBUSTLY

        // get the number of node indices requested
        uint tNumIndices = aNodeIndices.size();

        // get the number of dv types requested
        uint tNumTypes = aDvTypes.size();

        // set size for list of dv values
        aDvIds.resize( tNumTypes );

        // loop over the requested dv types
        for ( uint tPDVTypeIndex = 0; tPDVTypeIndex < tNumTypes; tPDVTypeIndex++ )
        {
            aDvIds( tPDVTypeIndex ).resize( tNumIndices, -1 );

            // loop over the node indices
            for ( uint iN = 0; iN < tNumIndices; iN++ )
            {
                aDvIds( tPDVTypeIndex )( iN ) = iN * tNumTypes + static_cast< uint >( aDvTypes( tPDVTypeIndex ) );
            }
        }
    }

    //------------------------------------------------------------------------------

    void QI_Manager_STK::get_ip_requested_dv_types( Vector< enum gen::PDV_Type >& aDvTypes ) const
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    //------------------------------------------------------------------------------

    void QI_Manager_STK::get_ig_requested_dv_types( Vector< enum gen::PDV_Type >& aDvTypes )
    {
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
    }

    //------------------------------------------------------------------------------

    bool QI_Manager_STK::is_gen_workflow() const
    {
        return false;
    }
}    // namespace moris::MSI