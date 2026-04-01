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
        return mDummyMap;
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
        MORIS_ERROR( false, "QI_Manager_STK - function not implemented." );
        return;
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