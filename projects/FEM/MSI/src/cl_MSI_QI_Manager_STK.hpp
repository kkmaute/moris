/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_QI_Manager_STK.hpp
 *
 */

#pragma once

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris
{
    namespace MSI
    {
        class Equation_Model;
        class QI_Manager_STK : public Design_Variable_Interface
        {

          private:
            Matrix< DDSMat > mDummyMap;
            uint             mNumNodes = 0;
            Matrix< DDSMat > mPDVIds;
            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            QI_Manager_STK( std::shared_ptr< mtk::Integration_Mesh > aIgMesh, Vector< std::string >& aRequestedQIs );

            //------------------------------------------------------------------------------

            /**
             * get unique dv types for set
             * @param[ in ] aIntegrationMeshSetIndex
             * @param[ in ] aDvTypes
             */
            virtual void get_ip_unique_dv_types_for_set(
                    const moris_index             aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type >& aDvTypes ) const override;

            virtual void get_ig_unique_dv_types_for_set(
                    const moris_index             aIntegrationMeshSetIndex,
                    Vector< enum gen::PDV_Type >& aDvTypes ) const override;

            //------------------------------------------------------------------------------

            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ] aIntegrationMeshSetIndex  integration Mesh index
             * @param[ in ] aDvTypes                  list of group of dv types
             */
            virtual void get_ip_dv_types_for_set(
                    const moris_index                       aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) const override;

            virtual void get_ig_dv_types_for_set(
                    const moris_index                       aIntegrationMeshSetIndex,
                    Vector< Vector< enum gen::PDV_Type > >& aDvTypes ) const override;


            virtual void get_ig_pdv_value(
                    const Vector< moris_index >&        aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< Matrix< DDRMat > >&         aDvValues,
                    Vector< Vector< bool > >&           aIsActiveDv ) const override;

            //------------------------------------------------------------------------------

            /**
             * get pdv values for requested vertex indices and dv types
             * @param[ in ]     aNodeIndices list of vertex indices
             * @param[ in ]     aDvTypes     list of dv types
             * @param[ in/out ] aDvValues    list of dv values
             */
            virtual void get_ip_pdv_value(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< Matrix< DDRMat > >&         aDvValues ) const override;

            //------------------------------------------------------------------------------

            /**
             * return local to global dv map
             */
            virtual const Matrix< DDSMat >& get_my_local_global_map() override;
            // BRENDAN to make this work in parallel

            //------------------------------------------------------------------------------

            /**
             * return local to global dv type map
             * @param[ in ] aVertexIndex   List of vertex indices
             * @param[ in ] aDvType        List of Dv types
             * @param[ in ] aDvIds         List of Dv Ids
             */
            virtual void get_ip_dv_ids_for_type_and_ind(
                    const Matrix< IndexMat >&           aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< Matrix< IdMat > >&          aDvIds ) const override;

            virtual void get_ig_dv_ids_for_type_and_ind(
                    const Vector< moris_index >&        aNodeIndices,
                    const Vector< enum gen::PDV_Type >& aDvTypes,
                    Vector< Vector< moris_index > >&    aDvIds ) const override;

            //------------------------------------------------------------------------------

            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ip_requested_dv_types( Vector< enum gen::PDV_Type >& aDvTypes ) const override;

            //------------------------------------------------------------------------------

            /**
             * get requested dv types for sensitivity analysis
             * @param[ in ] aDvTypes list of dv types to fill
             */
            virtual void get_ig_requested_dv_types( Vector< enum gen::PDV_Type >& aDvTypes ) override;

            //------------------------------------------------------------------------------

            /**
             * Tells FEM to ensure we create PDV IDs for GEN nodes for optimization
             * FIXME @bc: There should be a cleaner way to handle this through inheritance - the FEM model, relies on knowing this info
             * to avoid creating PDV IDs for nodes if we are doing STK. Thus if body fitted shape optimization was implemented this would
             * be unnecessary.
             */
            virtual bool is_gen_workflow() const final;

        };    // class QI_Manager_STK

        //------------------------------------------------------------------------------

    }    // namespace MSI
}    // namespace moris