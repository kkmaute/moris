/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_BSpline_Field.hpp
 *
 */

#pragma once

#include "cl_GEN_Field_Discrete_Integration.hpp"
#include "cl_MTK_Field.hpp"

namespace moris::ge
{
    class BSpline_Field : public Field_Discrete_Integration
    {

    private:
        uint mADVOffsetID;
        mtk::Mesh_Pair mMeshPair;
        uint mDiscretizationIndex;
        sol::Dist_Vector* mOwnedNodalValues = nullptr;
        sol::Dist_Vector* mSharedNodalValues = nullptr;

    public:
        /**
         * Constructor where ADVs are added based on an input field and a B-spline mesh.
         *
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         * @param aMeshPair The mesh pair where the discretization information can be obtained
         * @param aField Field for initializing the B-spline level set discretization
         */
        BSpline_Field(
                mtk::Mesh_Pair           aMeshPair,
                sol::Dist_Vector*        aOwnedADVs,
                const Matrix<DDSMat>&    aSharedADVIds,
                uint                     aADVOffsetID,
                uint                     aDiscretizationIndex,
                std::shared_ptr< Field > aField );

        BSpline_Field(
                sol::Dist_Vector*             aOwnedADVs,
                const Matrix<DDSMat>&         aSharedADVIds,
                uint                          aADVOffsetID,
                uint                          aDiscretizationIndex,
                std::shared_ptr< mtk::Field > aMTKField );

        /**
         * Destructor
         */
        ~BSpline_Field() override;

        /**
         * Given a node index, returns the field value.
         *
         * @param aNodeIndex Node index
         * @return Distance to this geometry
         */
        real get_field_value( uint aNodeIndex ) override;

        /**
         * Given a node index, evaluates the sensitivity of the geometry field with respect to all of the
         * geometry variables.
         *
         * @param aNodeIndex Node index
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs( uint aNodeIndex ) override;

        /**
         * Gets the IDs of ADVs which this field depends on for evaluations.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Matrix< DDSMat > get_determining_adv_ids( uint aNodeIndex ) override;

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector, and recomputes nodal values.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        void import_advs( sol::Dist_Vector* aOwnedADVs ) override;

        /**
         * Gets an MTK field, if this field needs to be remapped to a new mesh. In this implementation, B-spline coefficients
         * are set on the MTK field instead of just nodal values, as B-spline coefficients actually exist for a GEN B-spline field.
         *
         * @return MTK field
         */
        std::shared_ptr< mtk::Field > get_mtk_field() override;

    private:

        /**
         * Maps the level set field from nodes to B-splines for the given geometry.
         *
         * @return Target field
         */
        Matrix< DDRMat > map_to_bsplines( std::shared_ptr< Field > aField );

        void distribute_coeffs(
                const Matrix<DDRMat> & aTargetField,
                sol::Dist_Vector*      aOwnedADVs,
                uint                   aNumberOfCoefficients );

    };
}
