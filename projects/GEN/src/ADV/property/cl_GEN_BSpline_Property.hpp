/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_BSpline_Property.hpp
 *
 */

//#pragma once
//
//#include "cl_GEN_BSpline_Field.hpp"
//#include "cl_GEN_Property.hpp"
//
//namespace moris::ge
//{
//    class BSpline_Property : public BSpline_Field
//    {
//    public:
//        /**
//         * Constructor where ADVs are added based on an input field and a B-spline mesh.
//         *
//         * @param aOwnedADVs Pointer to the owned distributed ADVs
//         * @param aCoefficientIndices Coefficient indices to be mapped to
//         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
//         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
//         * @param aMeshPair The mesh pair where the discretization information can be obtained
//         * @param aProperty Property for initializing the B-spline level set discretization
//         */
//        BSpline_Property(
//                sol::Dist_Vector*         aOwnedADVs,
//                const Matrix<DDUMat>&     aCoefficientIndices,
//                const Matrix<DDSMat>&     aSharedADVIds,
//                uint                      aADVOffsetID,
//                mtk::Mesh_Pair            aMeshPair,
//                std::shared_ptr<Property> aProperty);
//
//        //--------------------------------------------------------------------------------------
//        /**
//         * Constructor where ADVs are added based on an input field and a B-spline mesh.
//         *
//         * @param aOwnedADVs Pointer to the owned distributed ADVs
//         * @param aCoefficientIndices Coefficient indices to be mapped to
//         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
//         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
//         * @param aMeshPair The mesh pair where the discretization information can be obtained
//         * @param aProperty Property for initializing the B-spline level set discretization
//         * @param aFiled Field that will be used to obtain the coefficients directly
//         */
//        BSpline_Property(
//                sol::Dist_Vector*           aOwnedADVs,
//                const Matrix< DDUMat >&     aCoefficientIndices,
//                const Matrix< DDSMat >&     aSharedADVIds,
//                uint                        aADVOffsetID,
//                mtk::Mesh_Pair              aMeshPair,
//                std::shared_ptr< Property > aProperty,
//                std::shared_ptr<mtk::Field> aField );
//    };
//}
