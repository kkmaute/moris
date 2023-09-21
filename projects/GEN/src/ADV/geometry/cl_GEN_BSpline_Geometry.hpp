/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_BSpline_Geometry.hpp
 *
 */

//#ifndef MORIS_CL_GEN_BSPLINE_GEOMETRY_HPP
//#define MORIS_CL_GEN_BSPLINE_GEOMETRY_HPP
//
//#include "cl_GEN_BSpline_Field.hpp"
//#include "cl_GEN_Level_Set_Geometry.hpp"
//
//namespace moris
//{
//    namespace ge
//    {
//        class BSpline_Geometry : public BSpline_Field, public Level_Set_Geometry
//        {
//        public:
//
//            /**
//             * Constructor where ADVs are added based on an input field and a B-spline mesh.
//             *
//             * @param aOwnedADVs Pointer to the owned distributed ADVs
//             * @param aCoefficientIndices Coefficient indices to be mapped to
//             * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
//             * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
//             * @param aMeshPair The mesh pair where the discretization information can be obtained
//             * @param aGeometry Geometry for initializing the B-spline level set discretization
//             */
//            BSpline_Geometry(
//                    sol::Dist_Vector*         aOwnedADVs,
//                    const Matrix<DDUMat>&     aCoefficientIndices,
//                    const Matrix<DDSMat>&     aSharedADVIds,
//                    uint                      aADVOffsetID,
//                    mtk::Mesh_Pair            aMeshPair,
//                    std::shared_ptr< Level_Set_Geometry > aGeometry);
//
//            BSpline_Geometry(
//                    sol::Dist_Vector*         aOwnedADVs,
//                    const Matrix<DDUMat>&     aCoefficientIndices,
//                    const Matrix<DDSMat>&     aSharedADVIds,
//                    uint                      aADVOffsetID,
//                    mtk::Mesh_Pair            aMeshPair,
//                    std::shared_ptr< Level_Set_Geometry > aGeometry,
//                    std::shared_ptr<mtk::Field> aField );
//        };
//    }
//}
//
//#endif //MORIS_CL_GEN_BSPLINE_GEOMETRY_HPP

