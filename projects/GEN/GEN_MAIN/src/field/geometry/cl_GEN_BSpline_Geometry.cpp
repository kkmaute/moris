#include "cl_GEN_BSpline_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Geometry::BSpline_Geometry(
                sol::Dist_Vector*         aOwnedADVs,
                const Matrix<DDUMat>&     aCoefficientIndices,
                const Matrix<DDSMat>&     aSharedADVIds,
                uint                      aADVOffsetID,
                mtk::Mesh_Pair            aMeshPair,
                std::shared_ptr<Geometry> aGeometry)
                : Field(aCoefficientIndices, aSharedADVIds, aMeshPair, aGeometry)
                , BSpline_Field(aOwnedADVs, aCoefficientIndices, aSharedADVIds, aADVOffsetID, aMeshPair, aGeometry)
                , Geometry(aGeometry)
        {

        }

        BSpline_Geometry::BSpline_Geometry(
                sol::Dist_Vector*         aOwnedADVs,
                const Matrix<DDUMat>&     aCoefficientIndices,
                const Matrix<DDSMat>&     aSharedADVIds,
                uint                      aADVOffsetID,
                mtk::Mesh_Pair            aMeshPair,
                std::shared_ptr<Geometry> aGeometry,
                std::shared_ptr<mtk::Field> aField )
                : Field(aCoefficientIndices, aSharedADVIds, aMeshPair, aGeometry)
                , BSpline_Field(aOwnedADVs, aCoefficientIndices, aSharedADVIds, aADVOffsetID, aMeshPair, aGeometry)
                , Geometry(aGeometry)
        {

        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
