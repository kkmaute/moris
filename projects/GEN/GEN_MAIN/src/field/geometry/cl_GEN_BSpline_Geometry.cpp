#include "cl_GEN_BSpline_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Geometry::BSpline_Geometry(
                sol::Dist_Vector*         aOwnedADVs,
                const Matrix<DDUMat>&     aCoefficientIndices,
                const Matrix<DDSMat>&     aOwnedADVIds,
                const Matrix<DDSMat>&     aSharedADVIds,
                uint                      aOwnedADVIdsOffset,
                mtk::Interpolation_Mesh*  aMesh,
                std::shared_ptr<Geometry> aGeometry)
                : Field(aCoefficientIndices, aSharedADVIds, aGeometry)
                , BSpline_Field(aOwnedADVs, aCoefficientIndices, aOwnedADVIds, aSharedADVIds, aOwnedADVIdsOffset, aMesh, aGeometry)
                , Geometry(aGeometry)
        {

        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
