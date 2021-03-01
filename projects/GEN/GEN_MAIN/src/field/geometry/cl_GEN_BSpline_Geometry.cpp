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
                mtk::Interpolation_Mesh*  aMesh,
                std::shared_ptr<Geometry> aGeometry)
                : Field(aCoefficientIndices, aSharedADVIds, aGeometry)
                , BSpline_Field(aOwnedADVs, aCoefficientIndices, aSharedADVIds, aADVOffsetID, aMesh, aGeometry)
                , Geometry(aGeometry)
        {

        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
