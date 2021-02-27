#include "cl_GEN_BSpline_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Property::BSpline_Property(
                sol::Dist_Vector*         aOwnedADVs,
                const Matrix<DDUMat>&     aCoefficientIndices,
                const Matrix<DDSMat>&     aOwnedADVIds,
                const Matrix<DDSMat>&     aSharedADVIds,
                uint                      aOwnedADVIdsOffset,
                mtk::Interpolation_Mesh*  aMesh,
                std::shared_ptr<Property> aProperty)
                : Field(aCoefficientIndices, aSharedADVIds, aProperty)
                , BSpline_Field(aOwnedADVs, aCoefficientIndices, aOwnedADVIds, aSharedADVIds, aOwnedADVIdsOffset, aMesh, aProperty)
                , Property(aProperty)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
