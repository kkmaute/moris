#include "cl_GEN_BSpline_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Property::BSpline_Property(
                sol::Dist_Vector*         aOwnedADVs,
                const Matrix<DDUMat>&     aCoefficientIndices,
                const Matrix<DDSMat>&     aSharedADVIds,
                uint                      aADVOffsetID,
                mtk::Interpolation_Mesh*  aMesh,
                std::shared_ptr<Property> aProperty)
                : Field(aCoefficientIndices, aSharedADVIds, aProperty)
                , BSpline_Field(aOwnedADVs, aCoefficientIndices, aSharedADVIds, aADVOffsetID, aMesh, aProperty)
                , Property(aProperty)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
