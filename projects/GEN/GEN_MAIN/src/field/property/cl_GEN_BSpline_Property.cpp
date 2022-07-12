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
                mtk::Mesh_Pair            aMeshPair,
                std::shared_ptr<Property> aProperty)
                : Field(aCoefficientIndices, aSharedADVIds, aMeshPair, aProperty)
                , BSpline_Field(aOwnedADVs, aCoefficientIndices, aSharedADVIds, aADVOffsetID, aMeshPair, aProperty)
                , Property(aProperty)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Property::BSpline_Property(
                sol::Dist_Vector*           aOwnedADVs,
                const Matrix< DDUMat >&     aCoefficientIndices,
                const Matrix< DDSMat >&     aSharedADVIds,
                uint                        aADVOffsetID,
                mtk::Mesh_Pair              aMeshPair,
                std::shared_ptr< Property > aProperty,
                 std::shared_ptr<mtk::Field> aField )
                : Field( aCoefficientIndices, aSharedADVIds, aMeshPair, aProperty )
                , BSpline_Field( aOwnedADVs, aCoefficientIndices, aSharedADVIds, aADVOffsetID, aMeshPair, aProperty,aField )
                , Property( aProperty )
        {
        }
    }
}
