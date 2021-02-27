#ifndef MORIS_CL_GEN_BSPLINE_PROPERTY_HPP
#define MORIS_CL_GEN_BSPLINE_PROPERTY_HPP

#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {
        class BSpline_Property : public BSpline_Field, public Property
        {
        public:
            /**
             * Constructor where ADVs are added based on an input field and a B-spline mesh.
             *
             * @param aOwnedADVs Pointer to the owned distributed ADVs
             * @param aOwnedADVIds All owned ADV IDs on this processor
             * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
             * @param aOwnedADVIdsOffset Offset in the owned ADV IDs for pulling ADV IDs
             * @param aMesh The mesh pointer where the B-spline information can be obtained
             * @param aProperty Property for initializing the B-spline level set discretization
             */
            BSpline_Property(
                    sol::Dist_Vector*         aOwnedADVs,
                    const Matrix<DDUMat>&     aCoefficientIndices,
                    const Matrix<DDSMat>&     aOwnedADVIds,
                    const Matrix<DDSMat>&     aSharedADVIds,
                    uint                      aOwnedADVIdsOffset,
                    mtk::Interpolation_Mesh*  aMesh,
                    std::shared_ptr<Property> aProperty);
        };
    }
}

#endif //MORIS_CL_GEN_BSPLINE_PROPERTY_HPP
