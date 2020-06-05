#include "cl_GEN_Pdv_Intersection.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Intersection::Pdv_Intersection(GEN_Geometry_Object* aIntersection, uint aDimension)
        : mIntersection(aIntersection),
          mDimension(aDimension)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Pdv_Intersection::get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return aCoordinates(mDimension);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Intersection::get_sensitivity(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
            aSensitivities = mIntersection->get_sensitivity_dx_dp().get_row(mDimension);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
