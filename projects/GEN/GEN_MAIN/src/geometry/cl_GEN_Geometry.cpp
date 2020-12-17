#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry(Intersection_Interpolation aIntersectionInterpolation)
                : mIntersectionInterpolation(aIntersectionInterpolation)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry::set_intersection_interpolation(std::string aInterpolationName)
        {
            if (aInterpolationName == "linear")
            {
                mIntersectionInterpolation = Intersection_Interpolation::LINEAR;
            }
            else if (aInterpolationName == "multilinear")
            {
                mIntersectionInterpolation = Intersection_Interpolation::MULTILINEAR;
            }
            else
            {
                MORIS_ERROR(false, aInterpolationName.append(" is not recognized as an intersection interpolation type.").c_str());
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Interpolation Geometry::get_intersection_interpolation()
        {
            return mIntersectionInterpolation;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
