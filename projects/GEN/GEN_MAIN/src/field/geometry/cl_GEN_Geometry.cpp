#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry(Geometry_Field_Parameters aParameters)
                : mParameters(aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry(std::shared_ptr<Geometry> aGeometry)
                : mParameters(aGeometry->mParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Geometry::set_intersection_interpolation(std::string aInterpolationName)
        {
            if (aInterpolationName == "linear")
            {
                mParameters.mIntersectionInterpolation = Intersection_Interpolation::LINEAR;
            }
            else if (aInterpolationName == "multilinear")
            {
                mParameters.mIntersectionInterpolation = Intersection_Interpolation::MULTILINEAR;
            }
            else
            {
                MORIS_ERROR(false, aInterpolationName.append(" is not recognized as an intersection interpolation type.").c_str());
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Interpolation Geometry::get_intersection_interpolation()
        {
            return mParameters.mIntersectionInterpolation;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
