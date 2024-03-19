//
// Created by frank on 1/11/24.
//

#include "cl_MTK_Ray_Line_Intersection.hpp"
#include "armadillo"

namespace moris::mtk
{
    /**
     * @brief Implementation according to https://stackoverflow.com/a/2932601
     */
    void Ray_Line_Intersection::perform_raytracing()
    {
        // Nomenclature:
        //     mOrigin: r
        //     mDirection: dr
        //        -> p(u) = r + u * dr
        //     mTargetOrigin: s
        //     mTargetSpan: ds = e - s (between line segment vertices s and e)
        //        -> q(v) = s + v * ds

        // distance between origins r and s
        // dOrigins = s - r
        Matrix< DDRMat > dOrigins = mTargetOrigin - mOrigin;

        // determinant of the matrix D = [dr, ds]
        MORIS_ASSERT( mSpatialDimension == 2, "Ray_Line_Intersection::perform: ray line intersection only implemented in 2D" );
        real const detD = mTargetSpan( 0 ) * mDirection( 1 ) - mTargetSpan( 1 ) * mDirection( 0 );

        // if detD is zero, the lines are parallel and no calculation is necessary
        if ( std::abs( detD ) > 1e-16 )
        {
            // scaling factor u and v for the ray and the line segment, respectively
            real u = NAN;
            real v = NAN;

            // u = (dOrigins_y * ds_x - dOrigins_x * ds_y) / detD
            u = ( dOrigins( 1 ) * mTargetSpan( 0 ) - dOrigins( 0 ) * mTargetSpan( 1 ) ) / detD;

            // v = (dOrigins_y * dr_x - dOrigins_x * dr_y) / detD
            v = ( dOrigins( 1 ) * mDirection( 0 ) - dOrigins( 0 ) * mDirection( 1 ) ) / detD;

            // if v is between 0 and 1, the intersection point is on the line segment
            if ( v >= 0.0 && v <= 1.0 )
            {
                mHasIntersection = true;

                // calculate the (signed) ray length between the origin and the line segment
                // negative if the intersection point is in the negative direction of the ray
                mSignedRayLength = norm( u * mDirection ) * ( u > 0.0 ? 1.0 : -1.0 );

                // set the factors by which the span of the line segment has to be multiplied to get the intersection point
                // q(v) = s + v * ds
                this->set_intersection_factors( { v } );

                // the parametric coordinate goes from -1 to 1 and has the center in the middle of the line segment
                mParametricCoordinate = 2.0 * ( v - 0.5 );
                return;    // early exit
            }
        }
        mHasIntersection = false;
    }

    Matrix< DDRMat > const Ray_Line_Intersection::get_intersection_parametric() const
    {
        return Matrix< DDRMat >( { { mParametricCoordinate } } );
    }

    real Ray_Line_Intersection::get_signed_ray_length() const
    {
        return mSignedRayLength;
    }
}    // namespace moris::mtk