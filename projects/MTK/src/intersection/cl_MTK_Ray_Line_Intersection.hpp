//
// Created by frank on 1/11/24.
//

#ifndef CL_MTK_RAY_LINE_INTERSECTION_HPP
#define CL_MTK_RAY_LINE_INTERSECTION_HPP

#include "cl_MTK_Ray_Intersection.hpp"
#include "cl_MTK_Space_Interpolator.hpp"

namespace moris::mtk
{
  class Space_Interpolator;

    class Ray_Line_Intersection : public Ray_Intersection
    {
      public:
        explicit Ray_Line_Intersection( uint const aSpatialDimension )
                : Ray_Intersection( Geometry_Type::LINE, aSpatialDimension ) {};

        ~Ray_Line_Intersection() override = default;

        void perform_raytracing() override;

        void perform_nonlinear_raytracing(
                Space_Interpolator&     tGeomSpaceInterpolator,
                Space_Interpolator&     tFieldSpaceInterpolator,
                Matrix< DDRMat >&       mOrigin,
                Matrix< DDRMat >&       mDirection,
                const Matrix< DDRMat >& tTargetLocalCoordinates );

        Matrix< DDRMat > get_intersection_parametric() const override;

        Matrix< DDRMat > get_ray_direction_param() const
        {
            return mRayDirectionParam;
        }

        real get_signed_ray_length() const override;

      private:
        Matrix< DDRMat >  mParametricCoordinate;
        real              mSignedRayLength = 0.0;
        Matrix< DDRMat >  mRayDirectionParam;
    };

}    // namespace moris::mtk

#endif    // CL_MTK_RAY_LINE_INTERSECTION_HPP
