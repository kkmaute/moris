//
// Created by frank on 1/11/24.
//

#ifndef CL_MTK_RAY_LINE_INTERSECTION_HPP
#define CL_MTK_RAY_LINE_INTERSECTION_HPP

#include "cl_MTK_Ray_Intersection.hpp"

namespace moris::mtk
{

    class Ray_Line_Intersection : public Ray_Intersection
    {
      public:
        explicit Ray_Line_Intersection( uint const aSpatialDimension )
                : Ray_Intersection( Geometry_Type::LINE, aSpatialDimension ){};

        ~Ray_Line_Intersection() override = default;

        void perform_raytracing() override;

        Matrix< DDRMat > const get_intersection_parametric() const override;

        real get_signed_ray_length() const override;

      private:
        real mParametricCoordinate = 0.0;
        real mSignedRayLength      = 0.0;
    };


}    // namespace moris::mtk

#endif    // CL_MTK_RAY_LINE_INTERSECTION_HPP
