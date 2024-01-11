//
// Created by frank on 1/11/24.
//

#ifndef CL_MTK_RAY_INTERSECTION_HPP
#define CL_MTK_RAY_INTERSECTION_HPP

#include "cl_Vector.hpp"
#include <cl_Matrix.hpp>
#include <fn_norm.hpp>
#include <linalg_typedefs.hpp>
#include <moris_typedefs.hpp>
#include <cl_MTK_Enums.hpp>

namespace moris::mtk
{
    class Ray_Intersection
    {
      public:
        Ray_Intersection(
                Geometry_Type const aGeometryType,
                uint const          aSpatialDimension )
                : mGeometryType( aGeometryType )
                , mSpatialDimension( aSpatialDimension )
                , mOrigin( Matrix< DDRMat >( mSpatialDimension, 1, 0.0 ) )
                , mDirection( Matrix< DDRMat >( mSpatialDimension, 1, 0.0 ) )
                , mTargetOrigin( Matrix< DDRMat >( mSpatialDimension, 1, 0.0 ) )
                , mTargetSpan( Matrix< DDRMat >( mSpatialDimension, get_target_span_dimension( aGeometryType ), 0.0 ) )
        {
            MORIS_ASSERT( aSpatialDimension == 2 || aSpatialDimension == 3, "Ray_Intersection::Ray_Intersection: spatial dimension must be 2 or 3" );
        };

        virtual ~Ray_Intersection() = default;

        void set_ray_origin( const Matrix< DDRMat > &aOrigin )
        {
            MORIS_ASSERT( aOrigin.n_cols() == 1, "Ray_Intersection::set_origin: origin must be a column vector" );
            MORIS_ASSERT( aOrigin.n_rows() == mSpatialDimension, "Ray_Intersection::set_origin: origin must be a 3x1 vector" );
            mOrigin = aOrigin;
        }

        void set_ray_direction( const Matrix< DDRMat > &aDirection )
        {
            MORIS_ASSERT( aDirection.n_cols() == 1, "Ray_Intersection::set_direction: direction must be a column vector" );
            MORIS_ASSERT( aDirection.n_rows() == mSpatialDimension, "Ray_Intersection::set_direction: direction must be a 3x1 vector" );
            mDirection = aDirection;
        }

        void set_target_origin( const Matrix< DDRMat > &aTargetOrigin )
        {
            MORIS_ASSERT( aTargetOrigin.n_cols() == 1, "Ray_Intersection::set_target_origin: target origin must be a column vector" );
            MORIS_ASSERT( aTargetOrigin.n_rows() == mSpatialDimension, "Ray_Intersection::set_target_origin: target origin must be a 3x1 vector" );
            mTargetOrigin = aTargetOrigin;
        }

        void set_target_span( const Matrix< DDRMat > &aTargetSpan )
        {
            MORIS_ASSERT( aTargetSpan.n_cols() == mTargetSpan.n_cols(), "Ray_Intersection::set_target_span: target span has not the correct number of columns" );
            MORIS_ASSERT( aTargetSpan.n_rows() == mSpatialDimension, "Ray_Intersection::set_target_span: target span must be a 3xn matrix" );
            mTargetSpan = aTargetSpan;
        }

        virtual void perform_raytracing() = 0;

        bool has_intersection() const
        {
            return mHasIntersection;
        }

        virtual Matrix< DDRMat > const get_intersection_parametric() const = 0;

        Matrix< DDRMat > get_intersection_physical() const
        {
            MORIS_ASSERT( mHasIntersection, "Ray_Intersection::get_intersection_physical: no intersection" );
            Matrix< DDRMat > tIntersection = mTargetOrigin;
            for ( uint i = 0; i < mTargetSpan.n_cols(); ++i )
            {
                tIntersection += mIntersectionFactors( i ) * mTargetSpan.get_column( i );
            }
            return tIntersection;
        }

        virtual real get_ray_length() const
        {
            return norm( get_intersection_physical() - mOrigin );
        };

      protected:
        void set_intersection_factors( Vector< real > const &aIntersectionFactors )
        {
            MORIS_ASSERT( aIntersectionFactors.size() == mTargetSpan.n_cols(), "Intersection factors should have the same number of elements as the target span vectors" );
            mIntersectionFactors = aIntersectionFactors;
        }

        static int get_target_span_dimension( Geometry_Type const &aGeometryType )
        {
            switch ( aGeometryType )
            {
                case Geometry_Type::LINE:
                    return 1;
                case Geometry_Type::TRI:
                    return 2;
                case Geometry_Type::QUAD:
                    return 2;
                default:
                    MORIS_ERROR( false, "Ray_Intersection::get_target_span_dimension: invalid geometry type" );
                    return 0;
            }
        }

        Geometry_Type    mGeometryType;
        uint             mSpatialDimension;
        bool             mHasIntersection = false;
        Matrix< DDRMat > mOrigin;                 // origin of the ray
        Matrix< DDRMat > mDirection;              // direction of the ray
        Matrix< DDRMat > mTargetOrigin;           // reference point on the target geometry (e.g. corner of a tri cell)
        Matrix< DDRMat > mTargetSpan;             // vectors that span the target geometry (e.g. two edges of a tri cell)
        Vector< real >   mIntersectionFactors;    // factors to multiply the target span vectors with to get the intersection point
    };
}    // namespace moris::mtk


#endif    // CL_MTK_RAY_INTERSECTION_HPP
