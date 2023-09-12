/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry.cpp
 *
 */

#include "cl_GEN_Geometry.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry( Geometry_Field_Parameters aParameters )
                : mParameters( aParameters )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Geometry::Geometry( std::shared_ptr< Geometry > aGeometry )
                : mParameters( aGeometry->mParameters )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry::set_intersection_interpolation( std::string aInterpolationName )
        {
            if ( aInterpolationName == "linear" )
            {
                mParameters.mIntersectionInterpolation = Intersection_Interpolation::LINEAR;
            }
            else if ( aInterpolationName == "multilinear" )
            {
                mParameters.mIntersectionInterpolation = Intersection_Interpolation::MULTILINEAR;
            }
            else
            {
                MORIS_ERROR( false, "%s is not recognized as an intersection interpolation type.", aInterpolationName.c_str() );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Interpolation
        Geometry::get_intersection_interpolation()
        {
            return mParameters.mIntersectionInterpolation;
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void
        Geometry::set_intersection_mode( Intersection_Mode aIntersectionMode )
        {
            mParameters.mIntersectionMode = aIntersectionMode;
        }

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Mode
        Geometry::get_intersection_mode()
        {
            return mParameters.mIntersectionMode;
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void
        Geometry::set_isocontour_threshold( real aIsocontourThreshold )
        {
            mParameters.mIsocontourThreshold = aIsocontourThreshold;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Geometry::get_isocontour_threshold()
        {
            return mParameters.mIsocontourThreshold;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry::set_isocontour_tolerance( real aIsocontourTolerance )
        {
            MORIS_ERROR( aIsocontourTolerance > 1e-14,
                "Geometry_Engine::Geometry_Engine - Isocontour tolerance should be larger than 1e-14" );
            mParameters.mIsocontourTolerance = aIsocontourTolerance;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Geometry::get_isocontour_tolerance()
        {
            return mParameters.mIsocontourTolerance;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Geometry::set_intersection_tolerance( real aIntersectionTolerance )
        {
            MORIS_ERROR( aIntersectionTolerance > 1e-14,
                "Geometry_Engine::Geometry_Engine - Intersection tolerance should be larger than 1e-14" );
            mParameters.mIntersectionTolerance = aIntersectionTolerance;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Geometry::get_intersection_tolerance()
        {
            return mParameters.mIntersectionTolerance;
        }

        //--------------------------------------------------------------------------------------------------------------

        
    }    // namespace ge
}    // namespace moris
