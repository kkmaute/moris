#include "cl_MTK_Interpolation_Rule.hpp"             //MTK/src
#include "cl_MTK_Interpolation_Function_Base.hpp"    //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        Interpolation_Rule::Interpolation_Rule( const Geometry_Type       & aGeometryType,
                                                const Interpolation_Type       & aSpaceInterpolationType,
                                                const Interpolation_Order & aSpaceInterpolationOrder,
                                                const Interpolation_Type       & aTimeInterpolationType,
                                                const Interpolation_Order & aTimeInterpolationOrder)
                                              : mGeometryType( aGeometryType ),
                                                mSpaceInterpolationType( aSpaceInterpolationType ),
                                                mSpaceInterpolationOrder( aSpaceInterpolationOrder ),
                                                mTimeGeometryType( Geometry_Type::LINE ),
                                                mTimeInterpolationType( aTimeInterpolationType ),
                                                mTimeInterpolationOrder( aTimeInterpolationOrder )
        {}

        Interpolation_Rule::Interpolation_Rule( const Geometry_Type       & aGeometryType,
                                                const Interpolation_Type       & aSpaceInterpolationType,
                                                const Interpolation_Order & aSpaceInterpolationOrder,
                                                const Geometry_Type       & aTimeGeometryType,
                                                const Interpolation_Type       & aTimeInterpolationType,
                                                const Interpolation_Order & aTimeInterpolationOrder)
                                              : mGeometryType( aGeometryType ),
                                                mSpaceInterpolationType( aSpaceInterpolationType ),
                                                mSpaceInterpolationOrder( aSpaceInterpolationOrder ),
                                                mTimeGeometryType( aTimeGeometryType ),
                                                mTimeInterpolationType( aTimeInterpolationType ),
                                                mTimeInterpolationOrder( aTimeInterpolationOrder )
        {}

//------------------------------------------------------------------------------

        Interpolation_Function_Base *
        Interpolation_Rule::create_space_interpolation_function() const
        {
            // create the factory
            Interpolation_Function_Factory tFactory;

            // return new interpolation function pointer
            return tFactory.create_interpolation_function( mGeometryType,
                                                           mSpaceInterpolationType,
                                                           mSpaceInterpolationOrder );
        }

//------------------------------------------------------------------------------

        Interpolation_Function_Base *
        Interpolation_Rule::create_time_interpolation_function() const
        {
            // create the factory
            Interpolation_Function_Factory tFactory;

            // return new interpolation function pointer
            return tFactory.create_interpolation_function( mTimeGeometryType,
                                                           mTimeInterpolationType,
                                                           mTimeInterpolationOrder );
        }

//----------------------------------------------------------------------------
     } /* namespace mtk */
} /* namespace moris */