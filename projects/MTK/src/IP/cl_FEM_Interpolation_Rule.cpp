#include "cl_FEM_Interpolation_Rule.hpp"             //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp"    //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Factory.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Interpolation_Rule::Interpolation_Rule( const mtk::Geometry_Type       & aGeometryType,
                                                const Interpolation_Type       & aSpaceInterpolationType,
                                                const mtk::Interpolation_Order & aSpaceInterpolationOrder,
                                                const Interpolation_Type       & aTimeInterpolationType,
                                                const mtk::Interpolation_Order & aTimeInterpolationOrder)
                                              : mGeometryType( aGeometryType ),
                                                mSpaceInterpolationType( aSpaceInterpolationType ),
                                                mSpaceInterpolationOrder( aSpaceInterpolationOrder ),
                                                mTimeGeometryType( mtk::Geometry_Type::LINE ),
                                                mTimeInterpolationType( aTimeInterpolationType ),
                                                mTimeInterpolationOrder( aTimeInterpolationOrder )
        {}

        Interpolation_Rule::Interpolation_Rule( const mtk::Geometry_Type       & aGeometryType,
                                                const Interpolation_Type       & aSpaceInterpolationType,
                                                const mtk::Interpolation_Order & aSpaceInterpolationOrder,
                                                const mtk::Geometry_Type       & aTimeGeometryType,
                                                const Interpolation_Type       & aTimeInterpolationType,
                                                const mtk::Interpolation_Order & aTimeInterpolationOrder)
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
     } /* namespace fem */
} /* namespace moris */
