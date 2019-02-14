#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Factory.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Interpolation_Rule::Interpolation_Rule( const mtk::Geometry_Type       & aGeometryType,
                                                const Interpolation_Type       & aSpaceTimeInterpolationType,
                                                const mtk::Interpolation_Order & aSpaceTimeInterpolationOrder) : mGeometryType( aGeometryType ),
                                                                                                                 mSpaceInterpolationType( aSpaceTimeInterpolationType ),
                                                                                                                 mSpaceInterpolationOrder( aSpaceTimeInterpolationOrder ),
                                                                                                                 mTimeInterpolationType( Interpolation_Type::UNDEFINED ),
                                                                                                                 mTimeInterpolationOrder( mtk::Interpolation_Order::UNDEFINED ),
                                                                                                                 mSpaceTimeInterpolationType( aSpaceTimeInterpolationType ),
                                                                                                                 mSpaceTimeInterpolationOrder( aSpaceTimeInterpolationOrder ),
                                                                                                                 mHasTwoRulesFlag( false )
        {

        }

//------------------------------------------------------------------------------

        Interpolation_Rule::Interpolation_Rule(  const mtk::Geometry_Type       & aGeometryType,
                                                 const Interpolation_Type       & aSpaceInterpolationType,
                                                 const mtk::Interpolation_Order & aSpaceInterpolationOrder,
                                                 const Interpolation_Type       & aTimeInterpolationType,
                                                 const mtk::Interpolation_Order & aTimeInterpolationOrder) : mGeometryType( aGeometryType ),
                                                                                                             mSpaceInterpolationType( aSpaceInterpolationType ),
                                                                                                             mSpaceInterpolationOrder( aSpaceInterpolationOrder ),
                                                                                                             mTimeInterpolationType( aTimeInterpolationType ),
                                                                                                             mTimeInterpolationOrder( aTimeInterpolationOrder ),
                                                                                                             mSpaceTimeInterpolationType( Interpolation_Type::UNDEFINED ),
                                                                                                             mSpaceTimeInterpolationOrder( mtk::Interpolation_Order::UNDEFINED  ),
                                                                                                             mHasTwoRulesFlag( true )
        {

        }


//------------------------------------------------------------------------------

        Interpolation_Function_Base * Interpolation_Rule::create_space_time_interpolation_function() const
        {
            // create the factory
            Interpolation_Function_Factory tFactory;

            // return new interpolation function pointer
            return tFactory.create_interpolation_function( mGeometryType,
                                                           mSpaceTimeInterpolationType,
                                                           mSpaceTimeInterpolationOrder );
        }

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
            return tFactory.create_interpolation_function(
                    mtk::Geometry_Type::LINE,
                    mTimeInterpolationType,
                    mTimeInterpolationOrder );
        }

//------------------------------------------------------------------------------

        Interpolation_Type
        Interpolation_Rule::get_type_in_space() const
        {
            if( this->has_two_rules() )
            {
                return mSpaceInterpolationType;
            }
            else
            {
                return mSpaceTimeInterpolationType;
            }
        }

//------------------------------------------------------------------------------

        mtk::Interpolation_Order
        Interpolation_Rule::get_order_in_space() const
        {
            if( this->has_two_rules() )
            {
                return mSpaceInterpolationOrder;
            }
            else
            {
                return mSpaceTimeInterpolationOrder;
            }
        }

//----------------------------------------------------------------------------
     } /* namespace fem */
} /* namespace moris */
