#include "cl_FEM_Interpolation_Rule_Bis.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Factory.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        Interpolation_Rule_Bis::Interpolation_Rule_Bis(
                const mtk::Geometry_Type      		& aGeometryType,
                const Interpolation_Type      		& aSpaceInterpolationType,
                const mtk::Interpolation_Order     	& aSpaceInterpolationOrder) :
                    mGeometryType( aGeometryType ),
                    mSpaceInterpolationType( aSpaceInterpolationType),
                    mSpaceInterpolationOrder( aSpaceInterpolationOrder ),
                    mTimeInterpolationType( Interpolation_Type::UNDEFINED ),
                    mTimeInterpolationOrder( mtk::Interpolation_Order::UNDEFINED )
        {

        }

//------------------------------------------------------------------------------
        Interpolation_Rule_Bis::Interpolation_Rule_Bis(
                const mtk::Geometry_Type      	& aGeometryType,
                const Interpolation_Type      	& aSpaceInterpolationType,
                const mtk::Interpolation_Order  & aSpaceInterpolationOrder,
                const Interpolation_Type      	& aTimeInterpolationType,
                const mtk::Interpolation_Order  & aTimeInterpolationOrder) :
                    mGeometryType( aGeometryType ),
                    mSpaceInterpolationType( aSpaceInterpolationType ),
                    mSpaceInterpolationOrder( aSpaceInterpolationOrder ),
                    mTimeInterpolationType( aTimeInterpolationType ),
                    mTimeInterpolationOrder( aTimeInterpolationOrder )
        {

        }

//------------------------------------------------------------------------------
        Interpolation_Function_Base *
        Interpolation_Rule_Bis::create_space_interpolation_function() const
        {
            // create the factory
            Interpolation_Function_Factory tFactory;

            // return new interpolation function pointer
            return tFactory.create_interpolation_function(mGeometryType,
                    									  mSpaceInterpolationType,
														  mSpaceInterpolationOrder );
        }

//------------------------------------------------------------------------------
        Interpolation_Function_Base *
        Interpolation_Rule_Bis::create_time_interpolation_function() const
        {
            // create the factory
            Interpolation_Function_Factory tFactory;

            // return new interpolation function pointer
            return tFactory.create_interpolation_function(mtk::Geometry_Type::LINE,
                    									  mTimeInterpolationType,
														  mTimeInterpolationOrder );
        }

//------------------------------------------------------------------------------
        Interpolation_Function_Base *
        Interpolation_Rule_Bis::create_space_time_interpolation_function() const
        {
        	this->create_space_interpolation_function();
        	this->create_time_interpolation_function();


            //return both the space and the time interpolation function pointer
            MORIS_ERROR( false,"create_space_time_interpolation_function() not implemented " );
            return nullptr;
        }

//------------------------------------------------------------------------------
        Interpolation_Type
        Interpolation_Rule_Bis::get_space_interpolation_type() const
        {
        	return mSpaceInterpolationType;
        }

//------------------------------------------------------------------------------
        mtk::Interpolation_Order
        Interpolation_Rule_Bis::get_space_interpolation_order() const
        {
        	return mSpaceInterpolationOrder;
        }

//------------------------------------------------------------------------------
        Interpolation_Type
		Interpolation_Rule_Bis::get_time_interpolation_type() const
        {
        	return mTimeInterpolationType;
        }

//------------------------------------------------------------------------------
        mtk::Interpolation_Order
        Interpolation_Rule_Bis::get_time_interpolation_order() const
        {
        	return mTimeInterpolationOrder;
        }

//----------------------------------------------------------------------------
     } /* namespace fem */
} /* namespace moris */
